#include <protomol/integrator/openMM/OpenMMIntegrator.h>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/base/PMConstants.h>
#include <protomol/base/Zap.h>
#include <protomol/topology/LennardJonesParameters.h>

#include <vector>
#include <algorithm>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;
//____ OpenMMIntegrator

const string OpenMMIntegrator::keyword( "OpenMM" );

OpenMMIntegrator::OpenMMIntegrator() :
	STSIntegrator() {
#if defined (HAVE_OPENMM)
	system = 0;
	bonds = 0;
	angles = 0;
	nonbonded = 0;
	integrator = 0;
	context = 0;
#endif
}

OpenMMIntegrator::
OpenMMIntegrator( Real timestep, ForceGroup *overloadedForces ) :
	STSIntegrator( timestep, overloadedForces ) {
#if defined (HAVE_OPENMM)
	system = 0;
	bonds = 0;
	angles = 0;
	nonbonded = 0;
	integrator = 0;
	context = 0;
#endif
}

OpenMMIntegrator::~OpenMMIntegrator() {
#if defined (HAVE_OPENMM)
	zap( context );
	zap( integrator );
	//zap(nonbonded);
	//zap(angles);
	zap( system );
	//zap(bonds);
#endif
}

struct NBForce {
	int atom1;
	int atom2;
	Real charge;
	Real sigma;
	Real epsilon;

	NBForce( int a, int b, Real c, Real s, Real e ) {
		atom1 = a;
		atom2 = b;
		charge = c;
		sigma = s;
		epsilon = e;
	}

	bool operator< ( const NBForce &other ) const {
		if( atom1 < other.atom1 ) {
			return true;
		}

		if( atom1 == other.atom1 ) {
			if( atom2 < other.atom2 ) {
				return true;
			}
		}

		return false;
	}
};

void OpenMMIntegrator::initialize( ProtoMolApp *app ) {
	STSIntegrator::initialize( app );
	initializeForces();

	//openMM

#if defined (HAVE_OPENMM)
	OpenMM::Platform::loadPluginsFromDirectory(
		OpenMM::Platform::getDefaultPluginsDirectory()
	);

	//find system size
	unsigned int sz = app->positions.size();

	//find constraint size
	const std::vector<Bond::Constraint> *myListOfConstraints =
		&( app->topology->bondRattleShakeConstraints );

	unsigned int numConstraints = ( *myListOfConstraints ).size();

	//Initialize system
	system = new OpenMM::System();
	for( unsigned int i = 0; i < sz; ++i ) {
		system->addParticle( app->topology->atoms[i].scaledMass );
	}

	//remove common motion?
	if( myCommonMotionRate > 0 ) {
		system->addForce( new OpenMM::CMMotionRemover( myCommonMotionRate ) );
	}

	//openMM forces
	if( HarmonicBondForce ) {
		unsigned int numBonds = app->topology->bonds.size();

		unsigned int numConstBonds = 0;

		if( numConstraints ) {
			for( unsigned int i = 0; i < numBonds; ++i ) {
				if( ( app->topology->atoms[ app->topology->bonds[i].atom1 ].name[0] == 'H' ) ||
						( app->topology->atoms[ app->topology->bonds[i].atom2 ].name[0] == 'H' ) ) {
					numConstBonds++;
				}
			}
		}

		bonds = new OpenMM::HarmonicBondForce();
		system->addForce( bonds );

		unsigned int bondsIndex = 0;

		for( unsigned int i = 0; i < numBonds; ++i ) {
			unsigned int a1 = app->topology->bonds[i].atom1; unsigned int a2 = app->topology->bonds[i].atom2;
			Real r_0 = app->topology->bonds[i].restLength  * Constant::ANGSTROM_NM;
			Real k = app->topology->bonds[i].springConstant
					 * Constant::KCAL_KJ * Constant::INV_ANGSTROM_NM * Constant::INV_ANGSTROM_NM * 2.0; //times 2 as Amber is 1/2 k(b-b_0)^2
			if( numConstraints ) {
				if( ( app->topology->atoms[ app->topology->bonds[i].atom1 ].name[0] != 'H' ) &&
						( app->topology->atoms[ app->topology->bonds[i].atom2 ].name[0] != 'H' ) ) {

					//bonds->setBondParameters(bondsIndex++, a1, a2, r_0, k);
				}
			} else {
				bonds->addBond( a1, a2, r_0, k );
			}
		}
	}

	if( HarmonicAngleForce ) {
		unsigned int numAngles = app->topology->angles.size();

		angles = new OpenMM::HarmonicAngleForce();
		system->addForce( angles );
		for( unsigned int i = 0; i < numAngles; i++ ) {
			unsigned int a1 = app->topology->angles[i].atom1;
			unsigned int a2 = app->topology->angles[i].atom2;
			unsigned int a3 = app->topology->angles[i].atom3;
			Real theta0 = acos( cos( app->topology->angles[i].restAngle ) );
			Real k_t = app->topology->angles[i].forceConstant * Constant::KCAL_KJ * 2.0; //times 2 as Amber is 1/2 k(a-a_0)^2
			//report << hint << "rest angle " << theta0 << " k " << k_t << " ang size " << numAngles << endr;

			angles->addAngle( a1, a2, a3, theta0, k_t );
		}
	}

	if( PeriodicTorsion ) {
		unsigned int numPTor = app->topology->dihedrals.size();
		unsigned int totalNumPTor = 0;

		//currently openMM cannot do mutilicity > 1
		//for (unsigned int i = 0; i < numPTor; i++)
		//  totalNumPTor += app->topology->dihedrals[i].multiplicity;

		PTorsion = new OpenMM::PeriodicTorsionForce();//totalNumPTor);
		system->addForce( PTorsion );
		for( unsigned int i = 0; i < numPTor; i++ ) {
			unsigned int a1 = app->topology->dihedrals[i].atom1;
			unsigned int a2 = app->topology->dihedrals[i].atom2;
			unsigned int a3 = app->topology->dihedrals[i].atom3;
			unsigned int a4 = app->topology->dihedrals[i].atom4;

			//unsigned int multiplicity = app->topology->dihedrals[i].multiplicity;

			for( unsigned int j = 0; j < 1/*multiplicity*/; j++ ) {

				unsigned int mult = app->topology->dihedrals[i].periodicity[j];
				Real phiA = app->topology->dihedrals[i].phaseShift[j];
				Real cpA = app->topology->dihedrals[i].forceConstant[j] * Constant::KCAL_KJ;

				//idef.iparams[type].pdihs.mult, idef.iparams[type].pdihs.phiA*M_PI/180.0, idef.iparams[type].pdihs.cpA

				PTorsion->addTorsion( a1, a2, a3, a4, mult, phiA, cpA );
			}
		}
	}

	if( RBDihedralForce ) {
		unsigned int numRBDih = app->topology->rb_dihedrals.size();

		RBDihedral = new OpenMM::RBTorsionForce();

		system->addForce( RBDihedral );
		for( unsigned int i = 0; i < numRBDih; i++ ) {
			unsigned int a1 = app->topology->rb_dihedrals[i].atom1;
			unsigned int a2 = app->topology->rb_dihedrals[i].atom2;
			unsigned int a3 = app->topology->rb_dihedrals[i].atom3;
			unsigned int a4 = app->topology->rb_dihedrals[i].atom4;
			Real C0 = app->topology->rb_dihedrals[i].C0 * Constant::KCAL_KJ;
			Real C1 = app->topology->rb_dihedrals[i].C1 * Constant::KCAL_KJ;
			Real C2 = app->topology->rb_dihedrals[i].C2 * Constant::KCAL_KJ;
			Real C3 = app->topology->rb_dihedrals[i].C3 * Constant::KCAL_KJ;
			Real C4 = app->topology->rb_dihedrals[i].C4 * Constant::KCAL_KJ;
			Real C5 = app->topology->rb_dihedrals[i].C5 * Constant::KCAL_KJ;
			if( C0 != 0 || C1 != 0 || C2 != 0 || C3 != 0 || C4 != 0 || C5 != 0 ) {
				RBDihedral->addTorsion( a1, a2, a3, a4, C0, C1, C2, C3, C4, C5 );
			}
		}
	}

	if( NonbondedForce ) {

		//get 1-4 interaction size
		unsigned int exclSz = app->topology->exclusions.getTable().size();

		nonbonded = new OpenMM::NonbondedForce();//0);
		system->addForce( nonbonded );

		//normal interactions
		for( unsigned int i = 0; i < sz; i++ ) {
			int type1 = app->topology->atoms[i].type;
			Real sigma = app->topology->atomTypes[type1].sigma;
			//topo->atomTypes[i].sigma14 = par.nonbondeds[bi].sigma14;
			Real epsilon = app->topology->atomTypes[type1].epsilon;
			//topo->atomTypes[i].epsilon14 = par.nonbondeds[bi].epsilon14;
			Real charge = app->topology->atoms[i].scaledCharge / Constant::SQRTCOULOMBCONSTANT;
			Real mass = app->topology->atoms[i].scaledMass;

			nonbonded->addParticle( charge, sigma * Constant::ANGSTROM_NM , epsilon * Constant::KCAL_KJ );

		}

		//1-4 interactions, note fudgeQQ set to 0.6059 and fudgeLJ set to 0.33333.
		for( unsigned int i = 0; i < exclSz; i++ ) {

			unsigned int atom1 = ( app->topology->exclusions.getTable() )[i].a1;
			unsigned int atom2 = ( app->topology->exclusions.getTable() )[i].a2;

			//full exclusion?
			if( ( app->topology->exclusions.getTable() )[i].excl == EXCLUSION_FULL ) {
				nonbonded->addException( atom1, atom2, 0.0, 1.0, 0.0 );
			}

			//modified exclusion?
			if( ( app->topology->exclusions.getTable() )[i].excl == EXCLUSION_MODIFIED ) {

				unsigned int type1 = app->topology->atoms[atom1].type;
				unsigned int type2 = app->topology->atoms[atom2].type;
				Real sigma = 0.5 * ( app->topology->atomTypes[type1].sigma +
									 app->topology->atomTypes[type2].sigma );
				Real epsilon = app->topology->LJScalingFactor *
							   sqrt( app->topology->atomTypes[type1].epsilon *
									 app->topology->atomTypes[type2].epsilon );
				Real chargeij =  app->topology->coulombScalingFactor * //FudgeQQ
								 ( app->topology->atoms[atom1].scaledCharge / Constant::SQRTCOULOMBCONSTANT ) *
								 ( app->topology->atoms[atom2].scaledCharge / Constant::SQRTCOULOMBCONSTANT );

				nonbonded->addException( atom1, atom2, chargeij, sigma * Constant::ANGSTROM_NM, epsilon * Constant::KCAL_KJ );

			}
		}
	}

	// Add GBSA if needed.
	if( app->topology->implicitSolvent  == GBSA && GBForce ) {
		gbsa = new OpenMM::GBSAOBCForce();
		system->addForce( gbsa );

		gbsa->setSoluteDielectric( myGBSAEpsilon ); //ir->epsilon_r);
		gbsa->setSolventDielectric( myGBSASolvent ); //ir->gb_epsilon_solvent);

		vector<Real> scaleFactors;
		getObcScaleFactors( scaleFactors );

		for( unsigned int i = 0; i < sz; ++i ) {
			Real charge = app->topology->atoms[i].scaledCharge / Constant::SQRTCOULOMBCONSTANT;
			unsigned int type = app->topology->atoms[i].type;
			Real radius = app->topology->atoms[i].myGBSA_T->vanDerWaalRadius * Constant::ANGSTROM_NM; //0.1 factor in openMM, file in A

			gbsa->addParticle( charge, radius, scaleFactors[i] );

		}
	}

	// Set constraints.
	for( unsigned int i = 0; i < numConstraints; ++i ) {

		int atom1 = ( *myListOfConstraints )[i].atom1;
		int atom2 = ( *myListOfConstraints )[i].atom2;
		Real restLength = ( *myListOfConstraints )[i].restLength * Constant::ANGSTROM_NM;

		system->addConstraint( atom1, atom2, restLength );
	}

#ifndef HAVE_OPENMM_OLD
	integrator = new OpenMM::LangevinIntegrator( myLangevinTemperature, myGamma, getTimestep() * Constant::FS_PS );
#else
	if( myIntegratorType == 1 ) {
		integrator = new OpenMM::LangevinIntegrator( myLangevinTemperature, myGamma, getTimestep() * Constant::FS_PS );
	} else {
		integrator = new OpenMM::NMLIntegrator( myLangevinTemperature, myGamma, getTimestep() * Constant::FS_PS, &app->eigenInfo );
	}
#endif

	context = new OpenMM::Context( *system, *integrator );

	OpenMM::Vec3 openMMvecp, openMMvecv;
	for( unsigned int i = 0; i < sz; ++i ) {
		for( int j = 0; j < 3; j++ ) {
			openMMvecp[j] = app->positions[i].c[j] * Constant::ANGSTROM_NM;
			openMMvecv[j] = app->velocities[i].c[j] * Constant::ANGSTROM_NM
							* Constant::INV_TIMEFACTOR * Constant::PS_FS;
		}
		openMMpositions.push_back( openMMvecp );
		openMMvelocities.push_back( openMMvecv );
	}

	context->setPositions( openMMpositions );
	context->setVelocities( openMMvelocities );

	//print platform
	report << plain << "OpenMM platform is: '" << context->getPlatform().getName()
		   << "' Integrator " << myIntegratorType << "." << endr;
	const OpenMM::State state = context->getState( OpenMM::State::Positions |
								OpenMM::State::Velocities |
								OpenMM::State::Forces |
								OpenMM::State::Energy );

#else

	//print platform
	report << plain << "OpenMM platform is not available." << endr;

#endif

}

void OpenMMIntegrator::run( int numTimesteps ) {

	preStepModify();

#if defined (HAVE_OPENMM)
	unsigned int sz = app->positions.size();

	// do integration
	integrator->step( numTimesteps );

	// Retrive data
	const OpenMM::State state = context->getState( OpenMM::State::Positions |
								OpenMM::State::Velocities |
								OpenMM::State::Forces |
								OpenMM::State::Energy );
	openMMpositions = state.getPositions();
	openMMvelocities = state.getVelocities();
	openMMforces = state.getForces();
	report.precision( 5 );
	//report << plain << "x1 " << openMMpositions[0][0] << ",y1 " << openMMpositions[0][1] << ",z1 " << openMMpositions[0][2] << endr;
	//report << plain << "vx1 " << openMMvelocities[0][0] << ",vy1 " << openMMvelocities[0][1] << ",vz1 " << openMMvelocities[0][2] << endr;

	for( unsigned int i = 0; i < sz; ++i ) {
		for( int j = 0; j < 3; j++ ) {
			app->positions[i].c[j] = openMMpositions[i][j] * Constant::NM_ANGSTROM; //nm to A
			app->velocities[i].c[j] = openMMvelocities[i][j] * Constant::NM_ANGSTROM *
									  Constant::TIMEFACTOR * Constant::FS_PS; //nm/ps to A/fs?
			( *myForces )[i].c[j] = openMMforces[i][j] * Constant::INV_NM_ANGSTROM * Constant::KJ_KCAL; //KJ/nm to Kcal/A
		}
	}

	app->energies.clear();
	//calculateForces();

	//std::cout << "Forces " << (*myForces)[0].c[0] << std::endl;
	//clear old energies
	app->energies[ScalarStructure::COULOMB] =
		app->energies[ScalarStructure::LENNARDJONES] =
			app->energies[ScalarStructure::BOND] =
				app->energies[ScalarStructure::ANGLE] =
					app->energies[ScalarStructure::DIHEDRAL] =
						app->energies[ScalarStructure::IMPROPER] = 0.0;

	//save total potential energy
	app->energies[ScalarStructure::OTHER] = state.getPotentialEnergy() * Constant::KJ_KCAL;

	//std::cout << "Temp " << temperature(app->topology, &app->velocities) << std::endl;

	//state.getKineticEnergy();

	//fix time as no forces calculated
	app->topology->time += numTimesteps * getTimestep();
#endif

	postStepModify();
}

void OpenMMIntegrator::calcForces() {
#if defined (HAVE_OPENMM)
	const OpenMM::State state = context->getState( OpenMM::State::Forces );
#endif
}

void OpenMMIntegrator::getParameters( vector<Parameter> &parameters )
const {
	STSIntegrator::getParameters( parameters );
	parameters.push_back( Parameter( "temperature", Value( myLangevinTemperature, ConstraintValueType::NotNegative() ) ) );
	parameters.push_back( Parameter( "gamma", Value( myGamma, ConstraintValueType::NotNegative() ) ) ); // * (1000 * Constant::INV_TIMEFACTOR)
	parameters.push_back( Parameter( "seed", Value( mySeed, ConstraintValueType::NotNegative() ), 1234 ) );
	//OpenMM forces
	parameters.push_back( Parameter( "HarmonicBondForce", Value( HarmonicBondForce, ConstraintValueType::NoConstraints() ), false ) );
	parameters.push_back( Parameter( "HarmonicAngleForce", Value( HarmonicAngleForce, ConstraintValueType::NoConstraints() ), false ) );
	parameters.push_back( Parameter( "RBDihedralForce", Value( RBDihedralForce, ConstraintValueType::NoConstraints() ), false ) );
	parameters.push_back( Parameter( "PeriodicTorsion", Value( PeriodicTorsion, ConstraintValueType::NoConstraints() ), false ) );
	parameters.push_back( Parameter( "NonbondedForce", Value( NonbondedForce, ConstraintValueType::NoConstraints() ), false ) );
	parameters.push_back( Parameter( "IntegratorType", Value( myIntegratorType, ConstraintValueType::NotNegative() ), 1 ) );
	//Implicit solvent parameters
	parameters.push_back( Parameter( "GBSAEpsilon", Value( myGBSAEpsilon, ConstraintValueType::NotNegative() ), 1.0 ) );
	parameters.push_back( Parameter( "GBSASolvent", Value( myGBSASolvent, ConstraintValueType::NotNegative() ), 78.3 ) );
	parameters.push_back( Parameter( "commonmotion", Value( myCommonMotionRate, ConstraintValueType::NotNegative() ), 0.0 ) );
	parameters.push_back( Parameter( "GBForce", Value( GBForce, ConstraintValueType::NoConstraints() ), true ) );
}

STSIntegrator *OpenMMIntegrator::doMake( const vector<Value> &values,
		ForceGroup *fg ) const {
	OpenMMIntegrator *myIntegP = new OpenMMIntegrator( values[0], fg );

	std::vector<Value> myValues( values.begin() + 1, values.end() );

	myIntegP->setupValues( myValues );

	return ( STSIntegrator * )myIntegP;

}

void OpenMMIntegrator::setupValues( std::vector<Value> &values ) {

	//these must be in the same order as getParameters()
	myLangevinTemperature = values[0];
	myGamma = ( Real )values[1]; // / (1000.0 * Constant::INV_TIMEFACTOR);
	// gamma is in Kcal/ps, myGamma is in Kcal/(fs*INV_TIMEFACTOR)
	mySeed = values[2];
	HarmonicBondForce = values[3];
	HarmonicAngleForce = values[4];
	RBDihedralForce = values[5];
	PeriodicTorsion = values[6];
	NonbondedForce = values[7];
	myIntegratorType = values[8];
	myGBSAEpsilon = values[9];
	myGBSASolvent = values[10];
	myCommonMotionRate = values[11];
	GBForce = values[12];
}

/**
 * Figure out OBC scale factors based on the atomic masses.
 */

void OpenMMIntegrator::getObcScaleFactors( vector<Real>& scaleFactors ) {

	unsigned int numAtoms = app->positions.size();

	scaleFactors.resize( numAtoms );
	for( unsigned int atomI = 0; atomI < numAtoms; atomI++ ) {

		Real scaleFactor = 0.8;
		Real mass        = app->topology->atoms[atomI].scaledMass;

		if( mass < 1.2 && mass >= 1.0 ) {        // hydrogen
			scaleFactor  = 0.85;
		} else if( mass > 11.8 && mass < 12.2 ) { // carbon
			scaleFactor  = 0.72;
		} else if( mass > 14.0 && mass < 15.0 ) { // nitrogen
			scaleFactor  = 0.79;
		} else if( mass > 15.5 && mass < 16.5 ) { // oxygen
			scaleFactor  = 0.85;
		} else if( mass > 31.5 && mass < 32.5 ) { // sulphur
			scaleFactor  = 0.96;
		} else if( mass > 29.5 && mass < 30.5 ) { // phosphorus
			scaleFactor  = 0.86;
		} else {
			report << plain << " Warning: mass for atom with mass = " << mass << " not recognized." << endr;
		}

		scaleFactors[atomI] = scaleFactor;
	}
}


