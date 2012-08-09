#include <protomol/integrator/openMM/OpenMMFBMIntegrator.h>
#include <protomol/base/Report.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/ProtoMolApp.h>

#ifdef HAVE_OPENMM_FBM
#include <OpenMMFBM/FBMParameters.h>
#include <OpenMMFBM/FlexibleBlockMethod.h>
#endif

#include <iostream>

using namespace ProtoMol::Report;

namespace ProtoMol {
	const string OpenMMFBMIntegrator::keyword( "OpenMMFBM" );

	OpenMMFBMIntegrator::OpenMMFBMIntegrator() : OpenMMIntegrator() {
	
	}

	OpenMMFBMIntegrator::OpenMMFBMIntegrator( const std::vector<Value>& base, const std::vector<Value>& params, ForceGroup *forces )
		: OpenMMIntegrator( base, forces ) {
		isLTMD = false;
		
		mModes = params[0];
		mResiduesPerBlock = params[1];
		mBlockDOF = params[2];
		mBlockDelta = params[3];
		mBlockPlatform = params[4];
		mSDelta = params[5];
		mEigenvalueFilename = (std::string) params[6];
		mEigenvectorFilename = (std::string) params[7];
		
	}

	OpenMMFBMIntegrator::~OpenMMFBMIntegrator() {

	}

	void OpenMMFBMIntegrator::initialize( ProtoMolApp *app ) {
		#ifdef HAVE_OPENMM_FBM
		// Setup LTMD Parameters
		mFBMParameters.blockDelta = mBlockDelta * Constant::ANGSTROM_NM;
		mFBMParameters.sDelta = mSDelta * Constant::ANGSTROM_NM;
		mFBMParameters.bdof = mBlockDOF;
		mFBMParameters.res_per_block = mResiduesPerBlock;
		mFBMParameters.modes = mModes;
		switch( mBlockPlatform ){
	   			case 0: // Reference
					std::cout << "OpenMM Block Diagonalization Platform: Reference" << std::endl;
					mFBMParameters.blockPlatform = OpenMMFBM::Preference::Reference;
					break;
				case 1:	// OpenCL
					std::cout << "OpenMM Block Diagonalization Platform: OpenCL" << std::endl;
					mFBMParameters.blockPlatform = OpenMMFBM::Preference::OpenCL;
					break;
				case 2: // CUDA
					std::cout << "OpenMM Block Diagonalization Platform: CUDA" << std::endl;
					mFBMParameters.blockPlatform = OpenMMFBM::Preference::CUDA;
					break;
		}
		
		int current_res = app->topology->atoms[0].residue_seq;
		int res_size = 0;
		for( int i = 0; i < app->topology->atoms.size(); i++ ) {
			if( app->topology->atoms[i].residue_seq != current_res ) {
				mFBMParameters.residue_sizes.push_back( res_size );
				current_res = app->topology->atoms[i].residue_seq;
				res_size = 0;
			}
			res_size++;
		}
		mFBMParameters.residue_sizes.push_back( res_size );
		#endif
		
		//initialize base
		OpenMMIntegrator::initialize( app );

		initializeForces();
	}

  long OpenMMFBMIntegrator::run( const long numTimesteps ) {
		if( numTimesteps < 1 ) return 0;
		
		preStepModify();

		#ifdef HAVE_OPENMM
		#ifdef HAVE_OPENMM_FBM
		cout << "Diagonalizing" << endl;
		OpenMMFBM::FlexibleBlockMethod fbm(mFBMParameters);
		fbm.diagonalize(*context);

		vector<vector<double> > blockHessian;
		fbm.getBlockHessian(blockHessian);

		fstream blockHessianOut;
		blockHessianOut.open("blockHessian.txt", fstream::out);
		blockHessianOut.precision(10);
		for(unsigned int col = 0; col < blockHessian[0].size(); col++)
		  {
		    for(unsigned int row = 0; row < blockHessian.size(); row++)
		      {
			if(blockHessian[row][col] != 0.0)
			  {
			    blockHessianOut << (row + 1) << " " << (col + 1) << " " << blockHessian[row][col] << endl;
			  }
		      }
		  }
		blockHessianOut.close();

		vector<vector<double> > blockEigenvectors;
		fbm.getBlockEigenvectors(blockEigenvectors);

		fstream blockEigenvectorsOut;
		blockEigenvectorsOut.open("block_eigenvectors.txt", fstream::out);
		blockEigenvectorsOut.precision(10);
		for(unsigned int col = 0; col < blockEigenvectors[0].size(); col++)
		  {
		    for(unsigned int row = 0; row < blockEigenvectors.size(); row++)
		      {
			if(blockEigenvectors[row][col] != 0.0)
			  {
			    blockEigenvectorsOut << (row + 1) << " " << (col + 1) << " " << blockEigenvectors[row][col] << endl;
			  }
		      }
		  }
		blockEigenvectorsOut.close();

		
		cout << "Retrieving eigenpairs" << endl;

		const vector<vector<OpenMM::Vec3> > eigenvectors = fbm.getEigenvectors();
		const vector<Real> eigenvalues = fbm.getEigenvalues();

		cout << "Writing out eigenpairs" << endl;
		
		fstream eigenvalues_out;
		eigenvalues_out.open(mEigenvalueFilename.c_str(), fstream::out);
		eigenvalues_out << eigenvalues.size() << endl;
		eigenvalues_out << "! Eigenvalues from OpenMM FBM" << endl;
		for(unsigned int i = 0; i < eigenvalues.size(); i++)
		  {
		    eigenvalues_out << (i + 1) << " " << eigenvalues[i] << endl;
		  }
		eigenvalues_out.close();

		cout << "Wrote out eigenvalues" << endl;

		fstream eigenvectors_out;
		eigenvectors_out.open(mEigenvectorFilename.c_str(), fstream::out);
		eigenvectors_out << eigenvectors[0].size() << " " << eigenvectors.size() << " 0.0" << endl;
		eigenvectors_out << "! Eigenvectors from OpenMM FBM" << endl;
		for(unsigned int i = 0; i < eigenvectors.size(); i++)
		  {
		    for(unsigned int j = 0; j < eigenvectors[0].size(); j++)
		      {
			eigenvectors_out << (j + 1) << " " << (i + 1) << " ";
			for(unsigned int k = 0; k < 3; k++)
			  {
			    eigenvectors_out << eigenvectors[i][j][k] << " ";
			  }
			eigenvectors_out << endl;
		      }
		  }
		eigenvectors_out.close();

		cout << "Wrote out eigenpairs" << endl;

		const OpenMM::State state = context->getState(OpenMM::State::Positions |
							      OpenMM::State::Velocities |
							      OpenMM::State::Forces |
							      OpenMM::State::Energy);

		const std::vector<OpenMM::Vec3> positions(state.getPositions());
		const std::vector<OpenMM::Vec3> velocities(state.getVelocities());
		const std::vector<OpenMM::Vec3> forces(state.getForces());
		
		const unsigned int sz = app->positions.size();
		for(unsigned int i = 0; i < sz; i++)
		{
			for(unsigned int j = 0; j < 3; j++)
			{
			  app->positions[i].c[j] = positions[i][j] * Constant::NM_ANGSTROM;
			  app->velocities[i].c[j] = velocities[i][j] * Constant::NM_ANGSTROM * Constant::TIMEFACTOR * Constant::FS_PS;
			  (*myForces)[i].c[j] = forces[i][j] * Constant::INV_NM_ANGSTROM * Constant::KJ_KCAL;
			}
		}

		app->energies.clear();

		app->energies[ScalarStructure::COULOMB] = 
		  app->energies[ScalarStructure::LENNARDJONES] =
		  app->energies[ScalarStructure::BOND] =
		  app->energies[ScalarStructure::ANGLE] =
		  app->energies[ScalarStructure::DIHEDRAL] =
		  app->energies[ScalarStructure::IMPROPER] = 0.0;

		app->energies[ScalarStructure::OTHER] = state.getPotentialEnergy() * Constant::KJ_KCAL;
		app->topology->time += numTimesteps * getTimestep();

		#endif
		#endif
		
		postStepModify();
    
		return numTimesteps;
	}

	void OpenMMFBMIntegrator::getParameters( vector<Parameter>& parameters ) const {
		OpenMMIntegrator::getParameters( parameters );

		parameters.push_back( Parameter( "numbermodes", Value( mModes, ConstraintValueType::NoConstraints() ), 1, Text( "Number of modes propagated" ) ) );
		parameters.push_back( Parameter( "resPerBlock", Value( mResiduesPerBlock, ConstraintValueType::NotNegative() ), 1 ) );
		parameters.push_back( Parameter( "bdof", Value( mBlockDOF, ConstraintValueType::NotNegative() ), 12 ) );
		parameters.push_back( Parameter( "blockEpsilon", Value( mBlockDelta, ConstraintValueType::NotNegative() ), 1e-3 ) );
		parameters.push_back( Parameter( "blockHessianPlatform", Value( mBlockPlatform, ConstraintValueType::NoConstraints() ), 0 ) );
		parameters.push_back( Parameter( "sEpsilon", Value( mSDelta, ConstraintValueType::NotNegative() ), 1e-3 ) );
		parameters.push_back( Parameter( "eigvalFile", Value ( mEigenvalueFilename, ConstraintValueType::NoConstraints() ), string("") ) ); 
		parameters.push_back( Parameter( "eigvecFile", Value ( mEigenvectorFilename, ConstraintValueType::NoConstraints() ), string("") ) ); 

	}

	STSIntegrator *OpenMMFBMIntegrator::doMake( const vector<Value>& values, ForceGroup *fg ) const {
		const std::vector<Value> base( values.begin(), values.begin() + OpenMMIntegrator::getParameterSize() );
		const std::vector<Value> params( values.begin() + OpenMMIntegrator::getParameterSize(), values.end() );
		
		return ( STSIntegrator * ) new OpenMMFBMIntegrator( base, params, fg );
	}

	unsigned int OpenMMFBMIntegrator::getParameterSize() const {
	  return OpenMMFBMIntegrator::getParameterSize() + 8;
	}
}

