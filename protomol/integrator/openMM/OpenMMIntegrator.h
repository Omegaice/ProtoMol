#ifndef OPENMMINTEGRATOR_H
#define OPENMMINTEGRATOR_H

#include <protomol/integrator/STSIntegrator.h>

#if defined (HAVE_OPENMM)
#include <OpenMM.h>
#include "LTMD/Integrator.h"
#endif

namespace ProtoMol {
	class ScalarStructure;
	class ForceGroup;
	class Vector3DBlock;

	class OpenMMIntegrator : public STSIntegrator {

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Constructors, destructors, assignment
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		public:
			OpenMMIntegrator();
			OpenMMIntegrator( Real timestep, ForceGroup *overloadedForces );
			~OpenMMIntegrator();

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// From class Makeable
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		public:
			virtual std::string getIdNoAlias() const {
				return keyword;
			}
			virtual void getParameters( std::vector<Parameter> &parameters ) const;

			virtual unsigned int getParameterSize() const {
				return 13;
			}

		protected:
			virtual void setupValues( std::vector<Value> &params );

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// From class Integrator
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		public:
			virtual void initialize( ProtoMolApp *app );
			virtual void run( int numTimesteps );
			virtual void calcForces();
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// From class StandardIntegrator
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		protected:

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// From class STSIntegrator
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		private:
			virtual STSIntegrator *doMake( const std::vector<Value> &values,
										   ForceGroup *fg ) const;

		protected:

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// New methods of class OpenMMIntegrator
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		private:
			void getObcScaleFactors( vector<Real>& scaleFactors );

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// My data members
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		public:
			static const std::string keyword;
		protected:
			Real myLangevinTemperature;
			Real myGamma;
			int mySeed;
			
#if defined (HAVE_OPENMM)
			OpenMM::System *system;
			OpenMM::HarmonicBondForce *bonds;
			OpenMM::HarmonicAngleForce *angles;
			OpenMM::RBTorsionForce *RBDihedral;
			OpenMM::PeriodicTorsionForce *PTorsion;
			OpenMM::NonbondedForce *nonbonded;
			OpenMM::GBSAOBCForce *gbsa;
			OpenMM::Integrator *integrator;
			OpenMM::Context *context;
			vector<OpenMM::Vec3> openMMpositions, openMMvelocities, openMMforces;
#endif

			//OpenMM force/integrator parameters
			bool HarmonicBondForce, HarmonicAngleForce, NonbondedForce, RBDihedralForce, PeriodicTorsion, GBForce;
			int  myIntegratorType;
			Real myGBSAEpsilon, myGBSASolvent;
			int myCommonMotionRate;
			int resPerBlock;
			int bdof;
			double delta;
			int modes;
			int rediagFreq;
			double minLimit;
	};
}

#endif // OPENMMINTEGRATOR_H 
