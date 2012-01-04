#ifndef OPENMMINTEGRATOR_H
#define OPENMMINTEGRATOR_H

#include <protomol/integrator/STSIntegrator.h>

#if defined (HAVE_OPENMM)
#include <OpenMM.h>
#include <LTMD/Integrator.h>
#include <LTMD/Parameters.h>
#endif

namespace ProtoMol {
	class ScalarStructure;
	class ForceGroup;
	class Vector3DBlock;

	class OpenMMIntegrator : public STSIntegrator {
		public:
			OpenMMIntegrator();
			OpenMMIntegrator( const std::vector<Value>& params, ForceGroup *overloadedForces );
			~OpenMMIntegrator();
		public:
			virtual std::string getIdNoAlias() const { return keyword; }
			virtual void getParameters( std::vector<Parameter> &parameters ) const;
			virtual unsigned int getParameterSize() const;

			virtual void initialize( ProtoMolApp *app );
			virtual void run( int numTimesteps );
		private:
			virtual STSIntegrator *doMake( const std::vector<Value> &values, ForceGroup *fg ) const;
			void getObcScaleFactors( std::vector<Real>& scaleFactors );
		public:
			static const std::string keyword;
		protected:
#if defined (HAVE_OPENMM)
			OpenMM::System *system;
			OpenMM::Integrator *integrator;
			OpenMM::Context *context;
#endif
			// OpenMM Parameters
			int mPlatform, mMinSteps;
			double mTolerance;
			
			// Integrator Parameters
			int mSeed;
			Real mTemperature, mGamma;
			
			// Solvent Parameters
			int mCommonMotionRate;
			Real mGBSAEpsilon, mGBSASolvent;
			
			// Force Switches
			bool isUsingHarmonicBondForce;
			bool isUsingHarmonicAngleForce;
			bool isUsingNonBondedForce;
			bool isUsingRBDihedralForce;
			bool isUsingPeriodicTorsionForce;
			bool isUsingGBForce;
			
			// LTMD Data
			bool isLTMD;
			OpenMM::LTMD::Parameters mLTMDParameters;
			std::vector<OpenMM::LTMD::Force> mForceList;
	};
}

#endif // OPENMMINTEGRATOR_H 
