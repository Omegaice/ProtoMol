#ifndef OPENMMFBM_H
#define OPENMMFBM_H

#include <protomol/integrator/openMM/OpenMMIntegrator.h>

namespace ProtoMol {
	class ScalarStructure;
	class ForceGroup;

	class OpenMMFBMIntegrator : public OpenMMIntegrator {
		public:
			OpenMMFBMIntegrator();
			OpenMMFBMIntegrator( const std::vector<Value>& base, const std::vector<Value>& params, ForceGroup *forces );
			~OpenMMFBMIntegrator();
			
			virtual std::string getIdNoAlias() const { return keyword; }
			virtual void getParameters( std::vector<Parameter>& parameters ) const;
			virtual unsigned int getParameterSize() const;
			
			virtual void initialize( ProtoMolApp *app );
			virtual long run( const long numTimesteps );
		private:
			virtual STSIntegrator *doMake( const std::vector<Value>& values, ForceGroup *fg )const;
		public:
			static const std::string keyword;
		private:
			int mResiduesPerBlock;
			int  mBlockDOF;
			Real mBlockDelta;
			Real mSDelta;
			int mModes;
			int mBlockPlatform;
			std::string mEigenvalueFilename;
			std::string mEigenvectorFilename;

	};
}

#endif
