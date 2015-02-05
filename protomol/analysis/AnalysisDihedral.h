#ifndef PROTOMOL_ANALYSIS_DIHEDRAL_H
#define PROTOMOL_ANALYSIS_DIHEDRAL_H

#include <protomol/base/StringUtilities.h>
#include <protomol/analysis/Analysis.h>

namespace ProtoMol {
	class Configuration;

	class AnalysisDihedral : public Analysis {
		public:
			static const std::string keyword;

		public:
			AnalysisDihedral();

		public:
			void doIt (long step);

		private:
			Analysis *doMake (const std::vector<Value> &values) const;
			void doInitialize ();
			void doRun (long step);
			void doFinalize(long) {}
			bool isIdDefined (const Configuration *config) const;
			bool addDoKeyword() const { return false; }

		public:
			std::string getIdNoAlias() const { return keyword; }
			void getParameters (std::vector<Parameter> &) const;
	};
}

#endif	// PROTOMOL_ANALYSIS_DIHEDRAL_H
