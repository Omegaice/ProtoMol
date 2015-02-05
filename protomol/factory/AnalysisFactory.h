/* -*- c++ -*- */
#ifndef ANALYSIS_FACTORY_H
#define ANALYSIS_FACTORY_H

#include <protomol/base/Factory.h>

#include <protomol/analysis/Analysis.h>
#include <protomol/config/Value.h>

namespace ProtoMol {
	class AnalysisCollection;

	class AnalysisFactory : public Factory<Analysis> {
		public:
			virtual void registerHelpText () const;

		public:
			void registerAllExemplarsConfiguration (Configuration *config) const;
			Analysis *make (const std::string &id, const std::vector<Value> &values) const;
			AnalysisCollection *makeCollection (const Configuration *config) const;
	};
}

#endif	// ANALYSIS_FACTORY_H
