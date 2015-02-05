#ifndef PROTOMOL_ANALYSIS_COLLECTION_H
#define PROTOMOL_ANALYSIS_COLLECTION_H

#include <list>

namespace ProtoMol {
	class Analysis;
	class AnalysisFactory;
	class ProtoMolApp;

	// / Container class for Analysis objects invoked at application level.
	class AnalysisCollection {
		friend class AnalysisFactory;

		typedef std::list<Analysis *> Container;
		Container outputList;

		const ProtoMolApp *app;

		public:
			AnalysisCollection() : app(0) {}
			~AnalysisCollection();

			// / Initialize all Analysis object
			void initialize (const ProtoMolApp *app);

			// / Invoke all Analysis objects with run().  Returns true if an Analysis ran.
			bool run (long step);

			// / Finalize all Outout object
			void finalize (long step);

			// / Add new Analysis object to the collection
			void adoptAnalysis (Analysis *output);

			// / Iterators, const
			typedef Container::const_iterator const_iterator;
			const_iterator begin() const { return outputList.begin(); }
			const_iterator end()   const { return outputList.end(); }

		private:
			typedef Container::iterator iterator;
			iterator begin() { return outputList.begin(); }
			iterator end() { return outputList.end(); }
	};
}
#endif	// PROTOMOL_ANALYSIS_COLLECTION_H
