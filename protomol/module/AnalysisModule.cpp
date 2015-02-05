#include <protomol/module/AnalysisModule.h>

#include <protomol/ProtoMolApp.h>
#include <protomol/config/Configuration.h>
#include <protomol/factory/AnalysisFactory.h>

#include <protomol/analysis/AnalysisDihedral.h>

using namespace ProtoMol;

defineInputValue(InputAnalysis,"analysis");
defineInputValue(InputAnalysisfreq,"analysisfreq");

void AnalysisModule::init(ProtoMolApp *app) {
	InputAnalysis::registerConfiguration(&app->config, true);
	InputAnalysisfreq::registerConfiguration(&app->config, 1L);

	app->analysisFactory.registerExemplar(new AnalysisDihedral());
}
