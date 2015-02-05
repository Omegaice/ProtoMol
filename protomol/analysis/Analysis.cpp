#include <protomol/analysis/Analysis.h>
#include <protomol/config/Configuration.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/module/MainModule.h>
#include <protomol/ProtoMolApp.h>

using namespace std;
using namespace ProtoMol;
using namespace ProtoMol::Report;

const string Analysis::scope("Analysis");

Analysis::Analysis(long freq) :
	app(0), firstStep(0), nextStep(0), lastStep(0), outputFreq(freq) {}

void Analysis::initialize(const ProtoMolApp *app) {
	this->app = app;

	if( outputFreq <= 0 ) {	//  used for finalize only output
		if( app->config.valid(InputOutputfreq::keyword)) {
			outputFreq = app->config[InputOutputfreq::keyword];
		} else {
			outputFreq = 1;
		}
	}

	if( app->config.valid(InputFirststep::keyword)) {
		nextStep = app->config[InputFirststep::keyword];
		firstStep = app->config[InputFirststep::keyword];
		lastStep = firstStep;
	}

	if( app->config.valid(InputNumsteps::keyword)) {
		lastStep = lastStep + app->config[InputNumsteps::keyword].operator long();
	}

	doInitialize();
}

bool Analysis::run(long step) {
	if( step >= nextStep ) {
		long n = ( step - nextStep ) / outputFreq;
		nextStep += max(n, 1L) * outputFreq;

		doRun(step);
		return true;
	}

	return false;
}

void Analysis::finalize(long step) {
	doRun(step);
	doFinalize(step);
}

Analysis *Analysis::make(const vector<Value> &values) const {
	assertParameters(values);
	return adjustAlias(doMake(values));
}

bool Analysis::isIdDefined(const Configuration *config) const {
	if( !addDoKeyword()) { return config->valid(getId()); }
	string str("do" + getId());
	if( !config->valid(getId()) || config->empty(str)) { return false; }
	if( !config->valid(str)) { return true; }
	return (bool)( *config )[str];
}

void Analysis::getParameters(vector<Parameter> &parameter) const {
	parameter.push_back(Parameter(getId() + "AnalysisFreq", Value(outputFreq, ConstraintValueType::Positive()), Text("output frequency")));
}

bool Analysis:: adjustWithDefaultParameters(std::vector<Value> &values, const Configuration *config) const {
	if( !checkParameterTypes(values)) { return false; }
	if( config->valid(InputOutputfreq::keyword) && !values[1].valid()) {
		values[1] = ( *config )[InputOutputfreq::keyword];
	}

	return checkParameters(values);
}
