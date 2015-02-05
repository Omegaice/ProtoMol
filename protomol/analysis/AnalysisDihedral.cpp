#include <protomol/analysis/AnalysisDihedral.h>

#include <protomol/topology/TopologyUtilities.h>
#include <protomol/module/MainModule.h>
#include <protomol/ProtoMolApp.h>

#include <protomol/base/MathUtilities.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/base/SystemUtilities.h>

#include <protomol/io/XYZWriter.h>
#include <protomol/io/CheckpointConfigWriter.h>

#include <sstream>
#include <iostream>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

const string AnalysisDihedral::keyword("Dihedral");

AnalysisDihedral::AnalysisDihedral(int freq) : Analysis(freq) {}

void AnalysisDihedral::doInitialize() {}

void AnalysisDihedral::doIt(long step) {}

void AnalysisDihedral::doRun(long step) {}

Analysis *AnalysisDihedral::doMake(const vector<Value> &values) const {
	return new AnalysisDihedral(toInt(values[1]));
}

bool AnalysisDihedral::isIdDefined(const Configuration *config) const {
	return config->valid(getId());
}

void AnalysisDihedral::getParameters(vector<Parameter> &parameter) const {
	parameter.push_back( Parameter( getId(), Value( true ), true ) );
	parameter.push_back(Parameter(getId() + "Freq", Value(outputFreq, ConstraintValueType::Positive()), Text("output frequency")));
}

bool AnalysisDihedral::adjustWithDefaultParameters(vector<Value> &values, const Configuration *config) const {
	if( !checkParameterTypes(values)) { return false; }

	if( config->valid(InputOutputfreq::keyword) && !values[1].valid()) {
		values[1] = ( *config )[InputOutputfreq::keyword];
	}

	return checkParameters(values);
}
