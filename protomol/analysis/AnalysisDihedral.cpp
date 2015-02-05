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

const string AnalysisDihedral::keyword("AnalyzeDihedral");

AnalysisDihedral::AnalysisDihedral() : Analysis(true) {}

void AnalysisDihedral::doInitialize() {}

void AnalysisDihedral::doIt(long step) {}

void AnalysisDihedral::doRun(long step) {
	const size_t dihedralCount = app->topology->dihedrals.size();
	std::cout << dihedralCount << std::endl;
	for( int i = 0; i < dihedralCount; i++ ) {
		if( i == 1 || i == 2 || i == 3 ) {
			const Real angle = app->outputCache.getDihedralPhi( i );
			std::cout << i << " " << angle << std::endl;
			if( angle >= -180.0 && angle <= 0.0 ) {
				std::cout << "Folded" << std::endl;
				//app->finalize();
			}
		}
	}
}

Analysis *AnalysisDihedral::doMake(const vector<Value> &values) const {
	return new AnalysisDihedral();
}

bool AnalysisDihedral::isIdDefined(const Configuration *config) const {
	return config->valid(getId());
}

void AnalysisDihedral::getParameters(vector<Parameter> &parameter) const {
	parameter.push_back( Parameter( getId(), Value( true ), true ) );
}
