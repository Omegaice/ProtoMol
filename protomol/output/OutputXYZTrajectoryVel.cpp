#include <protomol/output/OutputXYZTrajectoryVel.h>
#include <protomol/config/Configuration.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/module/MainModule.h>
#include <protomol/base/Exception.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/io/XYZTrajectoryWriter.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____ OutputXYZTrajectoryVel
const string OutputXYZTrajectoryVel::keyword("XYZVelFile");

OutputXYZTrajectoryVel::OutputXYZTrajectoryVel() :
  Output(), myXYZ() {}

OutputXYZTrajectoryVel::OutputXYZTrajectoryVel(const string &filename,
                                               int freq) :
  Output(freq), myXYZ(new XYZTrajectoryWriter(filename)) {}

OutputXYZTrajectoryVel::~OutputXYZTrajectoryVel() {
  if (myXYZ != NULL) delete myXYZ;
}

void OutputXYZTrajectoryVel::doInitialize() {
  if (myXYZ == NULL || !myXYZ->open())
    THROW(string(" Can not open '") + (myXYZ ? myXYZ->getFilename() : "") +
                 "' for " + getId() + ".");
}

void OutputXYZTrajectoryVel::doRun(int) {
  if (!myXYZ->write(*&app->velocities, app->topology->atoms,
                    app->topology->atomTypes))
    THROW(string("Could not write ") + getId() + " '" +
          myXYZ->getFilename() + "'.");
}

void OutputXYZTrajectoryVel::doFinalize(int) {
  myXYZ->close();
}

Output *OutputXYZTrajectoryVel::doMake(const vector<Value> &values) const {
  return new OutputXYZTrajectoryVel(values[0], values[1]);
}

void OutputXYZTrajectoryVel::getParameters(vector<Parameter> &parameter)
const {
  parameter.push_back
    (Parameter(getId(), Value(myXYZ != NULL ? myXYZ->getFilename() : "",
                              ConstraintValueType::NotEmpty())));

  parameter.push_back
    (Parameter(keyword + "OutputFreq",
               Value(getOutputFreq(), ConstraintValueType::Positive())));
}

bool OutputXYZTrajectoryVel::adjustWithDefaultParameters(
  vector<Value> &values, const Configuration *config) const {

  if (!checkParameterTypes(values)) return false;

  if (config->valid(InputOutputfreq::keyword) && !values[1].valid())
    values[1] = (*config)[InputOutputfreq::keyword];

  return checkParameters(values);
}
