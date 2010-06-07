#include <protomol/output/OutputXYZTrajectoryForce.h>
#include <protomol/module/MainModule.h>
#include <protomol/config/Configuration.h>
#include <protomol/integrator/Integrator.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/io/XYZTrajectoryWriter.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/base/Exception.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____ OutputXYZTrajectoryForce
const string OutputXYZTrajectoryForce::keyword("XYZForceFile");

OutputXYZTrajectoryForce::OutputXYZTrajectoryForce() :
  Output(), myXYZ() {}

OutputXYZTrajectoryForce::OutputXYZTrajectoryForce(const string &filename,
                                                   int freq) :
  Output(freq), myXYZ(new XYZTrajectoryWriter(filename)) {}

OutputXYZTrajectoryForce::~OutputXYZTrajectoryForce() {
  if (myXYZ != NULL) delete myXYZ;
}

void OutputXYZTrajectoryForce::doInitialize() {
  if (myXYZ == NULL || !myXYZ->open())
    THROW(string("Can not open '") + (myXYZ ? myXYZ->getFilename() : "") +
          "' for " + getId() + ".");
}

void OutputXYZTrajectoryForce::doRun(int) {
  if (!myXYZ->write(*(app->integrator->getForces()), app->topology->atoms,
                    app->topology->atomTypes))
    THROW(string("Could not write ") + getId() + " '" +
          myXYZ->getFilename() + "'.");
}

void OutputXYZTrajectoryForce::doFinalize(int) {
  myXYZ->close();
}

Output *OutputXYZTrajectoryForce::doMake(const vector<Value> &values) const {
  return new OutputXYZTrajectoryForce(values[0], values[1]);
}

void OutputXYZTrajectoryForce::getParameters(vector<Parameter> &parameter)
const {
  parameter.push_back
    (Parameter(getId(), Value(myXYZ ? myXYZ->getFilename() : "",
                              ConstraintValueType::NotEmpty())));
  parameter.push_back
    (Parameter(keyword + "OutputFreq",
               Value(getOutputFreq(), ConstraintValueType::Positive())));
}

bool OutputXYZTrajectoryForce::adjustWithDefaultParameters(
  vector<Value> &values, const Configuration *config) const {
  if (!checkParameterTypes(values)) return false;

  if (config->valid(InputOutputfreq::keyword) && !values[1].valid())
    values[1] = (*config)[InputOutputfreq::keyword];

  return checkParameters(values);
}
