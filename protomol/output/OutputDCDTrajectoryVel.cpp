#include <protomol/output/OutputDCDTrajectoryVel.h>
#include <protomol/config/Configuration.h>
#include <protomol/output/OutputCache.h>
#include <protomol/module/MainModule.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/io/DCDTrajectoryWriter.h>
#include <protomol/base/Exception.h>
#include <protomol/ProtoMolApp.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____ OutputDCDTrajectoryVel
const string OutputDCDTrajectoryVel::keyword("DCDVELFile");

OutputDCDTrajectoryVel::OutputDCDTrajectoryVel() :
  Output(), myDCD(NULL), myMinimalImage(false) {}

OutputDCDTrajectoryVel::OutputDCDTrajectoryVel(const string &filename,
                                               int freq, bool minimal) :
  Output(freq), myDCD(new DCDTrajectoryWriter(filename)),
  myMinimalImage(minimal) {}

OutputDCDTrajectoryVel::~OutputDCDTrajectoryVel() {
  if (myDCD != NULL) delete myDCD;
}

void OutputDCDTrajectoryVel::doInitialize() {
  if (myDCD == NULL || !myDCD->open())
    THROW(string("Can not open '") + (myDCD ? myDCD->getFilename() : "") +
          "' for " + getId() + ".");
}

void OutputDCDTrajectoryVel::doRun(int) {
  const Vector3DBlock *vel = &app->velocities;
  if (!myDCD->write(*vel))
    THROW(string("Could not write ") + getId() + " '" + myDCD->getFilename()  +
          "'.");
}

void OutputDCDTrajectoryVel::doFinalize(int) {
  myDCD->close();
}

Output *OutputDCDTrajectoryVel::doMake(const vector<Value> &values) const {
  return new OutputDCDTrajectoryVel(values[0], values[1], values[2]);
}

void OutputDCDTrajectoryVel::getParameters(vector<Parameter> &parameter)
const {
  parameter.push_back
    (Parameter(getId(), Value(myDCD ? myDCD->getFilename() : "",
                              ConstraintValueType::NotEmpty())));
  parameter.push_back
    (Parameter(keyword + "OutputFreq",
               Value(getOutputFreq(), ConstraintValueType::Positive())));
  parameter.push_back
    (Parameter(keyword + "MinimalImage", Value(myMinimalImage),
               Text("whether the coordinates should be transformed to minimal "
                    "image or not (NA for Vel)")));
}

bool OutputDCDTrajectoryVel::adjustWithDefaultParameters(
  vector<Value> &values, const Configuration *config) const {
  if (!checkParameterTypes(values)) return false;

  if (config->valid(InputOutputfreq::keyword) && !values[1].valid())
    values[1] = (*config)[InputOutputfreq::keyword];

  if (config->valid(InputMinimalImage::keyword) && !values[2].valid())
    values[2] = (*config)[InputMinimalImage::keyword];

  return checkParameters(values);
}
