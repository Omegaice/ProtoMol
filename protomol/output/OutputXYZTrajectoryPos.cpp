#include <protomol/output/OutputXYZTrajectoryPos.h>
#include <protomol/module/MainModule.h>
#include <protomol/output/OutputCache.h>
#include <protomol/config/Configuration.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/base/Exception.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/io/XYZTrajectoryWriter.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____ OutputXYZTrajectoryPos
const string OutputXYZTrajectoryPos::keyword("XYZPosFile");

OutputXYZTrajectoryPos::OutputXYZTrajectoryPos() :
  Output(), myXYZ(NULL), myMinimalImage(false) {}

OutputXYZTrajectoryPos::OutputXYZTrajectoryPos(const string &filename, int freq,
                                               bool minimal) :
  Output(freq), myXYZ(new XYZTrajectoryWriter(filename)),
  myMinimalImage(minimal) {}

OutputXYZTrajectoryPos::~OutputXYZTrajectoryPos() {
  if (myXYZ != NULL) delete myXYZ;
}

void OutputXYZTrajectoryPos::doInitialize() {
  if (myXYZ == NULL || !myXYZ->open())
    THROW(string("Can not open '") + (myXYZ ? myXYZ->getFilename() : "") +
          "' for " + getId() + ".");
}

void OutputXYZTrajectoryPos::doRun(int) {
  const Vector3DBlock *pos =
    (myMinimalImage ? app->outputCache.minimalPositions() : &app->positions);

  if (!myXYZ->write(*pos, app->topology->atoms, app->topology->atomTypes))
    THROW(string("Could not write ") + getId() + " '" + myXYZ->getFilename() +
          "'.");
}

void OutputXYZTrajectoryPos::doFinalize(int) {
  myXYZ->close();
}

Output *OutputXYZTrajectoryPos::doMake(const vector<Value> &values) const {
  return new OutputXYZTrajectoryPos(values[0], values[1], values[2]);
}

void OutputXYZTrajectoryPos::getParameters(vector<Parameter> &parameter)
const {
  parameter.push_back
    (Parameter(getId(), Value(myXYZ != NULL ? myXYZ->getFilename() : "",
                              ConstraintValueType::NotEmpty())));

  parameter.push_back
    (Parameter(keyword + "OutputFreq",
               Value(getOutputFreq(), ConstraintValueType::Positive())));

  parameter.push_back
    (Parameter(keyword + "MinimalImage", Value(myMinimalImage),
               Text("whether the coordinates should be transformed to minimal "
                    "image or not")));
}

bool OutputXYZTrajectoryPos::adjustWithDefaultParameters(
  vector<Value> &values, const Configuration *config) const {

  if (!checkParameterTypes(values)) return false;

  if (config->valid(InputOutputfreq::keyword) && !values[1].valid())
    values[1] = (*config)[InputOutputfreq::keyword];

  if (config->valid(InputMinimalImage::keyword) && !values[2].valid())
    values[2] = (*config)[InputMinimalImage::keyword];

  return checkParameters(values);
}
