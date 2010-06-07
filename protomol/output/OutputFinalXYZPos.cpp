#include <protomol/output/OutputFinalXYZPos.h>
#include <protomol/config/Configuration.h>
#include <protomol/output/OutputCache.h>
#include <protomol/module/MainModule.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/io/XYZWriter.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/base/Exception.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____ OutputFinalXYZPos
const string OutputFinalXYZPos::keyword("finXYZPosFile");

OutputFinalXYZPos::OutputFinalXYZPos() :
  Output(-1), filename(""), myMinimalImage(false) {}

OutputFinalXYZPos::OutputFinalXYZPos(const string &filename, bool minimal) :
  Output(-1), filename(filename), myMinimalImage(minimal) {}

void OutputFinalXYZPos::doFinalize(int step) {
  XYZWriter writer;
  if (!writer.open(filename))
    THROW(string("Can't open ") + getId() + " '" + filename + "'.");

  const Vector3DBlock *pos =
    (myMinimalImage ? app->outputCache.minimalPositions() : &app->positions);
  writer.setComment("Time : " + toString(app->outputCache.time()) +
                    ", step : " + toString(step) +
                    (myMinimalImage ? ", minimal Image" : "") + ".");

  if (!writer.write(*pos, app->topology->atoms, app->topology->atomTypes))
    THROW(string("Could not write ") + getId() + " '" + filename + "'.");
}

Output *OutputFinalXYZPos::doMake(const vector<Value> &values) const {
  return new OutputFinalXYZPos(values[0], values[1]);
}

void OutputFinalXYZPos::getParameters(vector<Parameter> &parameter) const {
  parameter.push_back
    (Parameter(getId(), Value(filename, ConstraintValueType::NotEmpty())));
  parameter.push_back
    (Parameter(keyword + "MinimalImage", Value(myMinimalImage),
               Text("whether the coordinates should be transformed to minimal"
                    " image or not")));
}

bool OutputFinalXYZPos::adjustWithDefaultParameters(vector<Value> &values,
                                                    const Configuration *config)
const {
  if (!checkParameterTypes(values)) return false;

  if (config->valid(InputMinimalImage::keyword) && !values[1].valid())
    values[1] = (*config)[InputMinimalImage::keyword];

  return checkParameters(values);
}
