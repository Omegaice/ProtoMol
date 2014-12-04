#include <protomol/output/OutputDCDTrajectory.h>
#include <protomol/config/Configuration.h>
#include <protomol/output/OutputCache.h>
#include <protomol/module/MainModule.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/io/DCDTrajectoryWriter.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/base/Exception.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

const string OutputDCDTrajectory::keyword("DCDFile");


OutputDCDTrajectory::OutputDCDTrajectory() :
  dCD(0), minimalImage(false), frameOffset(0), cachesize(1), firstWrite(true), cacheoffset(0) {}


OutputDCDTrajectory::OutputDCDTrajectory(const string &filename, int freq,
                                         bool minimal, int frameoffs) :
  Output(freq), dCD(0), minimalImage(minimal), frameOffset(frameoffs),
  filename(filename), cachesize(10), firstWrite(true), cacheoffset(0) {

  report << plain << "DCD FrameOffset parameter set to "
         << frameOffset << "." << endr;
}


OutputDCDTrajectory::~OutputDCDTrajectory() {
  if (dCD){
    if(cachedCoords.size() > 0){
      if (!dCD->write(cachedCoords))
        THROWS("Could not write " << getId() << " '" << dCD->getFilename() << "'.");
      cachedCoords.clear();
    }
    delete dCD;
  }
}


void OutputDCDTrajectory::doInitialize() {
  // Get first frame (must exist or error)
  const int firstframe = toInt(app->config["firststep"]);

  report << debug(2) << "Firstframe " << firstframe << "." << endr;

  // open fil
  // if frameOffset is zero default to overwrite dat
  // note: now include "firstframe" data.
  if (!frameOffset)
    dCD = new DCDTrajectoryWriter(filename, 1.0, firstframe); // original call

  else { // else append to fil

    // file mod
    const std::ios::openmode mode = ios::binary | ios::ate | ios_base::in;

    // open DCD for writing with flags 'mode'
    dCD = new DCDTrajectoryWriter(mode, frameOffset, filename);

  }

  if (!dCD || !dCD->open())
    THROWS("Can not open '" << (dCD ? dCD->getFilename() : "")
           << "' for " << getId() << ".");
}


void OutputDCDTrajectory::doRun(long) {
  
  //don't write first frame if checkpoint re-start
  if(firstWrite && frameOffset != 0 ){
    //one shot
    firstWrite = false;
    return;
  }
  
  if(firstWrite){
    cacheoffset = 1;
    firstWrite = false;
  }

  const Vector3DBlock *pos =
    (minimalImage ? app->outputCache.getMinimalPositions() : &app->positions);
  
  //cache data
  cachedCoords.push_back(*pos);
  
  if (cachedCoords.size() >= cachesize + cacheoffset){
    cacheoffset = 0;
    if (!dCD->write(cachedCoords))
      THROWS("Could not write " << getId() << " '" << dCD->getFilename() << "'.");
    cachedCoords.clear();
  }
}


void OutputDCDTrajectory::doFinalize(long) {
  if(cachedCoords.size() > 0){
    if (!dCD->write(cachedCoords))
      THROWS("Could not write " << getId() << " '" << dCD->getFilename() << "'.");
    cachedCoords.clear();
  }
  dCD->close();
}


Output *OutputDCDTrajectory::doMake(const vector<Value> &values) const {
  return new OutputDCDTrajectory(values[0], values[1], values[2], values[3]);//, values[4]);
}


void OutputDCDTrajectory::getParameters(vector<Parameter> &parameter) const {
  parameter.push_back
    (Parameter(getId(), Value(dCD ? dCD->getFilename() : "",
                              ConstraintValueType::NotEmpty())));
  Output::getParameters(parameter);
  parameter.push_back
    (Parameter(keyword + "MinimalImage", Value(minimalImage),
               Text("whether the coordinates should be transformed to minimal "
                    "image or not")));
  parameter.push_back
    (Parameter(keyword + "FrameOffset",
                Value(frameOffset, ConstraintValueType::NotNegative()), 0 ));
  /*parameter.push_back
    (Parameter(keyword + "CacheSize",
             Value(cachesize, ConstraintValueType::NotNegative()), 1 ));*/
}


bool OutputDCDTrajectory::adjustWithDefaultParameters(
  vector<Value> &values, const Configuration *config) const {

  if (!checkParameterTypes(values)) return false;

  if (config->valid(InputOutputfreq::keyword) && !values[1].valid())
    values[1] = (*config)[InputOutputfreq::keyword];

  if (config->valid(InputMinimalImage::keyword) && !values[2].valid())
    values[2] = (*config)[InputMinimalImage::keyword];

  return checkParameters(values);
}

