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

//____ OutputDCDTrajectory
const string OutputDCDTrajectory::keyword("DCDFile");

OutputDCDTrajectory::OutputDCDTrajectory() :
  Output(), myDCD(NULL), myMinimalImage(false), myFrameOffset(0), myFilename("") {}

OutputDCDTrajectory::OutputDCDTrajectory(const string &filename, int freq,
                                         bool minimal, int frameoffs) :
  Output(freq), myDCD(NULL),//new DCDTrajectoryWriter(filename)),
  myMinimalImage(minimal), myFrameOffset(frameoffs), myFilename(filename) {

    
    //flag value
    report << plain << "DCD FrameOffset parameter set to " << myFrameOffset << "." << endr;
    
}

OutputDCDTrajectory::~OutputDCDTrajectory() {
  if (myDCD != NULL) delete myDCD;
}

void OutputDCDTrajectory::doInitialize() {
  
  //Get first frame (must exist or error)
  const int firstframe = toInt(app->config["firststep"]);
  
  report << debug(2) << "Firstframe " << firstframe << "." << endr;

  //open file
  //if myFrameOffset is zero default to overwrite data
  //note: now include "firstframe" data.
  if( myFrameOffset == 0){
    myDCD = new DCDTrajectoryWriter(myFilename, 1.0, firstframe);  //original call
   
  //else append to file
  }else{
   
    //file mode
    const std::ios::openmode mode = ios::binary | ios::ate | ios_base::in;

    //open DCD for writing with flags 'mode'
    myDCD = new DCDTrajectoryWriter(mode, myFrameOffset, myFilename);
   
  }
  
  if (myDCD == NULL || !myDCD->open())
    THROW(string("Can not open '") + (myDCD ? myDCD->getFilename() : "") +
          "' for " + getId() + ".");
}

void OutputDCDTrajectory::doRun(int) {
  const Vector3DBlock *pos =
    (myMinimalImage ? app->outputCache.minimalPositions() : &app->positions);

  if (!myDCD->write(*pos))
    THROW(string("Could not write ") + getId() + " '" +
          myDCD->getFilename() + "'.");
}

void OutputDCDTrajectory::doFinalize(int) {
  myDCD->close();
}

Output *OutputDCDTrajectory::doMake(const vector<Value> &values) const {
  return new OutputDCDTrajectory(values[0], values[1], values[2], values[3]);
}

void OutputDCDTrajectory::getParameters(vector<Parameter> &parameter) const {
  parameter.push_back
    (Parameter(getId(), Value(myDCD ? myDCD->getFilename() : "",
                              ConstraintValueType::NotEmpty())));
  parameter.push_back
    (Parameter(keyword + "OutputFreq",
               Value(getOutputFreq(), ConstraintValueType::Positive()) ));
  parameter.push_back
    (Parameter(keyword + "MinimalImage", Value(myMinimalImage),
               Text("whether the coordinates should be transformed to minimal "
                    "image or not")));
  parameter.push_back
    (Parameter(keyword + "FrameOffset", 
                Value(myFrameOffset, ConstraintValueType::NotNegative()), 0 ));
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

