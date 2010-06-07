#include <protomol/output/OutputFile.h>
#include <protomol/config/Configuration.h>
#include <protomol/module/MainModule.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____ Output
OutputFile::OutputFile() :
  Output(), myFilename(""), myCacheFreq(0), myCount(0), myCacheSize(0),
  myCloseTime(0.0) {}

OutputFile::OutputFile(const string &filename, int freq,
                       int cacheFreq, int cacheSize, Real closeTime) :
  Output(freq), myFilename(filename), myCacheFreq(cacheFreq), myCount(0),
  myCacheSize(cacheSize), myCloseTime(closeTime) {}

OutputFile::~OutputFile() {
  close();
}

bool OutputFile::open() {
  clearCache();
  myFile.open(myFilename.c_str(), ios::out | ios::trunc);
  if (!myFilename.empty())
    report << debug(5) << "Opening   file '" << myFilename << "'." << endr;
  myTimer.reset();
  return (bool)myFile;
}

bool OutputFile::close() {
  flushCache();
  if (myFile.is_open()) myFile.close();
  myTimer.reset();
  myTimer.start();
  if (!myFilename.empty())
    report << debug(5) << "Closing   file '" << myFilename << "'." << endr;
  return (bool)myFile;
}

bool OutputFile::reopen() {
  if (!myFile.is_open()) {
    if (!myFilename.empty())
      report << debug(5) << "Reopening file '" << myFilename << "'." << endr;
    myFile.open(myFilename.c_str(), ios::out | ios::app);
  }
  return (bool)myFile;
}

bool OutputFile::reclose() {
  Real t = myTimer.getTime().getRealTime();
  myTimer.reset();
  myTimer.start();
  if (t >= myCloseTime && myFile.is_open()) {
    if (!myFilename.empty())
      report << debug(5) << "Reclosing file '" << myFilename << "'." << endr;
    myFile.close();
  }
  return (bool)myFile;
}

void OutputFile::flushCache() {
  if (myBuffer.str().size() > 0) {
    if (!reopen())
      report << error << " Can not open \'" << myFilename << "\' for " <<
      this->getId() << "." << endr;
    myFile << myBuffer.str();
    reclose();
  }
  doFlushCache();
  clearCache();
}

void OutputFile::clearCache() {
  myCount = 0;
  myBuffer.str("");
}

void OutputFile::doRun(int step) {
  doRunCached(step);
  myCount = (myCount + 1) % myCacheFreq;
  if (myCount == 0 && static_cast<int>(myBuffer.str().size()) >= myCacheSize)
    flushCache();
}

void OutputFile::doFinalize(int) {
  close();
}

void OutputFile::getParameters(vector<Parameter> &parameter) const {
  parameter.push_back
    (Parameter(getId(), Value(myFilename, ConstraintValueType::NotEmpty())));
  parameter.push_back
    (Parameter(getId() + string("OutputFreq"),
               Value(getOutputFreq(), ConstraintValueType::Positive()),
               Text("output frequency")));
  parameter.push_back
    (Parameter(getId() + string("CacheFreq"),
               Value(myCacheFreq, ConstraintValueType::Positive()), 1,
               Text("frequency to flush the cache to file")));
  parameter.push_back
    (Parameter(getId() + string("CacheSize"),
               Value(myCacheSize, ConstraintValueType::NotNegative()), 0,
               Text("minimal size of cached data to flush to file")));
  parameter.push_back
    (Parameter(getId() + string("CloseTime"),
               Value(myCloseTime, ConstraintValueType::NotNegative()), 1.0,
               Text("minimal time interval between two writes to close the "
                    "file temporarily")));
}

bool OutputFile::adjustWithDefaultParameters(vector<Value> &values,
                                             const Configuration *config)
const {
  if (!checkParameterTypes(values))
    return false;
  if (config->valid(InputOutputfreq::keyword) && !values[1].valid())
    values[1] = (*config)[InputOutputfreq::keyword];
  return checkParameters(values);
}
