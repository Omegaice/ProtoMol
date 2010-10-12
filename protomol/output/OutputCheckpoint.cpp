#include <protomol/output/OutputCheckpoint.h>

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

#ifdef HAVE_LIBFAH
#include <fah/io/File.h>
typedef FAH::File fileStream;
#else
#include <fstream>
typedef std::fstream fileStream;
#endif

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

const string OutputCheckpoint::keyword("Checkpoint");

OutputCheckpoint::OutputCheckpoint() : mCurrent(0) {}

OutputCheckpoint::OutputCheckpoint(const std::string& name, int freq,
                                   int start, const std::string& posbase,
                                   const std::string& velbase) :
  Output(freq), mCurrent(start), mName(name),
  mPosBase(posbase), mVelBase(velbase) {
}

void OutputCheckpoint::doInitialize() {
  if (mPosBase == "") {
    const std::string temp = app->config["posfile"];
    
    mPosBase = temp.substr(0, temp.rfind('.') + 1);
  }
  
  if (mVelBase == "") {
    if (!app->config.valid("velfile")) {
      mVelBase = mPosBase;
    } else {
      const std::string temp = app->config["velfile"];
      
      mVelBase = temp.substr(0, temp.rfind('.') + 1);
    }
  }
}

void OutputCheckpoint::doIt(int step) {
  std::cout << "Checkpointing: Step " << step << ". . ." << std::flush;

  WritePositions(step);
  WriteVelocities(step);
  WriteConfig(step);

  // Remove old checkpoint files
  SystemUtilities::unlink(Append(Append(mPosBase, mCurrent - 1), ".pos"));
  SystemUtilities::unlink(Append(Append(mVelBase, mCurrent - 1), ".vel"));

  mCurrent += 1;
  
  std::cout << "done" << std::endl;
}

void OutputCheckpoint::doRun(int step) {
  const int firstStep = toInt(app->config["firststep"]);
  const int finalStep = firstStep + toInt(app->config["numsteps"]);
  
  if (step != firstStep && step != finalStep) {
    if (getOutputFreq() > 0 && (step % getOutputFreq()) == 0)
      doIt(step);
  }
}

void OutputCheckpoint::doFinalize(int step) {}

Output *OutputCheckpoint::doMake(const vector<Value> &values) const {
  return new OutputCheckpoint(values[0], values[1], values[2], values[3],
                              values[4]);
}

bool OutputCheckpoint::isIdDefined(const Configuration *config) const {
  return config->valid(getId());
}

void OutputCheckpoint::getParameters(vector<Parameter> &parameter) const {
  parameter.push_back
    (Parameter(getId(), Value(mName, ConstraintValueType::NotEmpty())));

  parameter.push_back
    (Parameter(keyword + "Freq",
               Value(getOutputFreq(), ConstraintValueType::Positive())));

  parameter.push_back
    (Parameter(keyword + "Start",
               Value(mCurrent, ConstraintValueType::NotNegative())));

  parameter.push_back
    (Parameter(keyword + "PosBase",
               Value(mPosBase, ConstraintValueType::NoConstraints())));

  parameter.push_back
    (Parameter(keyword + "VelBase",
               Value(mVelBase, ConstraintValueType::NoConstraints())));
}

bool OutputCheckpoint::
adjustWithDefaultParameters(vector<Value> &values,
                            const Configuration *config) const {
  if (!checkParameterTypes(values)) return false;

  if (!values[0].valid()) values[0] = mName;

  if (config->valid(InputOutputfreq::keyword) && !values[1].valid())
    values[1] = (*config)[InputOutputfreq::keyword];

  if (!values[2].valid()) values[2] = 0;
  if (!values[3].valid()) values[3] = "";
  if (!values[4].valid()) values[4] = "";

  return checkParameters(values);
}

void OutputCheckpoint::WritePositions(int step) {
  std::string posFile = Append(Append(mPosBase, mCurrent), ".pos");

  XYZWriter posWriter;
  if (!posWriter.open(posFile))
    THROWS("Can't open " << getId() << " '" << posFile + "'.");


  const Vector3DBlock *pos = &app->positions;
  posWriter.setComment("Time : " + toString(app->outputCache.time()) +
                       ", step : " + toString(step) +  ".");

  if (!posWriter.write(*pos, app->topology->atoms, app->topology->atomTypes))
    THROWS("Could not write " << getId() << " '" << posFile << "'.");
}

void OutputCheckpoint::WriteVelocities(int step) {
  std::string velFile = Append(Append(mVelBase, mCurrent), ".vel");

  XYZWriter velWriter;
  if (!velWriter.open(velFile))
    THROWS("Can't open " << getId() << " '" << velFile << "'.");

  velWriter.setComment("Time : " + toString(app->outputCache.time()) +
                       ", step : " + toString(step) + ".");

  if (!velWriter.write(*&app->velocities, app->topology->atoms,
                       app->topology->atomTypes))
    THROWS("Could not write " << getId() << " '" << velFile << "'.");
}

void OutputCheckpoint::WriteConfig(int step) {
  string confFile = mName + ".tmp";

  {
    CheckpointConfigWriter confWriter;
    if (!confWriter.open(confFile))
      THROWS("Can't open " << getId() << " '" << confFile << "'.");
    
    if (!confWriter.write(mCurrent, step, Random::Instance(), app->integrator))
      THROWS("Could not write " << getId() << " '" << confFile << "'.");

    confWriter.close();
  }

  SystemUtilities::rename(confFile, mName);
}
