#include "CheckpointModule.h"

#include "IOModule.h"

#include <protomol/ProtoMolApp.h>
#include <protomol/io/CheckpointConfigReader.h>
#include <protomol/output/OutputCheckpoint.h>
#include <protomol/factory/OutputFactory.h>

using namespace std;
using namespace ProtoMol;

void CheckpointModule::init(ProtoMolApp *app) {
  OutputFactory &f = app->outputFactory;
  f.registerExemplar(new OutputCheckpoint());
}

void CheckpointModule::configure(ProtoMolApp *app) {
  Configuration &config = app->config;

  if (!(enabled = config.valid("Checkpoint"))) return;

  config["CheckpointPosBase"] = "";
  config["CheckpointVelBase"] = "";

  if( config.valid("gromacstprfile") ) {
    const string temp = config["gromacstprfile"];
    config["CheckpointPosBase"] = temp.substr(0, temp.rfind('.')+1);
    config["CheckpointVelBase"] = temp.substr(0, temp.rfind('.')+1);
  }else{
    const string temp = config["posfile"];
    config["CheckpointPosBase"] = temp.substr(0, temp.rfind('.') + 1);

    if (!config.valid("velfile")) {
        config["CheckpointVelBase"] = config["CheckpointPosBase"];
    } else {
        const string tempv = config["velfile"];
        config["CheckpointVelBase"] = tempv.substr(0, tempv.rfind('.') + 1);
    }
  }
}

void CheckpointModule::read(ProtoMolApp *app) {
  if (!enabled) return;

  Configuration &config = app->config;

  std::cout << "Reading Checkpoint Base Data" << std::endl;

  CheckpointConfigReader confReader;
  if (confReader.open(config["Checkpoint"], ios::in)){
    confReader.readBase(config, Random::Instance());
  }
}

void CheckpointModule::postBuild(ProtoMolApp *app) {
  if (!enabled) return;

  Configuration &config = app->config;

  std::cout << "Reading Checkpoint Integrator Data" << std::endl;
  // Load integrator data
  CheckpointConfigReader confReader;
  if (confReader.open(config["Checkpoint"], ios::in)){
    confReader.readIntegrator(app->integrator);
  }
}

string CheckpointModule::WithoutExt(const string &path){
  return path.substr(0, path.rfind(".") + 1);
}
