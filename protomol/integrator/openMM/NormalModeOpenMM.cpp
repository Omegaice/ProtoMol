#include <protomol/integrator/openMM/NormalModeOpenMM.h>
#include <protomol/base/Report.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/ProtoMolApp.h>

//#include "ModifierForceProjection.h"

using namespace ProtoMol::Report;

namespace ProtoMol {
  //__________________________________________________ NormalModeOpenMM

  const string NormalModeOpenMM::keyword( "NormalModeOpenMM" );

  NormalModeOpenMM::NormalModeOpenMM() : OpenMMIntegrator(), NormalModeUtilities()
  {
  }

  NormalModeOpenMM::NormalModeOpenMM(Real timestep, int firstmode, int nummode, Real gamma, Real temperature,
      Real minimlim, ForceGroup *overloadedForces)
    : OpenMMIntegrator(timestep,overloadedForces), NormalModeUtilities( firstmode, nummode, gamma, 1234, temperature),
      minLim(minimlim)
  {
  }

  NormalModeOpenMM::~NormalModeOpenMM()
  {

  }

  void NormalModeOpenMM::initialize(ProtoMolApp* app){

    report << plain << "OpenMM NML Vector information: Vector number " << app->eigenInfo.myNumEigenvectors << ", length " << app->eigenInfo.myEigenvectorLength << "." << endr;
    //NM initialization
    NormalModeUtilities::initialize((int)app->positions.size(), app,
				    myForces, NO_NM_FLAGS);

    //Set up minimum limit
    app->eigenInfo.myMinimumLimit = minLim;
    
    //Set number of eigenvectors in use
    app->eigenInfo.myNumUsedEigenvectors = _rfM;

    //initialize base
    OpenMMIntegrator::initialize(app);
    
    initializeForces();


    
  }

  void NormalModeOpenMM::run(int numTimesteps) {

    if( numTimesteps < 1 ){
        return;
    }

    //check valid eigenvectors
    if(*Q == NULL){
        report << error << "No Eigenvectors for NormalMode integrator."<<endr;
    }
#if 0
      // Retrive data
      const OpenMM::State state1 = context->getState(OpenMM::State::Positions |
                                                    OpenMM::State::Velocities |
                                                    OpenMM::State::Forces |
                                                    OpenMM::State::Energy);
      openMMvelocities = state1.getVelocities();

      report.precision(5);
      report << plain <<  "Vstart ";
      for(int i=0; i<_N;i++){
        report << plain << openMMvelocities[i][0] << ", " << openMMvelocities[i][1] << ", " << openMMvelocities[i][2] << "." << endr;
      }
      
      tempV3DBlk.resize(_N);

      for (unsigned int i = 0; i < _N; ++i){
       for (int j = 0; j < 3; j++){
         tempV3DBlk[i].c[j] = openMMvelocities[i][j];
        }
      }

      //test, find modes
      std::cout << "Projv ";
      //subspaceForce(myForces, myForces);
      //subspaceVelocity(&tempV3DBlk,&tempV3DBlk);
      nonSubspacePosition(&tempV3DBlk,&tempV3DBlk);
      //subspaceForce(&tempV3DBlk,&tempV3DBlk);
      //nonSubspaceForce(&tempV3DBlk,&tempV3DBlk);
      for(int i=0; i<10; i++){
          std::cout << tmpC[i] << ", ";
      }
      std::cout << endl;

      report << plain <<  "Vproj  ";
      for(int i=0; i<_N;i++){
        report << plain << tempV3DBlk[i][0] << ", " << tempV3DBlk[i][1] << ", " << tempV3DBlk[i][2] << "." << endr;
      }
#endif
    OpenMMIntegrator::run(numTimesteps);
    
      //test, find modes
      //std::cout << "Proj ";// << (*Q)[0] << " " << OpenMMIntegrator::myForces->c[0] << ", " << myForces->c[0] <<": ";
      //subspaceForce(myForces, myForces);
     //subspaceForce(&app->velocities,&app->velocities);
      //for(int i=0; i<10; i++){
      //    std::cout << tmpC[i] << ", ";
      //}
      //std::cout << endl;
#if 0
      report << plain <<  "Vcuda  ";
      for(int i=0; i<_N;i++){
        report << plain << openMMvelocities[i][0] << ", " << openMMvelocities[i][1] << ", " << openMMvelocities[i][2] << "." << endr;
      }
#endif
  }

  void NormalModeOpenMM::getParameters(vector<Parameter>& parameters) const {
    OpenMMIntegrator::getParameters(parameters);

    parameters.push_back(Parameter("firstmode",Value(firstMode,ConstraintValueType::NoConstraints()),1,Text("First mode to use in set")));
    parameters.push_back(Parameter("numbermodes",Value(numMode,ConstraintValueType::NoConstraints()),1,Text("Number of modes propagated")));
    parameters.push_back(Parameter("gamma",Value(NormalModeUtilities::myGamma*(1000 * Constant::INV_TIMEFACTOR),ConstraintValueType::NotNegative()),80.0,Text("Langevin Gamma")));
    parameters.push_back(Parameter("temperature",Value(myTemp,ConstraintValueType::NotNegative()),300.0,Text("Langevin temperature")));
    parameters.push_back(Parameter("minimlim",Value(minLim,ConstraintValueType::NotNegative()),0.1,Text("Minimizer target PE difference kcal mole^{-1}")));
  }

  STSIntegrator* NormalModeOpenMM::doMake(const vector<Value>& values,ForceGroup* fg) const {
    const unsigned int numPar = OpenMMIntegrator::getParameterSize();

    std::vector<Value> ommValues(values.begin() + 1, values.begin() + numPar);

    NormalModeOpenMM* myIntegP = new NormalModeOpenMM(
        values[0], values[numPar+0], values[numPar+1], values[numPar+2], values[numPar+3], values[numPar+4], fg
    );

    myIntegP->OpenMMIntegrator::setupValues(ommValues);

    return (STSIntegrator*)myIntegP;
  }

  unsigned int NormalModeOpenMM::getParameterSize() const {
    return OpenMMIntegrator::getParameterSize() + 5;
  }
}

