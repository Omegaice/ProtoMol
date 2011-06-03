/*  -*- c++ -*-  */
#ifndef NORMALMODEOPENMM_H
#define NORMALMODEOPENMM_H

#include <protomol/integrator/openMM/OpenMMIntegrator.h>
#include <protomol/integrator/normal/NormalModeUtilities.h>

namespace ProtoMol {

  class ScalarStructure;
  class ForceGroup;

  //__________________________________________________ NormalModeOpenMM
  class NormalModeOpenMM : public OpenMMIntegrator, public NormalModeUtilities {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    NormalModeOpenMM();
    NormalModeOpenMM(Real timestep, int firstmode, int nummode, Real gamma, Real temperature, 
                            Real minimlim, ForceGroup *overloadedForces);
    ~NormalModeOpenMM(); 

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class NormalModeOpenMM
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const{return keyword;}
    virtual void getParameters(std::vector<Parameter>& parameters) const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class OpenMMIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    virtual unsigned int getParameterSize() const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Integrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void initialize(ProtoMolApp* app);
    virtual void run(int numTimesteps);
  protected:
    //virtual void addModifierAfterInitialize();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class STSIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    virtual STSIntegrator* doMake(const std::vector<Value>& values, ForceGroup* fg)const;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // My data members
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;

  private:
    Real minLim;
    NormalModeUtilities *myPreviousNormalMode;

  };

}

#endif


