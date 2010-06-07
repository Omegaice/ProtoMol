#ifndef NONBONDEDSIMPLEFULLFORCE_MODULE_H
#define NONBONDEDSIMPLEFULLFORCE_MODULE_H

#include <protomol/base/Module.h>

#include <string>

namespace ProtoMol {
  class ProtoMolApp;

  class NonbondedSimpleFullForceModule : public Module {
  public:
    const std::string getName() const {return "NonbondedSimpleFullForce";}
    void registerForces(ProtoMolApp *app);
  };
};

#endif // NONBONDEDSIMPLEFULLFORCE_MODULE_H
