#include <protomol/base/Module.h>

namespace ProtoMol {
  int Module::getPriority() const {
    return 0;
  }
  
  const std::string Module::getHelp() const {
    return "";
  }
  
  void Module::getDependencies(module_deps_t &/*deps*/) const {
  
  }
    
  void Module::configure(ProtoMolApp */*app*/) {
  
  }
  
  void Module::read(ProtoMolApp */*app*/) {
  
  }
  
  void Module::registerForces(ProtoMolApp */*app*/) {
  
  }
  
  void Module::postBuild(ProtoMolApp */*app*/) {
  
  }
  
  void Module::addModifiers(ProtoMolApp */*app*/) {
  
  }
    
  void Module::init(ProtoMolApp */*app*/) {
  
  }
  
  bool moduleLess::operator()(const Module *m1, const Module *m2) const {
    if (!m1 || !m2) THROW("null pointer");
    
    if (m1->getPriority() < m2->getPriority()) return true;
    if (m1->getPriority() > m2->getPriority()) return false;
    return m1->getName() < m2->getName();
  }
}