#include <protomol/force/GB/GBPartialSum.h>

#include <protomol/parallel/Parallel.h>

using namespace ProtoMol;

const string GBPartialSum::keyword("GBPartialSum");

GBPartialSum::GBPartialSum(){}

void GBPartialSum::operator()(Real &energy, Real &force,
                Real distSquared, Real rDistSquared, const Vector3D &diff,
                const GenericTopology *topo,
                int atom1, int atom2, ExclusionClass excl) const {
  
  Real bornRad_i = topo->atoms[atom1].myGBSA_T->bornRad;
  Real bornRad_l = topo->atoms[atom2].myGBSA_T->bornRad;
  
  Real scaledCharge_i = topo->atoms[atom1].scaledCharge;
  Real scaledCharge_l = topo->atoms[atom2].scaledCharge;
  
  Real ril = sqrt(distSquared);
  
  Real expterm = (ril*ril)/(4.0*bornRad_i*bornRad_l);
  Real filGB = sqrt(ril*ril + bornRad_i*bornRad_l*exp(-expterm));
  
  topo->atoms[atom1].myGBSA_T->partialGBForceTerms += scaledCharge_i*scaledCharge_l*(1/(filGB*filGB))*0.5*(1/filGB)*exp(-expterm)*(bornRad_l + (ril*ril)/(4.0*bornRad_i));
  topo->atoms[atom2].myGBSA_T->partialGBForceTerms += scaledCharge_i*scaledCharge_l*(1/(filGB*filGB))*0.5*(1/filGB)*exp(-expterm)*(bornRad_i + (ril*ril)/(4.0*bornRad_l));
}

void GBPartialSum::accumulateEnergy(ScalarStructure *energies, Real energy) {
  
}

Real GBPartialSum::getEnergy(const ScalarStructure *energies) {
  return 0;
}

void GBPartialSum::postProcess(const GenericTopology *topo, ScalarStructure *energies) {
  const unsigned int atomnumber = topo->atoms.size();
  
  for( unsigned int i=0; i<atomnumber; i++)
    topo->atoms[i].myGBSA_T->havePartialGBForceTerms = true;
}

void GBPartialSum::parallelPostProcess(const GenericTopology *topo, ScalarStructure *energies) {
  const unsigned int atomnumber = topo->atoms.size();
  
  // Copy Radii
  Real *partial = new Real[ atomnumber ];
  
  //put radii into array
  for( unsigned int i = 0; i < atomnumber; i++ ){
    partial[i] = topo->atoms[i].myGBSA_T->partialGBForceTerms;
  }
  
  //sum accross nodes
  Parallel::reduce(partial, partial + atomnumber); //reduceSlaves only?
  
  //put radii back
  for( unsigned int i = 0; i < atomnumber; i++ ){
    topo->atoms[i].myGBSA_T->partialGBForceTerms = partial[i];
  }
  
  delete [] partial;
}

bool GBPartialSum::doParallelPostProcess() {
  return true;
}

// Parsing
std::string GBPartialSum::getId() {
  return keyword;
}

unsigned int GBPartialSum::getParameterSize() {
  return 0;
}

void GBPartialSum::getParameters(std::vector<Parameter> &) const {

}

GBPartialSum GBPartialSum::make(const std::vector<Value> &) {
  return GBPartialSum();
}