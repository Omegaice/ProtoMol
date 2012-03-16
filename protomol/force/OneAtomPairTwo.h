/* -*- c++ -*- */
#ifndef ONEATOMPAIRTWO_H
#define ONEATOMPAIRTWO_H

#include <protomol/topology/Topology.h>
#include <protomol/config/Parameter.h>
#include <protomol/force/OneAtomPair.h>
#include <protomol/force/OneAtomContraints.h>

namespace ProtoMol {
  template<typename Boundary, typename SwitchA,
           typename ForceA, typename SwitchB,
           typename ForceB, typename Constraint = NoConstraint>
  class OneAtomPairTwo : public OneAtomPair<Boundary,SwitchA,ForceA,Constraint> {
    typedef OneAtomPair<Boundary,SwitchA,ForceA,Constraint> Base;
    
  public:
    OneAtomPairTwo() : Base() {
      
    }
    
    OneAtomPairTwo(ForceA f1, SwitchA sF1, ForceB f2, SwitchB sF2) 
      : Base( f1, sF1 ), SwitchFunctionB(sF2), ForceFunctionB(f2) {
      
        Base::mySquaredCutoff = std::max
                        (Cutoff<ForceA::CUTOFF>::cutoff(sF1, f1),
                         Cutoff<ForceB::CUTOFF>::cutoff(sF2, f2));
    }
    
    void doOneAtomPair(const int i, const int j) {
      if (Constraint::PRE_CHECK)
        if (!Constraint::check(Base::realTopo, i, j))
          return;

      // Get atom distance.
      Real distSquared;
      Vector3D diff(Base::realTopo->boundaryConditions.
                    minimalDifference((*Base::positions)[i], (*Base::positions)[j],
                                      distSquared));
      // Do switching function rough test, if necessary.
      if ((SwitchA::USE || SwitchB::USE ||
           ForceA::CUTOFF ||
           ForceB::CUTOFF) && distSquared > Base::mySquaredCutoff)
        return;

      // Check for an exclusion.
      int mi = Base::realTopo->atoms[i].molecule;
      int mj = Base::realTopo->atoms[j].molecule;
      bool same = (mi == mj);
      ExclusionClass excl =
        (same ? Base::realTopo->exclusions.check(i, j) : EXCLUSION_NONE);
      if (excl == EXCLUSION_FULL)
        return;

      // Calculate the force and energy.
      Real rDistSquared =
        ((ForceA::DIST_R2 ||
          ForceB::DIST_R2) ? 1.0 / distSquared : 1.0);
      Real energy1, force1, energy2 = 0, force2 = 0;
      Base::ForceFunction(energy1, force1, distSquared, rDistSquared,
                                  diff, Base::realTopo, i, j, excl);
      ForceFunctionB(energy2, force2, distSquared, rDistSquared,
                                   diff, Base::realTopo, i, j, excl);
      
      // Calculate the switched force and energy.
      if (SwitchA::MODIFY || SwitchB::MODIFY) {
        Real switchingValue, switchingDeriv;

        Base::SwitchFunction(switchingValue, switchingDeriv, distSquared);
        force1 = force1 * switchingValue - energy1 * switchingDeriv;
        energy1 = energy1 * switchingValue;

        SwitchFunctionB(switchingValue, switchingDeriv, distSquared);
        force2 = force2 * switchingValue - energy2 * switchingDeriv;
        energy2 = energy2 * switchingValue;
      }

      // Add this energy into the total system energy.
      Base::ForceFunction.accumulateEnergy(Base::energies, energy1);
      ForceFunctionB.accumulateEnergy(Base::energies, energy2);
      
      // Add this force into the atom forces.
      Vector3D fij(diff * (force1 + force2));
      (*Base::forces)[i] -= fij;
      (*Base::forces)[j] += fij;

      // compute the vector between molecular centers of mass
      if (!same && Base::energies->molecularVirial())
        // Add to the atomic and molecular virials
        Base::energies->
          addVirial(fij, diff, Base::realTopo->boundaryConditions.
                    minimalDifference(Base::realTopo->molecules[mi].position,
                                      Base::realTopo->molecules[mj].position));
      else if (Base::energies->virial())
        Base::energies->addVirial(fij, diff);
      
      // End of force computation.
      if (Constraint::POST_CHECK)
        Constraint::check(Base::realTopo, i, j, diff, energy1 + energy2, fij);
    }

    void getParameters(std::vector<Parameter> &parameters) const {
      Base::ForceFunction.getParameters(parameters);
      Base::SwitchFunction.getParameters(parameters);
      ForceFunctionB.getParameters(parameters);
      SwitchFunctionB.getParameters(parameters);
    }
    
    void postProcess(const GenericTopology *apptopo, ScalarStructure *appenergies){
		  Base::ForceFunction.postProcess(apptopo, appenergies);
		  ForceFunctionB.postProcess(apptopo, appenergies);
	  }

    void parallelPostProcess(const GenericTopology *apptopo, ScalarStructure *appenergies){
		  Base::ForceFunction.parallelPostProcess(apptopo, appenergies);
		  ForceFunctionB.parallelPostProcess(apptopo, appenergies);
	  }

    static unsigned int getParameterSize() {
      return
        ForceA::getParameterSize() + SwitchA::getParameterSize() +
        ForceB::getParameterSize() + SwitchB::getParameterSize();
    }

    bool doParallelPostProcess(){
      return Base::ForceFunction.doParallelPostProcess() || ForceFunctionB.doParallelPostProcess();
	  }
    
    static OneAtomPairTwo make(std::vector<Value> values) {
      unsigned int l1 = ForceA::getParameterSize();
      unsigned int l2 = SwitchA::getParameterSize() + l1;
      unsigned int l3 = ForceB::getParameterSize() + l2;

      std::vector<Value> F1(values.begin(), values.begin() + l1);
      std::vector<Value> S1(values.begin() + l1, values.begin() + l2);
      std::vector<Value> F2(values.begin() + l2, values.begin() + l3);
      std::vector<Value> S2(values.begin() + l3, values.end());

      return OneAtomPairTwo
        (ForceA::make(F1), SwitchA::make(S1), ForceB::make(F2), SwitchB::make(S2));
    }

    static std::string getId() {
      return
        Constraint::getPrefixId() + headString(ForceA::getId()) +
        Constraint::getPostfixId() + " " + 
        Constraint::getPrefixId() + headString(ForceB::getId()) +
        Constraint::getPostfixId() +
      
        (tailString(ForceA::getId()).empty() ? "" : " ") + tailString(ForceA::getId()) +
        (tailString(ForceB::getId()).empty() ? "" : " ") + tailString(ForceB::getId()) +
      
        std::string((!SwitchA::USE) ? std::string("") : std::string(" -switchingFunction " + SwitchA::getId())) +
        std::string((!SwitchB::USE) ? std::string("") : std::string(" -switchingFunction " + SwitchB::getId()));
    }
    
  protected:
    SwitchB SwitchFunctionB;
    ForceB ForceFunctionB;
  };
}

#endif /* ONEATOMPAIRTWO_H */
