/* -*- c++ -*- */
#ifndef ONEATOMPAIRTHREE_H
#define ONEATOMPAIRTHREE_H

#include <protomol/topology/Topology.h>
#include <protomol/config/Parameter.h>
#include <protomol/force/OneAtomPair.h>
#include <protomol/force/OneAtomContraints.h>

namespace ProtoMol {
  template<typename Boundary, typename SwitchA,
           typename ForceA, typename SwitchB,
           typename ForceB, typename SwitchC,
           typename ForceC, typename Constraint = NoConstraint>
  class OneAtomPairThree : public OneAtomPair<Boundary,SwitchA,ForceA,Constraint> {
    typedef OneAtomPair<Boundary,SwitchA,ForceA,Constraint> Base;
    
  public:
    OneAtomPairThree() : Base() {
      
    }
    
    OneAtomPairThree(ForceA f1, SwitchA sF1,
                     ForceB f2, SwitchB sF2,
                     ForceC f3, SwitchC sF3 ) 
      : Base( f1, sF1 ), SwitchFunctionB(sF2), ForceFunctionB(f2),
        SwitchFunctionC(sF3), ForceFunctionC(f3) {
          
        Base::mySquaredCutoff = std::max(
                             std::max
                             (Cutoff<ForceA::CUTOFF>::cutoff(sF1, f1),
                                Cutoff<ForceB::CUTOFF>::cutoff(sF2, f2)),
                                         Cutoff<ForceC::CUTOFF>::cutoff(sF3, f3));
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
      if ((SwitchA::USE || SwitchB::USE || SwitchC::USE ||
           ForceA::CUTOFF || ForceB::CUTOFF || ForceC::CUTOFF ) 
           && distSquared > Base::mySquaredCutoff)
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
            ForceB::DIST_R2 ||
              ForceC::DIST_R2 ) ? 1.0 / distSquared : 1.0);
      Real energy1, force1, energy2 = 0, force2 = 0, energy3 = 0, force3 = 0;
      Base::ForceFunction(energy1, force1, distSquared, rDistSquared,
                                  diff, Base::realTopo, i, j, excl);
      ForceFunctionB(energy2, force2, distSquared, rDistSquared,
                                   diff, Base::realTopo, i, j, excl);
      ForceFunctionC(energy3, force3, distSquared, rDistSquared,
                                   diff, Base::realTopo, i, j, excl);
      
      // Calculate the switched force and energy.
      if (SwitchA::MODIFY || 
            SwitchB::MODIFY || 
              SwitchC::MODIFY) {
        Real switchingValue, switchingDeriv;

        Base::SwitchFunction(switchingValue, switchingDeriv, distSquared);
        force1 = force1 * switchingValue - energy1 * switchingDeriv;
        energy1 = energy1 * switchingValue;

        SwitchFunctionB(switchingValue, switchingDeriv, distSquared);
        force2 = force2 * switchingValue - energy2 * switchingDeriv;
        energy2 = energy2 * switchingValue;

        SwitchFunctionC(switchingValue, switchingDeriv, distSquared);
        force3 = force3 * switchingValue - energy3 * switchingDeriv;
        energy3 = energy3 * switchingValue;
      }

      // Add this energy into the total system energy.
      Base::ForceFunction.accumulateEnergy(Base::energies, energy1);
      ForceFunctionB.accumulateEnergy(Base::energies, energy2);
      ForceFunctionC.accumulateEnergy(Base::energies, energy3);
      
      // Add this force into the atom forces.
      Vector3D fij(diff * (force1 + force2 + force3));
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
        Constraint::check(Base::realTopo, i, j, diff, energy1 + energy2 + energy3, fij);
    }

    virtual void getParameters(std::vector<Parameter> &parameters) const {
      Base::ForceFunction.getParameters(parameters);
      Base::SwitchFunction.getParameters(parameters);
      ForceFunctionB.getParameters(parameters);
      SwitchFunctionB.getParameters(parameters);
      ForceFunctionC.getParameters(parameters);
      SwitchFunctionC.getParameters(parameters);
    }
    
    virtual void preProcess(const GenericTopology *apptopo, const Vector3DBlock *positions){
      Base::ForceFunction.preProcess(apptopo, positions);
		  ForceFunctionB.preProcess(apptopo, positions);
		  ForceFunctionC.preProcess(apptopo, positions);
    }
    
    virtual void postProcess(const GenericTopology *apptopo, ScalarStructure *appenergies){
		  Base::ForceFunction.postProcess(apptopo, appenergies);
		  ForceFunctionB.postProcess(apptopo, appenergies);
		  ForceFunctionC.postProcess(apptopo, appenergies);
	  }
    
    virtual void parallelPostProcess(const GenericTopology *apptopo, ScalarStructure *appenergies){
		  Base::ForceFunction.parallelPostProcess(apptopo, appenergies);
		  ForceFunctionB.parallelPostProcess(apptopo, appenergies);
		  ForceFunctionC.parallelPostProcess(apptopo, appenergies);
	  }

    virtual bool doParallelPostProcess(){
      return Base::ForceFunction.doParallelPostProcess() || ForceFunctionB.doParallelPostProcess()
          || ForceFunctionC.doParallelPostProcess();
	  }

    static unsigned int getParameterSize() {
      return
        ForceA::getParameterSize() +
        SwitchA::getParameterSize() +
        ForceB::getParameterSize() +
        SwitchB::getParameterSize() +
        ForceC::getParameterSize() +
        SwitchC::getParameterSize();
    }

    static OneAtomPairThree make(std::vector<Value> values) {
      unsigned int l1a = ForceA::getParameterSize();
      unsigned int l1b = SwitchA::getParameterSize() + l1a;
      unsigned int l2a = ForceB::getParameterSize() + l1b;
      unsigned int l2b = SwitchB::getParameterSize() + l2a;
      unsigned int l3a = ForceC::getParameterSize() + l2b;

      std::vector<Value> F1(values.begin(), values.begin() + l1a);
      std::vector<Value> S1(values.begin() + l1a, values.begin() + l1b);
      std::vector<Value> F2(values.begin() + l1b, values.begin() + l2a);
      std::vector<Value> S2(values.begin() + l2a, values.begin() + l2b);
      std::vector<Value> F3(values.begin() + l2b, values.begin() + l3a);
      std::vector<Value> S3(values.begin() + l3a, values.end());

      return OneAtomPairThree (ForceA::make(F1), SwitchA::make(S1),
          ForceB::make(F2), SwitchB::make(S2),
          ForceC::make(F3), SwitchC::make(S3)
      );
    }

    static std::string getId() {
      return
        Constraint::getPrefixId() + headString(ForceA::getId()) +
        Constraint::getPostfixId() + " " + 
        Constraint::getPrefixId() + headString(ForceB::getId()) +
        Constraint::getPostfixId() + " " +  
        Constraint::getPrefixId() + headString(ForceC::getId()) + 
        Constraint::getPostfixId() +

        (tailString(ForceA::getId()).empty() ? "" : " ") + tailString(ForceA::getId()) +
        (tailString(ForceB::getId()).empty() ? "" : " ") + tailString(ForceB::getId()) +
        (tailString(ForceC::getId()).empty() ? "" : " ") + tailString(ForceC::getId()) +

        std::string((!SwitchA::USE) ? std::string("") : std::string(" -switchingFunction " + SwitchA::getId())) +
        std::string((!SwitchB::USE) ? std::string("") : std::string(" -switchingFunction " + SwitchB::getId())) +
        std::string((!SwitchC::USE) ? std::string("") : std::string(" -switchingFunction " + SwitchC::getId()));

    }
  protected:
    SwitchB SwitchFunctionB;
    ForceB ForceFunctionB;
    SwitchC SwitchFunctionC;
    ForceC ForceFunctionC;
  };
}

#endif /* ONEATOMPAIRTHREE_H */
