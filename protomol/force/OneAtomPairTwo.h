/* -*- c++ -*- */
#ifndef ONEATOMPAIRTWO_H
#define ONEATOMPAIRTWO_H

#include <protomol/topology/Topology.h>
#include <protomol/config/Parameter.h>
#include <protomol/force/OneAtomContraints.h>

namespace ProtoMol {
  //____ OneAtomPairTwo

  /**
   * Computes the interaction for a given force between two atoms with the
   * template arguments defining the boundary conditions, two switching 
   * functions, two potentials and optional constraint.
   */
  template<typename Boundary, typename SwitchA,
           typename ForceA, typename SwitchB,
           typename ForceB, typename Constraint = NoConstraint>
  class OneAtomPairTwo {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Typedef & sub classes
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    typedef Boundary BoundaryConditions;
    // Make the boundary conditions visible

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    OneAtomPairTwo() {}

    OneAtomPairTwo(ForceA f1, SwitchA sF1,
                   ForceB f2, SwitchB sF2) :
      SwitchFunctionA(sF1), ForceFunctionA(f1),
      SwitchFunctionB(sF2), ForceFunctionB(f2),
      mySquaredCutoff(std::max
                      (Cutoff<ForceA::CUTOFF>::cutoff(sF1, f1),
                       Cutoff<ForceB::CUTOFF>::cutoff(sF2, f2)))
    {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class OneAtomPairTwo
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void initialize(const SemiGenericTopology<Boundary> *topo,
                    const Vector3DBlock *pos, Vector3DBlock *f,
                    ScalarStructure *e) {
      realTopo = topo;
      positions = pos;
      forces = f;
      energies = e;
    }

    // Computes the force and energy for atom i and j.
    void doOneAtomPair(const int i, const int j) {
      if (Constraint::PRE_CHECK)
        if (!Constraint::check(realTopo, i, j))
          return;

      // Get atom distance.
      Real distSquared;
      Vector3D diff(realTopo->boundaryConditions.
                    minimalDifference((*positions)[i], (*positions)[j],
                                      distSquared));
      // Do switching function rough test, if necessary.
      if ((SwitchA::USE || SwitchB::USE ||
           ForceA::CUTOFF ||
           ForceB::CUTOFF) && distSquared > mySquaredCutoff)
        return;

      // Check for an exclusion.
      int mi = realTopo->atoms[i].molecule;
      int mj = realTopo->atoms[j].molecule;
      bool same = (mi == mj);
      ExclusionClass excl =
        (same ? realTopo->exclusions.check(i, j) : EXCLUSION_NONE);
      if (excl == EXCLUSION_FULL)
        return;

      // Calculate the force and energy.
      Real rDistSquared =
        ((ForceA::DIST_R2 ||
          ForceB::DIST_R2) ? 1.0 / distSquared : 1.0);
      Real energy1, force1, energy2 = 0, force2 = 0;
      ForceFunctionA(energy1, force1, distSquared, rDistSquared,
                                  diff, realTopo, i, j, excl);
      ForceFunctionB(energy2, force2, distSquared, rDistSquared,
                                   diff, realTopo, i, j, excl);
      //Report::report << "\t"<<i << "\t"<<j<<Report::endr;
      // Calculate the switched force and energy.
      if (SwitchA::MODIFY || SwitchB::MODIFY) {
        Real switchingValue, switchingDeriv;

        SwitchFunctionA(switchingValue, switchingDeriv, distSquared);
        force1 = force1 * switchingValue - energy1 * switchingDeriv;
        energy1 = energy1 * switchingValue;

        SwitchFunctionB(switchingValue, switchingDeriv, distSquared);
        force2 = force2 * switchingValue - energy2 * switchingDeriv;
        energy2 = energy2 * switchingValue;
      }

      // Add this energy into the total system energy.
      ForceFunctionA.accumulateEnergy(energies, energy1);
      ForceFunctionB.accumulateEnergy(energies, energy2);
      // Add this force into the atom forces.
      Vector3D fij(diff * (force1 + force2));
      (*forces)[i] -= fij;
      (*forces)[j] += fij;

      // compute the vector between molecular centers of mass
      if (!same && energies->molecularVirial())
        // Add to the atomic and molecular virials
        energies->
          addVirial(fij, diff, realTopo->boundaryConditions.
                    minimalDifference(realTopo->molecules[mi].position,
                                      realTopo->molecules[mj].position));
      else if (energies->virial())
        energies->addVirial(fij, diff);
      // End of force computation.
      if (Constraint::POST_CHECK)
        Constraint::check(realTopo, i, j, diff, energy1 + energy2, fij);
    }

    void getParameters(std::vector<Parameter> &parameters) const {
      ForceFunctionA.getParameters(parameters);
      SwitchFunctionA.getParameters(parameters);
      ForceFunctionB.getParameters(parameters);
      SwitchFunctionB.getParameters(parameters);
    }
    
    void postProcess(const GenericTopology *apptopo, ScalarStructure *appenergies){
		  ForceFunctionA.postProcess(apptopo, appenergies);
		  ForceFunctionB.postProcess(apptopo, appenergies);
	  }

    void parallelPostProcess(const GenericTopology *apptopo, ScalarStructure *appenergies){
		  ForceFunctionA.parallelPostProcess(apptopo, appenergies);
		  ForceFunctionB.parallelPostProcess(apptopo, appenergies);
	  }

    static unsigned int getParameterSize() {
      return
        ForceA::getParameterSize() +
        SwitchA::getParameterSize() +
        ForceB::getParameterSize() +
        SwitchB::getParameterSize();
    }

    bool doParallelPostProcess(){
      return
        ForceFunctionA.doParallelPostProcess()
          || ForceFunctionB.doParallelPostProcess();
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
        (ForceA::make(F1), SwitchA::make(S1),
         ForceB::make(F2), SwitchB::make(S2));
    }

    static std::string getId() {
      return
        Constraint::getPrefixId() +
        headString(ForceA::getId()) +
        Constraint::getPostfixId() + " " + Constraint::getPrefixId() +
        headString(ForceB::getId()) +
        Constraint::getPostfixId() +
        (tailString(ForceA::getId()).empty() ? "" : " ") +
        tailString(ForceA::getId()) +
        (tailString(ForceB::getId()).empty() ? "" : " ") +
        tailString(ForceB::getId()) +
        std::string((!SwitchA::USE) ?
                    std::string("") :
                    std::string(" -switchingFunction " +
                                SwitchA::getId())) +
        std::string((!SwitchB::USE) ?
                    std::string("") :
                    std::string(" -switchingFunction " +
                                SwitchB::getId()));
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    const SemiGenericTopology<Boundary> *realTopo;
    const Vector3DBlock *positions;
    Vector3DBlock *forces;
    ScalarStructure *energies;
    SwitchA SwitchFunctionA;
    ForceA ForceFunctionA;
    SwitchB SwitchFunctionB;
    ForceB ForceFunctionB;
    Real mySquaredCutoff;
  };
}

#endif /* ONEATOMPAIRTWO_H */
