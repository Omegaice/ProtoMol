/* -*- c++ -*- */
#ifndef ONEATOMPAIRTHREE_H
#define ONEATOMPAIRTHREE_H

#include <protomol/topology/Topology.h>
#include <protomol/config/Parameter.h>
#include <protomol/force/OneAtomContraints.h>

namespace ProtoMol {
  //____ OneAtomPairThree

  /**
   * Computes the interaction for a given force between two atoms with the
   * template arguments defining the boundary conditions, three switching 
   * functions, three potentials and optional constraint.
   */
  template<typename Boundary, typename SwitchA,
           typename ForceA, typename SwitchB,
           typename ForceB, typename SwitchC,
           typename ForceC, typename Constraint = NoConstraint>
  class OneAtomPairThree {
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
    OneAtomPairThree() {}

    OneAtomPairThree(ForceA f1, SwitchA sF1,
                   ForceB f2, SwitchB sF2,
                   ForceC f3, SwitchC sF3 ) :
      switchingFunctionFirst(sF1), nonbondedForceFunctionFirst(f1),
      switchingFunctionSecond(sF2), nonbondedForceFunctionSecond(f2),
      switchingFunctionThird(sF3), nonbondedForceFunctionThird(f3),
      mySquaredCutoff(std::max(
                      std::max
                      (Cutoff<ForceA::CUTOFF>::cutoff(sF1, f1),
                        Cutoff<ForceB::CUTOFF>::cutoff(sF2, f2)),
                          Cutoff<ForceC::CUTOFF>::cutoff(sF3, f3))
                       )
    {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class OneAtomPairThree
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
      if ((SwitchA::USE || SwitchB::USE || SwitchC::USE ||
           ForceA::CUTOFF || ForceB::CUTOFF || ForceC::CUTOFF ) 
           && distSquared > mySquaredCutoff)
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
            ForceB::DIST_R2 ||
              ForceC::DIST_R2 ) ? 1.0 / distSquared : 1.0);
      Real energy1, force1, energy2 = 0, force2 = 0, energy3 = 0, force3 = 0;
      nonbondedForceFunctionFirst(energy1, force1, distSquared, rDistSquared,
                                  diff, realTopo, i, j, excl);
      nonbondedForceFunctionSecond(energy2, force2, distSquared, rDistSquared,
                                   diff, realTopo, i, j, excl);
      nonbondedForceFunctionThird(energy3, force3, distSquared, rDistSquared,
                                   diff, realTopo, i, j, excl);
      //Report::report << "\t"<<i << "\t"<<j<<Report::endr;
      // Calculate the switched force and energy.
      if (SwitchA::MODIFY || 
            SwitchB::MODIFY || 
              SwitchC::MODIFY) {
        Real switchingValue, switchingDeriv;

        switchingFunctionFirst(switchingValue, switchingDeriv, distSquared);
        force1 = force1 * switchingValue - energy1 * switchingDeriv;
        energy1 = energy1 * switchingValue;

        switchingFunctionSecond(switchingValue, switchingDeriv, distSquared);
        force2 = force2 * switchingValue - energy2 * switchingDeriv;
        energy2 = energy2 * switchingValue;

        switchingFunctionThird(switchingValue, switchingDeriv, distSquared);
        force3 = force3 * switchingValue - energy3 * switchingDeriv;
        energy3 = energy3 * switchingValue;
      }

      // Add this energy into the total system energy.
      nonbondedForceFunctionFirst.accumulateEnergy(energies, energy1);
      nonbondedForceFunctionSecond.accumulateEnergy(energies, energy2);
      nonbondedForceFunctionThird.accumulateEnergy(energies, energy3);
      // Add this force into the atom forces.
      Vector3D fij(diff * (force1 + force2 + force3));
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
        Constraint::check(realTopo, i, j, diff, energy1 + energy2 + energy3, fij);
    }

    void getParameters(std::vector<Parameter> &parameters) const {
      nonbondedForceFunctionFirst.getParameters(parameters);
      switchingFunctionFirst.getParameters(parameters);
      nonbondedForceFunctionSecond.getParameters(parameters);
      switchingFunctionSecond.getParameters(parameters);
      nonbondedForceFunctionThird.getParameters(parameters);
      switchingFunctionThird.getParameters(parameters);
    }
    
    void postProcess(const GenericTopology *apptopo, ScalarStructure *appenergies){
		  nonbondedForceFunctionFirst.postProcess(apptopo, appenergies);
		  nonbondedForceFunctionSecond.postProcess(apptopo, appenergies);
		  nonbondedForceFunctionThird.postProcess(apptopo, appenergies);
	  }
    
    void parallelPostProcess(const GenericTopology *apptopo, ScalarStructure *appenergies){
		  nonbondedForceFunctionFirst.parallelPostProcess(apptopo, appenergies);
		  nonbondedForceFunctionSecond.parallelPostProcess(apptopo, appenergies);
		  nonbondedForceFunctionThird.parallelPostProcess(apptopo, appenergies);
	  }

    bool doParallelPostProcess(){
      return
        nonbondedForceFunctionFirst.doParallelPostProcess()
        || nonbondedForceFunctionSecond.doParallelPostProcess()
          || nonbondedForceFunctionThird.doParallelPostProcess();
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

      return OneAtomPairThree
        (ForceA::make(F1), SwitchA::make(S1),
          ForceB::make(F2), SwitchB::make(S2),
            ForceC::make(F3), SwitchC::make(S3));
    }

    static std::string getId() {
      return
        Constraint::getPrefixId() +
        headString(ForceA::getId()) +
        Constraint::getPostfixId() + " " + 
        Constraint::getPrefixId() +
        headString(ForceB::getId()) +
        Constraint::getPostfixId() + " " + 
        Constraint::getPrefixId() +
        headString(ForceC::getId()) +
        Constraint::getPostfixId() +

        (tailString(ForceA::getId()).empty() ? "" : " ") +
        tailString(ForceA::getId()) +
        (tailString(ForceB::getId()).empty() ? "" : " ") +
        tailString(ForceB::getId()) +
        (tailString(ForceC::getId()).empty() ? "" : " ") +
        tailString(ForceC::getId()) +

        std::string((!SwitchA::USE) ?
                    std::string("") :
                    std::string(" -switchingFunction " +
                                SwitchA::getId())) +
        std::string((!SwitchB::USE) ?
                    std::string("") :
                    std::string(" -switchingFunction " +
                                SwitchB::getId())) +
        std::string((!SwitchC::USE) ?
                    std::string("") :
                    std::string(" -switchingFunction " +
                                SwitchC::getId()));

    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    const SemiGenericTopology<Boundary> *realTopo;
    const Vector3DBlock *positions;
    Vector3DBlock *forces;
    ScalarStructure *energies;
    SwitchA switchingFunctionFirst;
    ForceA nonbondedForceFunctionFirst;
    SwitchB switchingFunctionSecond;
    ForceB nonbondedForceFunctionSecond;
    SwitchC switchingFunctionThird;
    ForceC nonbondedForceFunctionThird;
    Real mySquaredCutoff;
  };
}

#endif /* ONEATOMPAIRTHREE_H */
