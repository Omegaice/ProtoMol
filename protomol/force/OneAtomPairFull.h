/* -*- c++ -*- */
#ifndef ONEATOMPAIRFULL_H
#define ONEATOMPAIRFULL_H

#include <protomol/topology/Topology.h>
#include <protomol/config/Parameter.h>
#include <protomol/force/OneAtomContraints.h>

namespace ProtoMol {
  //____ OneAtomPairFull

  template<typename Boundary, typename Switch,
           typename Force, typename Constraint = NoConstraint>
  class OneAtomPairFull {
    // Computes the interaction for a given force between two atoms.

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Typedef
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    typedef Boundary BoundaryConditions;
    // Make the boundary conditions visible

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    OneAtomPairFull() : switchingFunction(), nonbondedForceFunction() {};
    OneAtomPairFull(Force nF, Switch sF) :
      switchingFunction(sF), nonbondedForceFunction(nF) {};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class OneAtomPairFull
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void initialize(const SemiGenericTopology<Boundary> *topo,
                    const Vector3DBlock *pos, Vector3DBlock *f,
                    ScalarStructure *e, const std::vector<Vector3D> *l) {
      realTopo = topo;
      positions = pos;
      forces = f;
      energies = e;
      lattice = l;
    }

    // Computes the force and energy for atom i and j.
    void doOneAtomPair(const int i, const int j) {
      if (Constraint::PRE_CHECK)
        if (!Constraint::check(realTopo, i, j))
          return;

      // Do we have something to do?
      bool same = (i == j);
      if (same && lattice->empty())
        return;

      Vector3D diffMinimal
        (realTopo->boundaryConditions.minimalDifference((*positions)[i],
                                                        (*positions)[j]));
      if (!same) {
        // Get atom distance.
        Real distSquared = diffMinimal.normSquared();
        // Do switching function rough test, if necessary.
        if (Switch::USE &&
            !switchingFunction.roughTest(distSquared))
          return;

        // Check for an exclusion.
        ExclusionClass excl = realTopo->exclusions.check(i, j);
        if (excl != EXCLUSION_FULL) {
          // Calculate the force and energy.
          Real rawEnergy, rawForce;
          Real rDistSquared = 1.0 / distSquared;
          nonbondedForceFunction(rawEnergy, rawForce, distSquared, rDistSquared,
                                 diffMinimal, realTopo, i, j, excl);
          // Calculate the switched force and energy.
          Real energy, force;
          if (Switch::USE) {
            Real switchingValue, switchingDeriv;
            switchingFunction(switchingValue, switchingDeriv, distSquared);
            energy = rawEnergy * switchingValue;
            // This has a - sign because the force is the negative of the
            // derivative of the energy (divided by the distance between the
            // atoms).
            force = rawForce * switchingValue - rawEnergy * switchingDeriv;
          } else {
            energy = rawEnergy;
            force = rawForce;
          }
          // Add this energy into the total system energy.
          nonbondedForceFunction.accumulateEnergy(energies, energy);
          // Add this force into the atom forces.
          Vector3D fij = -diffMinimal * force;
          (*forces)[i] += fij;
          (*forces)[j] -= fij;

          // compute the vector between molecular centers of mass
          int mi = realTopo->atoms[i].molecule;
          int mj = realTopo->atoms[j].molecule;
          if (mi != mj) {
            Vector3D molDiff =
              realTopo->boundaryConditions.minimalDifference
              (realTopo->molecules[mi].position,
               realTopo->molecules[mj].position);

            // Add to the atomic and molecular virials
            energies->addVirial(fij, -diffMinimal, -molDiff);
          } else
            energies->addVirial(fij, -diffMinimal);
          if (Constraint::POST_CHECK)
            Constraint::check(realTopo, i, j, diffMinimal, energy, fij);
        }
      }

      for (unsigned int k = 0; k < lattice->size(); k++) {
        Vector3D diff(diffMinimal + (*lattice)[k]);
        // Get atom distance.
        Real distSquared = diff.normSquared();
        // Do switching function rough test, if necessary.
        if (Switch::USE &&
            !switchingFunction.roughTest(distSquared))
          continue;

        // Calculate the force and energy.
        Real rawEnergy, rawForce;
        Real rDistSquared = 1.0 / distSquared;
        nonbondedForceFunction(rawEnergy, rawForce, distSquared, rDistSquared,
                               diff, realTopo, i, j, EXCLUSION_NONE);
        // Calculate the switched force and energy.
        Real energy, force;
        if (Switch::USE) {
          Real switchingValue, switchingDeriv;
          switchingFunction(switchingValue, switchingDeriv, distSquared);
          energy = rawEnergy * switchingValue;
          // This has a - sign because the force is the negative of the
          // derivative of the energy (divided by the distance between the
          // atoms).
          force = rawForce * switchingValue - rawEnergy * switchingDeriv;
        } else {
          energy = rawEnergy;
          force = rawForce;
        }
        // Correct the energy by factor 1/2 when same atom since
        // there is only one pair (i,j) with i==j, where
        // there are two pairs with same contribution with i !=j
        if (same) energy /= 2;
        else {
          // Add this force into the atom forces.
          Vector3D fij = -diff * force;
          (*forces)[i] += fij;
          (*forces)[j] -= fij;

          // compute the vector between molecular centers of mass
          int mi = realTopo->atoms[i].molecule;
          int mj = realTopo->atoms[j].molecule;
          if (mi != mj) {
            Vector3D molDiff =
              realTopo->boundaryConditions.minimalDifference
              (realTopo->molecules[mi].position,
               realTopo->molecules[mj].position);

            // Add to the atomic and molecular virials
            energies->addVirial(fij, -diff, -molDiff);
          } else
            energies->addVirial(fij, -diff);
        }
        // Add this energy into the total system energy.
        nonbondedForceFunction.accumulateEnergy(energies, energy);
        if (Constraint::POST_CHECK)
          Constraint::check(realTopo, i, j, diff, energy, -diff * force);
      }
    }

    void getParameters(std::vector<Parameter> &parameters) const {
      nonbondedForceFunction.getParameters(parameters);
      switchingFunction.getParameters(parameters);
    }
    
    void postProcess(const GenericTopology *apptopo, ScalarStructure *appenergies){
		  nonbondedForceFunction.postProcess(apptopo, appenergies);
	  }

    void parallelPostProcess(const GenericTopology *apptopo, ScalarStructure *appenergies){
		  nonbondedForceFunction.parallelPostProcess(apptopo, appenergies);
	  }

    bool doParallelPostProcess(){
      return
        nonbondedForceFunction.doParallelPostProcess();
	  }
    
    static unsigned int getParameterSize() {
      return Force::getParameterSize() +
        Switch::getParameterSize();
    }

    static OneAtomPairFull make(std::vector<Value> values) {
      unsigned int n = Force::getParameterSize();

      std::vector<Value> F(values.begin(), values.begin() + n);
      std::vector<Value> S(values.begin() + n, values.end());

      return OneAtomPairFull(Force::make(F),
                             Switch::make(S));
    }

    static std::string getId() {
      return Constraint::getPrefixId() + Force::getId() +
        Constraint::getPostfixId() +
        std::string((!Switch::USE) ? std::string("") :
                    std::string(" -switchingFunction " +
                                Switch::getId()));
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    const SemiGenericTopology<Boundary> *realTopo;
    const Vector3DBlock *positions;
    Vector3DBlock *forces;
    ScalarStructure *energies;
    Switch switchingFunction;
    Force nonbondedForceFunction;
    const std::vector<Vector3D> *lattice;
  };
}
#endif /* ONEATOMPAIR_H */
