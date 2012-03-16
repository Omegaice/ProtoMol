/* -*- c++ -*- */
#ifndef ONEATOMPAIRFULL_H
#define ONEATOMPAIRFULL_H

#include <protomol/topology/Topology.h>
#include <protomol/config/Parameter.h>
#include <protomol/force/OneAtomPair.h>
#include <protomol/force/OneAtomContraints.h>

namespace ProtoMol {
  template<typename Boundary, typename Switch, typename Force, typename Constraint = NoConstraint>
  class OneAtomPairFull : public OneAtomPair<Boundary,Switch,Force,Constraint> {
    typedef OneAtomPair<Boundary,Switch,Force,Constraint> Base;
    
    public:
      OneAtomPairFull() : Base() {
      
      }
      
      OneAtomPairFull(Force nF, Switch sF) : Base( nF, sF ){
        
      }
      
    public:
      void initialize(const SemiGenericTopology<Boundary> *topo,
                      const Vector3DBlock *pos, Vector3DBlock *f,
                      ScalarStructure *e, const std::vector<Vector3D> *l) {
        Base::realTopo = (SemiGenericTopology<Boundary> *)topo;
        Base::positions = pos;
        Base::forces = f;
        Base::energies = e;
        lattice = l;
      }

      // Computes the force and energy for atom i and j.
      void doOneAtomPair(const int i, const int j) {
        if (Constraint::PRE_CHECK)
          if (!Constraint::check(Base::realTopo, i, j))
            return;

        // Do we have something to do?
        bool same = (i == j);
        if (same && lattice->empty())
          return;

        Vector3D diffMinimal
          (Base::realTopo->boundaryConditions.minimalDifference((*Base::positions)[i],
                                                          (*Base::positions)[j]));
        if (!same) {
          // Get atom distance.
          Real distSquared = diffMinimal.normSquared();
          // Do switching function rough test, if necessary.
          if (Switch::USE &&
              !Base::SwitchFunction.roughTest(distSquared))
            return;

          // Check for an exclusion.
          ExclusionClass excl = Base::realTopo->exclusions.check(i, j);
          if (excl != EXCLUSION_FULL) {
            // Calculate the force and energy.
            Real rawEnergy, rawForce;
            Real rDistSquared = 1.0 / distSquared;
            ForceFunction(rawEnergy, rawForce, distSquared, rDistSquared,
                                   diffMinimal, Base::realTopo, i, j, excl);
            // Calculate the switched force and energy.
            Real energy, force;
            if (Switch::USE) {
              Real switchingValue, switchingDeriv;
              Base::SwitchFunction(switchingValue, switchingDeriv, distSquared);
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
            Base::ForceFunction.accumulateEnergy(Base::energies, energy);
            
            // Add this force into the atom forces.
            Vector3D fij = -diffMinimal * force;
            (*Base::forces)[i] += fij;
            (*Base::forces)[j] -= fij;

            // compute the vector between molecular centers of mass
            int mi = Base::realTopo->atoms[i].molecule;
            int mj = Base::realTopo->atoms[j].molecule;
            if (mi != mj) {
              Vector3D molDiff =
                Base::realTopo->boundaryConditions.minimalDifference
                (Base::realTopo->molecules[mi].position,
                 Base::realTopo->molecules[mj].position);

              // Add to the atomic and molecular virials
              Base::energies->addVirial(fij, -diffMinimal, -molDiff);
            } else
              Base::energies->addVirial(fij, -diffMinimal);
            if (Constraint::POST_CHECK)
              Constraint::check(Base::realTopo, i, j, diffMinimal, energy, fij);
          }
        }

        for (unsigned int k = 0; k < lattice->size(); k++) {
          Vector3D diff(diffMinimal + (*lattice)[k]);
          // Get atom distance.
          Real distSquared = diff.normSquared();
          // Do switching function rough test, if necessary.
          if (Switch::USE &&
              !Base::SwitchFunction.roughTest(distSquared))
            continue;

          // Calculate the force and energy.
          Real rawEnergy, rawForce;
          Real rDistSquared = 1.0 / distSquared;
          Base::ForceFunction(rawEnergy, rawForce, distSquared, rDistSquared,
                                 diff, Base::realTopo, i, j, EXCLUSION_NONE);
          // Calculate the switched force and energy.
          Real energy, force;
          if (Switch::USE) {
            Real switchingValue, switchingDeriv;
            Base::SwitchFunction(switchingValue, switchingDeriv, distSquared);
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
            (*Base::forces)[i] += fij;
            (*Base::forces)[j] -= fij;

            // compute the vector between molecular centers of mass
            int mi = Base::realTopo->atoms[i].molecule;
            int mj = Base::realTopo->atoms[j].molecule;
            if (mi != mj) {
              Vector3D molDiff =
                Base::realTopo->boundaryConditions.minimalDifference
                (Base::realTopo->molecules[mi].position,
                 Base::realTopo->molecules[mj].position);

              // Add to the atomic and molecular virials
              Base::energies->addVirial(fij, -diff, -molDiff);
            } else
              Base::energies->addVirial(fij, -diff);
          }
          // Add this energy into the total system energy.
          Base::ForceFunction.accumulateEnergy(Base::energies, energy);
          if (Constraint::POST_CHECK)
            Constraint::check(Base::realTopo, i, j, diff, energy, -diff * force);
        }
      }
    
      static unsigned int getParameterSize() {
        return Force::getParameterSize() + Switch::getParameterSize();
      }

      static OneAtomPairFull make(std::vector<Value> values) {
        unsigned int n = Force::getParameterSize();

        std::vector<Value> F(values.begin(), values.begin() + n);
        std::vector<Value> S(values.begin() + n, values.end());

        return OneAtomPairFull(Force::make(F), Switch::make(S));
      }

      static std::string getId() {
        return Constraint::getPrefixId() + Force::getId() + Constraint::getPostfixId() +
          std::string((!Switch::USE) ? std::string("") : std::string(" -switchingFunction " + Switch::getId()));
      }
    private:
      const std::vector<Vector3D> *lattice;
  };
}
#endif /* ONEATOMPAIR_H */
