/* -*- c++ -*- */
#ifndef ONEATOMPAIRNOEXCLUSION_H
#define ONEATOMPAIRNOEXCLUSION_H

#include <protomol/topology/Topology.h>
#include <protomol/config/Parameter.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/force/OneAtomPair.h>
#include <protomol/force/OneAtomContraints.h>

namespace ProtoMol {
  template<typename Boundary, typename Switch, typename Force, typename Constraint = NoConstraint>
  class OneAtomPairNoExclusion : public OneAtomPair<Boundary,Switch,Force,Constraint> {
    typedef OneAtomPair<Boundary,Switch,Force,Constraint> Base;
    
    public:
      OneAtomPairNoExclusion() : Base() {
        
      }
      
      OneAtomPairNoExclusion(Force nF, Switch sF) : Base( nF, sF ){
        
      }
    
      void doOneAtomPair(const int i, const int j) {
        if (Constraint::PRE_CHECK){
          if (!Constraint::check(Base::realTopo, i, j)) return;
        }

        // Get atom distance.
        Real distSquared;

        Vector3D diff = Base::realTopo-> boundaryConditions.minimalDifference(
          (*Base::positions)[i], (*Base::positions)[j], distSquared
        );
        
        // Do switching function rough test, if necessary.
        if (Switch::USE || Force::CUTOFF){
          if (distSquared > Base::mySquaredCutoff) return;
        }
        
        // Don't Check for an exclusion.
        int mi = Base::realTopo->atoms[i].molecule;
        int mj = Base::realTopo->atoms[j].molecule;
        bool same = (mi == mj);
        ExclusionClass excl =
          (same ? Base::realTopo->exclusions.check(i, j) : EXCLUSION_NONE);

        // Calculate the force and energy.
        Real energy = 0, force = 0;
        Real rDistSquared = (Force::DIST_R2 ? 1.0 / distSquared : 1.0);
        Base::ForceFunction(energy, force, distSquared, rDistSquared, diff,
                               Base::realTopo, i, j, excl);
        
        // Calculate the switched force and energy.
        if (Switch::MODIFY) {
          Real switchingValue, switchingDeriv;
          Base::SwitchFunction(switchingValue, switchingDeriv, distSquared);
          // This has a - sign because the force is the negative of the
          // derivative of the energy (divided by the distance between the atoms).
          force = force * switchingValue - energy * switchingDeriv;
          energy = energy * switchingValue;
        }

        // Add this energy into the total system energy.
        Base::ForceFunction.accumulateEnergy(Base::energies, energy);
        
        // Add this force into the atom forces.
        Vector3D fij(diff * force);
        (*Base::forces)[i] -= fij;
        (*Base::forces)[j] += fij;

        // compute the vector between molecular centers of mass
        if (!same && Base::energies->molecularVirial()){
          // Add to the atomic and molecular virials
          Base::energies->
            addVirial(fij, diff, Base::realTopo->boundaryConditions.
                      minimalDifference(Base::realTopo->molecules[mi].position,
                                        Base::realTopo->molecules[mj].position));
        } else if (Base::energies->virial()) {
          Base::energies->addVirial(fij, diff);
        }
        
        // End of force computation.
        if (Constraint::POST_CHECK){
          Constraint::check(Base::realTopo, i, j, diff, energy, fij);
        }
      }
      
      virtual void preProcess(const GenericTopology *apptopo, const Vector3DBlock *positions){
        Base::ForceFunction.preProcess(apptopo, positions);
      }
      
      virtual void postProcess(const GenericTopology *apptopo, ScalarStructure *appenergies){
        Base::ForceFunction.postProcess(apptopo, appenergies);
      }
      
      virtual void parallelPostProcess(const GenericTopology *apptopo, ScalarStructure *appenergies){
        Base::ForceFunction.parallelPostProcess(apptopo, appenergies);
      }
    
      static OneAtomPairNoExclusion make(std::vector<Value> values) {
        unsigned int n = Force::getParameterSize();

        std::vector<Value> parmsNF(values.begin(), values.begin() + n);
        std::vector<Value> parmsSF(values.begin() + n, values.end());

        return OneAtomPairNoExclusion(Force::make(parmsNF), Switch::make(parmsSF));
      }

      static std::string getId() {
        return Constraint::getPrefixId() + Force::getId() + Constraint::getPostfixId() +
          std::string((!Switch::USE) ? std::string("") : std::string(" -switchingFunction " + Switch::getId()));
      }
  };
}
#endif /* ONEATOMPAIRNOEXCLUSION_H */
