/* -*- c++ -*- */
#ifndef ONEATOMPAIR_H
#define ONEATOMPAIR_H

#include <protomol/topology/Topology.h>
#include <protomol/config/Parameter.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/force/OneAtomContraints.h>

namespace ProtoMol {
  //____ OneAtomPair

  /**
   * Computes the interaction for a given force between two atoms with the
   * template arguments defining the boundary conditions, switching function,
   * potential and optional constraint.
   */
  template<typename Boundary, typename Switch,
           typename Force, typename Constraint = NoConstraint>
  class OneAtomPair {
    public:
      typedef Boundary BoundaryConditions;
      typedef SemiGenericTopology<Boundary> TopologyType;
    
    public:
      OneAtomPair() : SwitchFunction(), ForceFunction() {};
      OneAtomPair(Force nF, Switch sF) :
        SwitchFunction(sF), ForceFunction(nF),
        mySquaredCutoff(Cutoff<Force::CUTOFF>::cutoff(sF, nF)) {}
    virtual ~OneAtomPair() {} // Compiler needs this

    public:
      void initialize(const TopologyType *topo, const Vector3DBlock *pos,
                      Vector3DBlock *f, ScalarStructure *e) {
        realTopo = (TopologyType *)topo;
        positions = pos;
        forces = f;
        energies = e;
      }

      void initialize(TopologyType *topo, const Vector3DBlock *pos,
                      Vector3DBlock *f, ScalarStructure *e) {
        initialize(static_cast<const TopologyType *>(topo), pos, f, e);
      }

      // Computes the force and energy for atom i and j.
      virtual void doOneAtomPair(const int i, const int j) {
        if (Constraint::PRE_CHECK)
          if (!Constraint::check(realTopo, i, j))
            return;

        // Get atom distance.
        Real distSquared;

        Vector3D diff = realTopo->
          boundaryConditions.minimalDifference((*positions)[i], (*positions)[j],
                                               distSquared);
        //      cout << "DIFF: " << diff << endl;
        // Do switching function rough test, if necessary.
        if (Switch::USE || Force::CUTOFF)
          if (distSquared > mySquaredCutoff)
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
        Real energy = 0, force = 0;
        Real rDistSquared = (Force::DIST_R2 ? 1.0 / distSquared : 1.0);
        ForceFunction(energy, force, distSquared, rDistSquared, diff,
                               realTopo, i, j, excl);
        //      cout << "EN: " << energy << " FO: " << force << endl;
        // Calculate the switched force and energy.
        if (Switch::MODIFY) {
          Real switchingValue, switchingDeriv;
          SwitchFunction(switchingValue, switchingDeriv, distSquared);
          // This has a - sign because the force is the negative of the
          // derivative of the energy (divided by the distance between the atoms).
          force = force * switchingValue - energy * switchingDeriv;
          energy = energy * switchingValue;
        }

        // Add this energy into the total system energy.
        ForceFunction.accumulateEnergy(energies, energy);
        // Add this force into the atom forces.
        Vector3D fij(diff * force);
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
          Constraint::check(realTopo, i, j, diff, energy, fij);
      }

      virtual void getParameters(std::vector<Parameter> &parameters) const {
        ForceFunction.getParameters(parameters);
        SwitchFunction.getParameters(parameters);
      }
      
      virtual void preProcess(const GenericTopology *apptopo, const Vector3DBlock *positions){
        ForceFunction.preProcess(apptopo, positions);
      }
    
      virtual void postProcess(const GenericTopology *apptopo, ScalarStructure *appenergies, Vector3DBlock *forces){
        ForceFunction.postProcess(apptopo, appenergies, forces);
      }

      virtual void parallelPostProcess(const GenericTopology *apptopo, ScalarStructure *appenergies){
        ForceFunction.parallelPostProcess(apptopo, appenergies);
      }

      virtual bool doParallelPostProcess(){
        return
          ForceFunction.doParallelPostProcess();
      }

      static OneAtomPair make(std::vector<Value> values) {
        unsigned int n = Force::getParameterSize();

        std::vector<Value> parmsNF(values.begin(), values.begin() + n);
        std::vector<Value> parmsSF(values.begin() + n, values.end());

        return OneAtomPair(Force::make(parmsNF), Switch::make(parmsSF));
      }

      static std::string getId() {
        return Constraint::getPrefixId() + Force::getId() + Constraint::getPostfixId() +
          std::string((!Switch::USE) ? std::string("") : std::string(" -switchingFunction " + Switch::getId()));
      }
    protected:
      mutable TopologyType *realTopo;
      const Vector3DBlock *positions;
      Vector3DBlock *forces;
      ScalarStructure *energies;
      Switch SwitchFunction;
      Force ForceFunction;
      Real mySquaredCutoff;
  };
}
#endif /* ONEATOMPAIR_H */
