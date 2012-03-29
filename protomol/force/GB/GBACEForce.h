/* -*- c++ -*- */
#ifndef GBACEFORCE_H
#define GBACEFORCE_H

#include <string>

#include <protomol/type/ScalarStructure.h>
#include <protomol/topology/ExclusionTable.h>
#include <protomol/topology/GenericTopology.h>

namespace ProtoMol {
  class GBACEForce {
    public:
      enum {DIST_R2 = 1};
      enum {CUTOFF = 0};
    public:
      GBACEForce();
      GBACEForce(Real s, Real r_s);
      
      void operator()(Real &energy, Real &force, Real distSquared, Real rDistSquared, 
                      const Vector3D &diff, const GenericTopology *topo, int atom1, 
                      int atom2, ExclusionClass excl) const;

      static void accumulateEnergy(ScalarStructure *energies, Real energy);
      static Real getEnergy(const ScalarStructure *energies);
      
      static void postProcess(const GenericTopology *topo, ScalarStructure *energies);
      static void parallelPostProcess(const GenericTopology *topo, ScalarStructure *energies);
      static bool doParallelPostProcess();
      
      // Parsing
      static std::string getId();
      static unsigned int getParameterSize();
      void getParameters(std::vector<Parameter> &) const;

      static GBACEForce make(const std::vector<Value> &);
    public:
      static const std::string keyword;
    private:
      Real sigma;
      Real rho_s;
  };
}

#endif
