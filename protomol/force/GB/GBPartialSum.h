/* -*- c++ -*- */
#ifndef GBPARTIALSUM_H
#define GBPARTIALSUM_H


#include <protomol/topology/GenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/topology/ExclusionTable.h>
#include <protomol/config/Parameter.h>
#include <protomol/base/MathUtilities.h>
#include <protomol/parallel/Parallel.h>

#include <protomol/base/Report.h>

using namespace ProtoMol::Report;

namespace ProtoMol {
  
  class GBPartialSum {
    
  public:
    
    enum {DIST_R2 = 1};
    enum {CUTOFF = 0};
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    GBPartialSum(){}
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class GBBornRadii
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
  public:
    
    void operator()(Real &energy, Real &force,
                    Real distSquared, Real rDistSquared, const Vector3D &diff,
                    const GenericTopology *topo,
                    int atom1, int atom2, ExclusionClass excl) const {
      
      Real scaledCharge_i, scaledCharge_l;
      
      Real bornRad_i, bornRad_l;
      
      Real filGB, expterm;
      
      bornRad_i = topo->atoms[atom1].myGBSA_T->bornRad;
      bornRad_l = topo->atoms[atom2].myGBSA_T->bornRad;
      
      scaledCharge_i = topo->atoms[atom1].scaledCharge;
      scaledCharge_l = topo->atoms[atom2].scaledCharge;
      
      Real ril = sqrt(distSquared);
        
      expterm = (ril*ril)/(4.0*bornRad_i*bornRad_l);
      filGB = sqrt(ril*ril + bornRad_i*bornRad_l*exp(-expterm));
        
      topo->atoms[atom1].myGBSA_T->partialGBForceTerms += scaledCharge_i*scaledCharge_l*(1/(filGB*filGB))*0.5*(1/filGB)*exp(-expterm)*(bornRad_l + (ril*ril)/(4.0*bornRad_i));
      
      topo->atoms[atom2].myGBSA_T->partialGBForceTerms += scaledCharge_i*scaledCharge_l*(1/(filGB*filGB))*0.5*(1/filGB)*exp(-expterm)*(bornRad_i + (ril*ril)/(4.0*bornRad_l));

    }
    
    static void accumulateEnergy(ScalarStructure *energies, Real energy) {
    }
    
    static Real getEnergy(const ScalarStructure *energies) {
      return 0;
    }
    
    static void postProcess(const GenericTopology *topo, ScalarStructure *energies) {
      const unsigned int atomnumber = topo->atoms.size();
      
      for( unsigned int i=0; i<atomnumber; i++)
        topo->atoms[i].myGBSA_T->havePartialGBForceTerms = true;
    }
    
    static void parallelPostProcess(const GenericTopology *topo, ScalarStructure *energies) {
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
    
    //do parallel post process?
    static bool doParallelPostProcess() {
      return true;
    }
    
    // Parsing
    static std::string getId() {return keyword;}
    static unsigned int getParameterSize() {return 0;}
    void getParameters(std::vector<Parameter> &) const {}
    
    static GBPartialSum make(const std::vector<Value> &) {
      return GBPartialSum();
    }
    
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;
    
  };
  
}

#endif
