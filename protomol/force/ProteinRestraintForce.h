/* -*- c++ -*- */
#ifndef PROTEINRESTRAINTFORCE_H
#define PROTEINRESTRAINTFORCE_H

#include <protomol/topology/Topology.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/config/Parameter.h>
#include <protomol/topology/TopologyUtilities.h>
#include <string>
#include <protomol/base/Report.h>

using namespace ProtoMol::Report;

namespace ProtoMol {
  //____ ProteinRestraintForce
  class ProteinRestraintForce {
  public:
    enum {DIST_R2 = 1};
    enum {CUTOFF = 0};
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class LennardJonesForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ProteinRestraintForce();
    ProteinRestraintForce(Real k, int atom);
    
    void operator()(Real &energy, Real &force,
                    Real /*distSquared*/, Real rDistSquared, const Vector3D &,
                    const GenericTopology *topo, int atom1, int atom2,
                    ExclusionClass excl)  {

      if( !forceCalculated  && myAtom == atom1){
        
        forceCalculated = true;
       
        Vector3D diff = (*pos)[atom1] - centerofmass;

        Real distance = diff.norm();
        
        //calc energy
        energy = sphereK * distance * distance;
        
        //Real myForce = 
        const Real dpotdr = 2.0 * sphereK * distance;  // Calculate dpot/dr
        
        const Real scalarForce = -dpotdr / distance;
        
        // Calculate force on atom1 due to atom2.
        Vector3D force1(diff * scalarForce);
        
        // Add to the total force.
        myForces = force1;
        
      }
    }

    static void accumulateEnergy(ScalarStructure *energies, Real energy) {
      (*energies)[ScalarStructure::OTHER] += energy;
    }

    static Real getEnergy(const ScalarStructure *energies) {
      return (*energies)[ScalarStructure::OTHER];
    }
    
    void preProcess(const GenericTopology *apptopo, const Vector3DBlock *positions) {
      
      if( firstEvaluate ){
        
        firstEvaluate = false;
        
        //get center of mass
        centerofmass = centerOfMass(positions, apptopo);        
      }

      //reset flags
      forceCalculated = false;
      
      //save position pointer
      pos = positions;
		  
    }
    
    //add in forces
    void postProcess(const GenericTopology *topo, ScalarStructure *energies, Vector3DBlock *forces){
      if(forceCalculated) (*forces)[myAtom] += myForces;
    }

    static void parallelPostProcess(const GenericTopology *topo, ScalarStructure *energies) {
      
    }

    //do parallel post process?
    static bool doParallelPostProcess() {
      return false;
    }

    // Parsing
    static std::string getId() {return keyword;}
    static unsigned int getParameterSize();
    void getParameters(std::vector<Parameter> &) const;

    static ProteinRestraintForce make(const std::vector<Value> &);
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;
  private:
    Vector3D myForces;
    Vector3D centerofmass;
    int myAtom;
    const Vector3DBlock *pos;
    bool forceCalculated;
    Real sphereK;
    
    //"first" flag
    bool firstEvaluate;
    
  
  };
  //____ INLINES
}

#endif /* PROTEINRESTRAINTFORCE_H */
