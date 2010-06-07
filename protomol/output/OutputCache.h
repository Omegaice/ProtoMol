/*  -*- c++ -*-  */
#ifndef OUTPUTCACHE_H
#define OUTPUTCACHE_H

#include <protomol/type/Real.h>
#include <protomol/type/PDB.h>
#include <protomol/type/PAR.h>
#include <protomol/type/PSF.h>
#include <protomol/type/Vector3D.h>

namespace ProtoMol {
  class ProtoMolApp;
  class Output;
  class Configuration;
  class GenericTopology;
  class ScalarStructure;
  class Vector3DBlock;
  class OutputFactory;
  class Integrator;
  struct PDB;
  struct Atom;

  /**
     OutputCache caches all kind of values, which may be needed
     by Output objects and simplifies the access to values of interest.
     Add new cached values, if needed ..
     There are some (feature) values, which will only change when the
     Topology changes
   */
  //____ OutputCache
  class OutputCache  {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    OutputCache();
    ~OutputCache();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class OutputCache
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void initialize(const ProtoMolApp *app);

    // Methods to add additional data for output objects
    void add(const std::vector<PDB::Atom> &pdbAtoms) {myAtoms = pdbAtoms;}
    void add(const PSF &psf) {myPSF = psf;}
    void add(const PAR &par) {myPAR = par;}
    void add(const std::vector<Real> &REMExchangeRate) {
      myREMRates = REMExchangeRate;
    }
    void add(const Real *replicaHistory) {myReplicaHistory = replicaHistory;}

    Real totalEnergy() const;
    Real potentialEnergy() const;
    Real kineticEnergy() const;
    Real temperature() const;
    Real volume() const;
    Real time() const;
    Real pressure() const;
    Real molecularPressure() const;
    Real molecularTemperature() const;
    Real molecularKineticEnergy() const;
    Vector3D linearMomentum() const;
    Vector3D angularMomentum() const;
    Vector3D centerOfMass() const;
    Real diffusion() const;
    Real density() const;
    Real mass() const;
    Real dihedralPhi(int index) const;
    Real brent(Real ax, Real bx, Real cx, Real tol, Real &xmin, int dihindex,
               bool max) const;

    std::vector<Real> dihedralPhis(std::vector<int> ) const;
    std::vector<std::vector<Real> > brentMaxima(std::vector<int> , bool) const;
    const Vector3DBlock *minimalPositions() const;
    const std::vector<Real> &REMRates() const {return myREMRates;}
    const Real        *replicaHistory() const {return myReplicaHistory;}

    const std::vector<PDB::Atom> &pdb() const {return myAtoms;}
    const PSF &psf() const {return myPSF;}
    const PAR &par() const {return myPAR;}

    /// To be called before every run() or finialize()
    void uncache() const;

    void setRestore() {myRestore = true;}
    void clearRestore() {myRestore = false;}
    bool restore() const {return myRestore;}

    const ProtoMolApp * getApp() const {return app;}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    const ProtoMolApp *app;

    Vector3DBlock *myInitialPositions;
    mutable Vector3DBlock * myMinimalPositions;

    // Additional data
    std::vector<PDB::Atom> myAtoms;
    std::vector<Real> myREMRates;
    const Real *myReplicaHistory;
    PSF myPSF;
    PAR myPAR;

    mutable bool myCachedKE;
    mutable Real myKE;
    mutable Real myT;

    mutable bool myCachedPE;
    mutable Real myPE;

    mutable bool myCachedV;
    mutable Real myV;

    mutable bool myCachedP;
    mutable Real myP;
    mutable bool myCachedMolP;
    mutable Real myMolP;

    mutable bool myCachedLinearMomentum;
    mutable Vector3D myLinearMomentum;

    mutable bool myCachedAngularMomentum;
    mutable Vector3D myAngularMomentum;

    mutable bool myCachedCenterOfMass;
    mutable Vector3D myCenterOfMass;

    mutable bool myCachedDiffusion;
    mutable Real myDiffusion;

    mutable bool myCachedDensity;
    mutable Real myDensity;

    mutable bool myCachedMass;
    mutable Real myMass;

    mutable int myCachedDihedralPhi;
    mutable Real myDihedralPhi;

    mutable bool myCachedDihedralPhis;
    mutable std::vector<Real> *myDihedralPhis;

    mutable bool myCachedBrentMaxima;
    mutable std::vector<std::vector<Real> > *myBrentMaxima;

    mutable bool myCachedMolT;
    mutable Real myMolT;

    mutable bool myCachedMolKE;
    mutable Real myMolKE;

    mutable bool myCachedMinimalPositions;

    bool myRestore;
  };
}
#endif
