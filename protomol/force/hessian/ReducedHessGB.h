/*  -*- c++ -*-  */
#ifndef REDUCEDHESSGB_H
#define REDUCEDHESSGB_H

#include <protomol/type/Matrix3By3.h>
#include <protomol/topology/ExclusionTable.h>

namespace ProtoMol {
  class GenericTopology;

  class ReducedHessGB {
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // New methods of class ReducedHessGB
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
      Matrix3By3 operator()(Real distSquared,
                          const Vector3D &diff,
                          const GenericTopology *topo,
                          int atom1, int atom2,
                          int na,
                          Real soluteDielec,Real solventDielec,
                          ExclusionClass excl) const;

  private:
    Real FirstDerivativeFGB(
                          const GenericTopology *topo,
                          const int atom1, const int atom2,
                          const Real bderv_ij) const;


     Real SecondDerivativeFGB(
                          const GenericTopology *topo,
                          int atom1, int atom2,
                          Real d2Ridrij2, Real bderv_ij) const;


  };

}


#endif
