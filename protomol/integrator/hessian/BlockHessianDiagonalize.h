/* -*- c++ -*- */
#ifndef BLOCKHESSIANDIAGONALIZE_H
#define BLOCKHESSIANDIAGONALIZE_H

#include <protomol/integrator/normal/NormalModeUtilities.h>
#include <protomol/integrator/hessian/BlockHessian.h>
#include <protomol/type/BlockMatrix.h>
#include <protomol/base/Timer.h>

namespace ProtoMol {
  /**
   *
   * Diagonalizes the Block Hessians
   * for the current force field.
   *
   */
  class BlockHessianDiagonalize {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Types and Enums
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  enum {SYMHESS = 0};
  enum {OUTPUTIHESS = 0};
  enum {OUTPUTEGVAL = 0};
  enum {OUTPUTBHESS = 0};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    BlockHessianDiagonalize();
    ~BlockHessianDiagonalize();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class BlockHessianDiagonalize
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void initialize(BlockHessian * bHessIn, const int sz);
    void initialize(const int sz);
    Real findEigenvectors(const Vector3DBlock *myPositions,
                          GenericTopology *myTopo, double * mhQu, 
                          const int _3N, const int _rfM, 
                          const Real blockCutoffDistance, 
                          const Real eigenValueThresh,
                          const int blockVectorCols);
    int diagHessian(double *eigVecO, double *eigValO,
                    double *hsnhessM, int dim, int &numFound);
    void absSort(double *eigVec, double *eigVal, int *eigIndx, int dim);

  private:
    void innerHessian();
    Real findCoarseBlockEigs(const Real eigenValueThresh, const int blockVectorCols);
    void fullElectrostaticBlocks();
    void outputDiagnostics(int typ); 

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Friends of class BlockHessianDiagonalize
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // private data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    BlockHessian * bHess;
    //Residues
    int *blocks_num_eigs;
    double *rE;
    //Residue eigenvectors as Block Matrices
    vector<BlockMatrix> blockEigVect;
    BlockMatrix innerDiag, innerEigVec;
    vector<Real> blocVectCol;
    //
  public:
    //Diagnostic data
    Timer rediagTime, hessianTime;
    //Diag data
    double *eigVal;
    int *eigIndx;
    //Blocks
    int residues_total_eigs;
    int memory_footprint;

  };
}
#endif
