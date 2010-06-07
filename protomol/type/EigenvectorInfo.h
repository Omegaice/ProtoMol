/*  -*- c++ -*-  */
#ifndef EIGENVECTORINFO_H
#define EIGENVECTORINFO_H

#include <string>
#include <iostream>
#include <vector>

using std::cout;
using std::endl;
namespace ProtoMol
{
  //_________________________________________________________________XYZ
  /**
   * Container holding coordinates/Vector3D and names
   */
  struct EigenvectorInfo {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EigenvectorInfo() : myEigenvectorLength ( 0 ), myNumEigenvectors ( 0 ), 
        myNumUsedEigenvectors( 0 ),
        myEigenvectors ( 0 ), myOrigCEigval( 0.0 ), myNewCEigval( 0.0 ), myOrigTimestep( 0.0 ), 
        reDiagonalize( false ),
        mySingleEigs( 0 ), myEigVecChanged( true ), myMinimumLimit( 0.5 ), currentMode( -1 ) {};

    EigenvectorInfo( unsigned int n, unsigned int m ) : myEigenvectorLength( n ), 
        myNumEigenvectors( m ), myNumUsedEigenvectors( 0 ),
        myMaxEigenvalue( 0.0 ), myEigenvectors ( new double[n * m * 3] ), 
        myOrigCEigval( 0.0 ), myNewCEigval( 0.0 ), myOrigTimestep( 0.0 ), reDiagonalize( false ),
        mySingleEigs( 0 ),
        myEigVecChanged( true ), myMinimumLimit( 0.5 ), currentMode( -1 ) {}

    ~EigenvectorInfo() {
      if ( myEigenvectors ) {
        delete [] myEigenvectors;
        myEigenvectors = 0;
      }
      if ( mySingleEigs ) {
        delete [] mySingleEigs;
        mySingleEigs = 0;
      }
    };

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class EigenvectorInfo
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    bool initializeEigenvectors() {
      try {

        myEigenvectors = new double[myEigenvectorLength * myNumEigenvectors * 3];

      } catch ( std::bad_alloc& ) {

        return false;

      }

      return true;

    }

    float* getFloatEigPointer() {

      //no eigenvectors assigned?
      if(!myEigenvectors) return (float*)0;

      const unsigned int arrayLen = myEigenvectorLength * myNumEigenvectors * 3;

      //assign storage if required
      if(!mySingleEigs) {  

        try {
          mySingleEigs = new float[arrayLen];
        } catch ( std::bad_alloc& ) {
          return (float*)0;
        }

      }

      //update values if double array updated
      if(myEigVecChanged) {

        for( unsigned int i=0; i<arrayLen; i++)
          mySingleEigs[i] = (float)myEigenvectors[i];

        myEigVecChanged = false;

      }

      return mySingleEigs;

    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Eigenvector information
    unsigned int myEigenvectorLength;
    unsigned int myNumEigenvectors;
    unsigned int myNumUsedEigenvectors;
    double myMaxEigenvalue;
    double *myEigenvectors;
    
    //Current and original max subspace eigenvalue, 
    //for adaptive timestep
    double myOrigCEigval, myNewCEigval;
    double myOrigTimestep;

    //re-diagonalization flag
    bool reDiagonalize;
    
    //OpenMM single precision interface
    float *mySingleEigs;
    bool myEigVecChanged;

    double myMinimumLimit;

    //Analytic integrator
    std::vector< double > myEigenvalues;
    int currentMode;


  };
}
#endif /* EIGENVECTORINFO_H */
