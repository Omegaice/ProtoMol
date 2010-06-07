/* -*- c++ -*- */
#ifndef CUBICCELLMANAGER_H
#define CUBICCELLMANAGER_H

#include <protomol/topology/ArrayCellListStructure.h>
#include <protomol/config/Parameter.h>
#include <protomol/type/Vector3D.h>

#include <vector>

namespace ProtoMol {
  //________________________________________ CubicCellManager
  /**
   * The cell manager for equal-sized (cubic) cells. For optimization reasons
   * in case of periodic boundary conditions the cells are not cubic any more
   * in order to fit the system by multiples of the cell dimensions.
   */
  class CubicCellManager {
  public:
    /// topology and cell location structure of the cell
    typedef CubicCellLocation Cell;
    /// implementation of the cell list
    typedef ArrayCellListStructure CellListStructure;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    CubicCellManager() : myCellSize(0.0) {}
    CubicCellManager(Real r);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class CubicCellManager
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    /// Set the size of each cell.
    void setCellSize(Real newSize);

    Real getCellSize(void) const {return myCellSize;}
    Vector3D getRealCellSize(void) const {return myRealCellSize;}
    /// Get the volume of the cell.
    Real getCellVolume(void) const {return myRealCellSize.c[0] *
                                           myRealCellSize.c[1] * myRealCellSize.c[2];}

    /// Find the cell that one atom belongs to.
    Cell findCell(const Vector3D &position) const {
      return Cell((int)floor(position.c[0] * myRealRCellSize.c[0]),
        (int)floor(position.c[1] * myRealRCellSize.c[1]),
        (int)floor(position.c[2] * myRealRCellSize.c[2]));
    }

    void initialize(CellListStructure &cellList,
                    const Vector3D &min,
                    const Vector3D &max,
                    bool pbc) const;
    void updateCache(CellListStructure &cellList) const;

    const std::string &getKeyword() const {return keyword;}

    void getParameters(std::vector<Parameter> &parameters) const;
    static CubicCellManager make(std::vector<Value> values);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;
  private:
    Real myCellSize;
    mutable Vector3D myRealCellSize;
    mutable Vector3D myRealRCellSize;
  };

  //________________________________________ INLINES
}
#endif /* CUBICCELLMANAGER_H */
