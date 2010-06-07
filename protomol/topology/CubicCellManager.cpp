#include <protomol/topology/CubicCellManager.h>

#include <protomol/base/Exception.h>

using namespace std;
using namespace ProtoMol;

//____ CubicCellManager

const string CubicCellManager::keyword("Cubic");

CubicCellManager::CubicCellManager(Real r) :
  myCellSize(r), myRealCellSize(Vector3D(r, r, r)),
  myRealRCellSize(Vector3D(1.0 / r, 1.0 / r, 1.0 / r)) {}

void CubicCellManager::setCellSize(Real newSize) {
  myCellSize = newSize;
  myRealCellSize = Vector3D(newSize, newSize, newSize);
  myRealRCellSize = Vector3D(1.0 / newSize, 1.0 / newSize, 1.0 / newSize);
}

void CubicCellManager::initialize(CellListStructure &cellList,
                                  const Vector3D &min, const Vector3D &max,
                                  bool pbc) const {
  if (pbc) {
    Vector3D d(max - min);
    myRealCellSize = Vector3D(d.c[0] / std::max(1.0, floor(d.c[0] / myCellSize)),
      d.c[1] / std::max(1.0, floor(d.c[1] / myCellSize)),
      d.c[2] / std::max(1.0, floor(d.c[2] / myCellSize)));
    myRealRCellSize.c[0] = 1.0 / myRealCellSize.c[0];
    myRealRCellSize.c[1] = 1.0 / myRealCellSize.c[1];
    myRealRCellSize.c[2] = 1.0 / myRealCellSize.c[2];
    cellList.initialize(max - min, myRealCellSize);
  } else
    cellList.initialize(max - min + myRealCellSize, myRealCellSize);
}

void CubicCellManager::updateCache(CellListStructure &cellList) const {
  cellList.updateCache();
}

void CubicCellManager::getParameters(vector<Parameter> &parameters) const {
  parameters.push_back
  (Parameter("cellSize", Value(myCellSize, ConstraintValueType::Positive()),
      Text("For Periodic BC this must be < least cell basis vector."
           "  Typically 1/2 of the least cutoff value")));
}

CubicCellManager CubicCellManager::make(vector<Value> values) {
  Real r = 1.0;

  if (!values[0].get(r) || r <= 0.0)
    THROW(keyword + " cellmanager: cutoff > 0 (" + values[0].getString() +
      ").");

  return CubicCellManager(r);
}

