#include <protomol/output/OutputCache.h>
#include <protomol/output/Output.h>
#include <protomol/config/Configuration.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/integrator/Integrator.h>

#include <protomol/base/Report.h>
#include <protomol/base/MathUtilities.h>
#include <protomol/base/SystemUtilities.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>

#include <protomol/ProtoMolApp.h>

#include <protomol/base/Zap.h>
#include <protomol/base/Exception.h>

#include <limits>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____ OutputCache
OutputCache::OutputCache() :
  myInitialPositions(new Vector3DBlock()),
  myMinimalPositions(new Vector3DBlock()),
  myDihedralPhi(Constant::REAL_NAN),
  myDihedralPhis(new vector<Real>()),
  myBrentMaxima(new vector<vector<Real> >()) {
  uncache();
}

OutputCache::~OutputCache() {
  zap(myInitialPositions);
  zap(myMinimalPositions);
  zap(myDihedralPhis);
  zap(myBrentMaxima);
}

void OutputCache::initialize(const ProtoMolApp *app) {
  this->app = app;
  *myInitialPositions = app->positions;
}

Real OutputCache::totalEnergy() const {
  return potentialEnergy() + kineticEnergy();
}

Real OutputCache::potentialEnergy() const {
  if (!myCachedPE) {
    myPE = app->energies.potentialEnergy();
    myCachedPE = true;
  }
  return myPE;
}

Real OutputCache::kineticEnergy() const {
  if (!myCachedKE) {
    myKE = ProtoMol::kineticEnergy(app->topology, &app->velocities);
    myT = ProtoMol::temperature(myKE, app->topology->degreesOfFreedom);
    myCachedKE = true;
  }
  return myKE;
}

Real OutputCache::temperature() const {
  if (!myCachedKE) {
    myKE = ProtoMol::kineticEnergy(app->topology, &app->velocities);
    myT = ProtoMol::temperature(myKE, app->topology->degreesOfFreedom);
    myCachedKE = true;
  }
  return myT;
}

Real OutputCache::molecularTemperature() const {
  if (!myCachedMolT) {
    myMolT = ProtoMol::temperature(
      molecularKineticEnergy(), 3 * app->topology->molecules.size());
    myCachedMolT = true;
  }
  return myMolT;
}

Real OutputCache::molecularKineticEnergy() const {
  if (!myCachedMolKE) {
    myMolKE = ProtoMol::molecularKineticEnergy(app->topology, &app->velocities);
    myCachedMolKE = true;
  }
  return myMolKE;
}

Real OutputCache::pressure() const {
  if (!myCachedP) {
    if (!app->energies.virial())
      myP = 0.0;
    else if (volume() > 0.0)
      myP = ProtoMol::computePressure(&app->energies, volume(),
                                      kineticEnergy());
    else
      myP = Constant::REAL_INFINITY;
    myCachedP = true;
  }
  return myP;
}

Real OutputCache::molecularPressure() const {
  if (!myCachedMolP) {
    if (!app->energies.molecularVirial())
      myMolP = 0.0;
    else if (volume() > 0.0)
      myMolP = ProtoMol::computeMolecularPressure(&app->energies,
                                                  volume(),
                                                  molecularKineticEnergy());
    else
      myMolP = Constant::REAL_INFINITY;
    myCachedMolP = true;
  }
  return myMolP;
}

Real OutputCache::volume() const {
  if (!myCachedV) {
    myV = app->topology->getVolume(app->positions);
    myCachedV = true;
  }
  return myV;
}

Vector3D OutputCache::linearMomentum() const {
  if (!myCachedLinearMomentum) {
    myLinearMomentum =
      ProtoMol::linearMomentum(&app->velocities, app->topology);
    myCachedLinearMomentum = true;
  }
  return myLinearMomentum;
}

Vector3D OutputCache::angularMomentum() const {
  if (!myCachedAngularMomentum) {
    myAngularMomentum =
      ProtoMol::angularMomentum(&app->positions, &app->velocities,
                                app->topology, OutputCache::centerOfMass());
    myCachedAngularMomentum = true;
  }
  return myAngularMomentum;
}

Vector3D OutputCache::centerOfMass() const {
  if (!myCachedCenterOfMass) {
    myCenterOfMass = ProtoMol::centerOfMass(&app->positions, app->topology);
    myCachedCenterOfMass = true;
  }
  return myCenterOfMass;
}

Real OutputCache::diffusion() const {
  if (!myCachedDiffusion) {
    myDiffusion = 0.0;
    unsigned int numberOfAtoms = app->positions.size();
    for (unsigned int i = 0; i < numberOfAtoms; i++)
      myDiffusion +=
        (app->positions[i] - (*myInitialPositions)[i]).normSquared();

    myDiffusion /= (6.0 * numberOfAtoms);
    myCachedDiffusion = true;
  }
  return myDiffusion;
}

Real OutputCache::density() const {
  if (!myCachedDensity) {
    myDensity =
      (volume() > 0.0 ?
       (mass() / volume() * Constant::SI::AMU *
        power<3>(Constant::SI::LENGTH_AA) *
        1e-3) :
       Constant::REAL_NAN);
    myCachedDensity = true;
  }
  return myDensity;
}

Real OutputCache::mass() const {
  if (!myCachedMass) {
    myMass = 0.0;
    unsigned int numberOfAtoms = app->positions.size();
    for (unsigned int i = 0; i < numberOfAtoms; i++)
      myMass += app->topology->atoms[i].scaledMass;

    myCachedMass = true;
  }
  return myMass;
}

Real OutputCache::time() const {
  return app->topology->time;
}

const Vector3DBlock *OutputCache::minimalPositions() const {
  if (!myCachedMinimalPositions) {
    *myMinimalPositions = app->positions;
    (const_cast<GenericTopology *>(app->topology))->
      minimalImage(*myMinimalPositions);
  }
  myCachedMinimalPositions = true;
  return myMinimalPositions;
}

Real OutputCache::dihedralPhi(int index) const {
  if (index < 0 || index >= static_cast<int>(app->topology->dihedrals.size()))
    index = -1;

  if (index < 0) {
    myDihedralPhi = Constant::REAL_NAN;
    myCachedDihedralPhi = index;
    return myDihedralPhi;
  }

  myDihedralPhi = computePhiDihedral(app->topology, &app->positions, index);
  return myDihedralPhi;
}

vector<Real> OutputCache::dihedralPhis(vector<int> dihedralset) const {
  if (!myCachedDihedralPhis) {
    myDihedralPhis->resize(dihedralset.size());
    for (unsigned int i = 0; i < dihedralset.size(); ++i)
      (*myDihedralPhis)[i] =
        computePhiDihedral(app->topology, &app->positions, dihedralset[i]);

    // different functions require different dihedralset
    myCachedDihedralPhis = false;
  }
  return *myDihedralPhis;
}

// Brent's Maxima function goes here to use topology for dihedral well
// calculation
vector<vector<Real> > OutputCache::brentMaxima(vector<int> dihedralset,
                                               bool max) const {
  if (!myCachedBrentMaxima) {
    //The Brent algorithm gets the maxima if maxmin = -1
    int maxmin = -1;
    if (!max)
      maxmin = 1;

    myBrentMaxima->clear();
    myBrentMaxima->resize(dihedralset.size());

    for (unsigned int i = 0; i < dihedralset.size(); ++i) {
      //note the function evaluates one step past 2 pi
      for (unsigned int j = 0; j <= 99; j++) {
        Real lradangle = (M_PI * 2 / 100 * j);
        Real radangle = (M_PI * 2 / 100 * (j + 1));
        Real rradangle = (M_PI * 2 / 100 * (j + 2));

        Real valLangle = maxmin *
                         computePhiDihedralEnergy(app->topology, dihedralset[i],
                                                  lradangle);

        Real valRangle = maxmin *
                         computePhiDihedralEnergy(app->topology, dihedralset[i],
                                                  rradangle);

        Real valAngle = maxmin *
                        computePhiDihedralEnergy(app->topology, dihedralset[i],
                                                 radangle);

        Real xmax = 0.0;

        Real tol = 0.01;
        if ((valLangle > valAngle) && (valRangle > valAngle)) {
          brent(lradangle, radangle, rradangle, tol, xmax, dihedralset[i],
                max);

          ((*myBrentMaxima)[i]).push_back(xmax);
        }
      }

      // Throws Warning if no maxima were found
      if (((*myBrentMaxima)[i]).size() == 0) {
        report << warning
               << "No dihedral maxima found for dihedral index: "
               << dihedralset[i] << " Check dihedral energy equation"
               << endr;
        ((*myBrentMaxima)[i]).push_back(0.0);
      }
    }

    //temp hack that allows multiple calls but defeats the purpose of
    // the app.outputCache... please fix me!
    //myCachedBrentMaxima= true;
  }
  return *myBrentMaxima;
}

//BRENT FUNCTION
Real OutputCache::brent(Real ax, Real bx, Real cx, Real tol, Real &xmin,
                        int dihindex,
                        bool max) const {
  const int ITMAX = 100;
  const Real CGOLD = 0.3819660;
  const Real ZEPS = numeric_limits<Real>::epsilon() * 1.0e-3;
  Real a, b, d = 0.0, etemp, fu, fv, fw, fx;
  Real p, q, r, tol1, tol2, u, v, w, x, xm;
  Real e = 0.0;

  //The Brent algorithm gets the maxima if maxmin = -1
  int maxmin = -1;
  if (!max)
    maxmin = 1;

  a = (ax < cx ? ax : cx);
  b = (ax > cx ? ax : cx);
  x = w = v = bx;
  fw = fv = fx = maxmin * computePhiDihedralEnergy(app->topology, dihindex, x);
  for (int iter = 0; iter < ITMAX; iter++) {
    xm = 0.5 * (a + b);
    tol2 = 2.0 * (tol1 = tol * fabs(x) + ZEPS);
    if (fabs(x - xm) <= (tol2 - 0.5 * (b - a))) {
      xmin = x;
      return fx;
    }
    if (fabs(e) > tol1) {
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);
      if (q < 0.0) p = -p;
      q = fabs(q);
      etemp = e;
      e = d;
      if (fabs(p) >=
          fabs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x))
        d = CGOLD * (e = (x >= xm ? a - x : b - x));
      else {
        d = p / q;
        u = x + d;
        if (u - a < tol2 || b - u < tol2)
          d = sign(tol1, xm - x);
      }
    } else
      d = CGOLD * (e = (x >= xm ? a - x : b - x));
    u = (fabs(d) >= tol1 ? x + d : x + sign(tol1, d));
    fu = maxmin * computePhiDihedralEnergy(app->topology, dihindex, u);
    if (fu <= fx) {
      if (u >= x) a = x;else b = x;
      shift(v, w, x, u);
      shift(fv, fw, fx, fu);
    } else {
      if (u < x) a = u;else b = u;
      if (fu <= fw || w == x) {
        v = w;
        w = u;
        fv = fw;
        fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u;
        fv = fu;
      }
    }
  }

  //too many iterations in brent
  xmin = x;
  return fx;
}

void OutputCache::uncache() const {
  myCachedKE = false;
  myCachedPE = false;
  myCachedV = false;
  myCachedP = false;
  myCachedMolP = false;
  myCachedLinearMomentum = false;
  myCachedAngularMomentum = false;
  myCachedCenterOfMass = false;
  myCachedDiffusion = false;
  myCachedDensity = false;
  myCachedMass = false;
  myCachedDihedralPhis = false;
  myCachedDihedralPhi = -1;
  myCachedBrentMaxima = false;
  myCachedMolT = false;
  myCachedMolKE = false;
  myCachedMinimalPositions = false;
}
