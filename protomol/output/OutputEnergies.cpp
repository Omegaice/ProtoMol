#include <protomol/output/OutputEnergies.h>
#include <protomol/config/Configuration.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/output/OutputCache.h>
#include <protomol/module/MainModule.h>
#include <protomol/ProtoMolApp.h>

#include <iomanip>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____ Output
const string OutputEnergies::keyword("allEnergiesFile");

OutputEnergies::OutputEnergies() :
  OutputFile(), myDoMolecularTemperature(false) {}

OutputEnergies::OutputEnergies(const string &filename, int freq,
                               int cacheFreq, int cacheSize, Real closeTime,
                               bool doMolTemp) :
  OutputFile(filename, freq, cacheFreq, cacheSize,
             closeTime), myDoMolecularTemperature(doMolTemp) {}

void OutputEnergies::doInitialize() {
#ifdef HAVE_LIBFAH
    FAH::ChecksummedFile allEnergiesHeaderFile;
#else
    ofstream allEnergiesHeaderFile;
#endif

  allEnergiesHeaderFile.open(string(myFilename + ".header").c_str(),
                             ios::out | ios::trunc);
  if (!allEnergiesHeaderFile)
    report << error << " Can not open \'" << myFilename <<
    ".header\' for " << getId() << "." << endr;

  allEnergiesHeaderFile
    << setw(14) << "Time(fs)" << " "
    << setw(14) << "E_potential" << " "
    << setw(14) << "E_kinetic" << " "
    << setw(14) << "E_total" << " "
    << setw(14) << "Temperature" << " "
    << setw(14) << "E_bond" << " "
    << setw(14) << "E_angle" << " "
    << setw(14) << "E_dihedral" << " "
    << setw(14) << "E_improper" << " "
    << setw(14) << "E_VdW" << " "
    << setw(14) << "E_coulomb" << " "
    << setw(14) << "E_other" << " "
    << setw(14) << "Volume(A^3)";

  if (app->energies.virial())
    allEnergiesHeaderFile 
      << " " << setw(14) << "Pressure(bar)";
  if (app->energies.molecularVirial())
    allEnergiesHeaderFile
      << " " << setw(14) << "Mol_Pres(bar)";
  if (myDoMolecularTemperature)
    allEnergiesHeaderFile
      << " " << setw(14) << "Mol_Temp(K)";
  allEnergiesHeaderFile
      << " " << setw(20) << "E_shadow";

  allEnergiesHeaderFile << endl;
  allEnergiesHeaderFile.close();
  open();
  close();
}

void OutputEnergies::doRunCached(int) {
  myBuffer << resetiosflags(
    ios::showpoint | ios::fixed | ios::floatfield)
           << setw(14)
           << setprecision(2)
           << setiosflags(ios::showpoint | ios::fixed)
           << app->outputCache.time() << " "
           << resetiosflags(
    ios::showpoint | ios::fixed | ios::floatfield)
           << setiosflags(ios::floatfield)
           << setprecision(8)
           << setw(14)
           << app->outputCache.potentialEnergy() << " "
           << setw(14)
           << app->outputCache.kineticEnergy() << " "
           << setw(14)
           << app->outputCache.totalEnergy() << " "
           << setw(14)
           << app->outputCache.temperature() << " "
           << setw(14)
           << app->energies[ScalarStructure::BOND] << " "
           << setw(14)
           << app->energies[ScalarStructure::ANGLE] << " "
           << setw(14)
           << app->energies[ScalarStructure::DIHEDRAL] << " "
           << setw(14)
           << app->energies[ScalarStructure::IMPROPER] << " "
           << setw(14)
           << app->energies[ScalarStructure::LENNARDJONES] << " "
           << setw(14)
           << app->energies[ScalarStructure::COULOMB] << " "
           << setw(14)
           << app->energies[ScalarStructure::OTHER] << " "
           << setw(14)
           << app->outputCache.volume();
  if (app->energies.virial())
    myBuffer << " " << setw(14)
             << app->outputCache.pressure();
  if (app->energies.molecularVirial())
    myBuffer << " " << setw(14)
             << app->outputCache.molecularPressure();
  if (myDoMolecularTemperature)
    myBuffer << " " << setw(14)
             << app->outputCache.molecularTemperature();
    myBuffer << " " << setw(20)
             << setprecision(16)    //  High precision needed.
             << app->energies[ScalarStructure::SHADOW];

  myBuffer << endl;
}

Output *OutputEnergies::doMake(const vector<Value> &values) const {
  return new OutputEnergies(values[0], values[1], values[2], values[3],
                            values[4], values[5]);
}

void OutputEnergies::getParameters(vector<Parameter> &parameter) const {
  OutputFile::getParameters(parameter);
  parameter.push_back(Parameter("molecularTemperature",
                                Value(myDoMolecularTemperature), false));
}
