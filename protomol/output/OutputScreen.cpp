#include <protomol/output/OutputScreen.h>
#include <protomol/output/OutputCache.h>
#include <protomol/integrator/Integrator.h>
#include <protomol/config/Configuration.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/module/MainModule.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/module/OutputModule.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____ OutputScreen
const string OutputScreen::keyword("Screen");

OutputScreen::OutputScreen() :
  Output(), myUnit("fs"), myFactor(1.0) {}

OutputScreen::OutputScreen(int freq) :
  Output(freq), myUnit("fs"), myFactor(1.0) {}

void OutputScreen::doInitialize() {
  Real step = app->integrator->getTimestep() *
    max(1, getOutputFreq());

  if (step >= 1e13) {
    myUnit = "s";
    myFactor = 1e-15;

  } else if (step >= 1e10) {
    myUnit = "ms";
    myFactor = 1e-12;

  } else if (step >= 1e7) {
    myUnit = "us";
    myFactor = 1e-9;

  } else if (step >= 1e4) {
    myUnit = "ns";
    myFactor = 1e-6;

  } else if (step >= 1e1) {
    myUnit = "ps";
    myFactor = 1e-3;
  }
}

void OutputScreen::doRun(int step) {
  report << plain << "Step : ";
  report.setf(ios::right);
  report.width(7);
  report << step << ", Time : ";
  report.width(10);
  report.setf(ios::showpoint | ios::fixed);
  report.precision(3);
  report << app->outputCache.time() * myFactor << " [" << myUnit << "], TE : ";
  report.precision(4);
  report.width(12);
  report << app->outputCache.totalEnergy() << " [kcal/mol]";
  report << ", T : ";
  report.precision(4);
  report.width(10);
  report << app->outputCache.temperature() << " [K]";
  report << ", V : ";
  report.precision(2);
  report.width(12);
  report << app->outputCache.volume() << " [AA^3]" << endr;
  report.reset();
}

Output *OutputScreen::doMake(const vector<Value> &values) const {
  return new OutputScreen(values[1]);
}

bool OutputScreen::isIdDefined(const Configuration *config) const {
  return config->valid("outputFreq") && !config->empty(getId()) &&
    (!config->valid(getId()) || ((*config)[getId()] == true));
}

void OutputScreen::getParameters(vector<Parameter> &parameter) const {
  parameter.push_back(Parameter(getId(), Value(true), true));
  parameter.push_back
    (Parameter(keyword + "OutputFreq",
               Value(getOutputFreq(), ConstraintValueType::Positive())));
}

bool OutputScreen::adjustWithDefaultParameters(vector<Value> &values,
                                               const Configuration *config)
const {
  if (!checkParameterTypes(values)) return false;

  if (config->valid(InputOutputfreq::keyword) && !values[1].valid())
    values[1] = (*config)[InputOutputfreq::keyword];

  return checkParameters(values);
}
