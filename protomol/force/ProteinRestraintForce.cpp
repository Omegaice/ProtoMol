#include <protomol/force/ProteinRestraintForce.h>

using namespace std;

//____ ProteinRestraintForce
using namespace ProtoMol;
using namespace ProtoMol::Report;

ProteinRestraintForce::ProteinRestraintForce() : sphereK(0), sphereradius(0), firstEvaluate(true) {}
ProteinRestraintForce::ProteinRestraintForce(Real k, Real r, Real switchon, Real cutoff, int atom) : sphereK(k), sphereradius(r), firstEvaluate(true),
          mySwitchon(switchon), mySwitchon2(switchon * switchon), myCutoff(cutoff),
          myCutoff2(cutoff * cutoff),
          mySwitch1(1.0 / power<3>(cutoff * cutoff - switchon * switchon)),
          mySwitch2(cutoff * cutoff - 3.0 * switchon * switchon),
	  mySwitch3(4.0 / power<3>(cutoff * cutoff - switchon * switchon)),								     myAtom(atom){}

const string ProteinRestraintForce::keyword("ProteinRestraint");

unsigned int ProteinRestraintForce::getParameterSize() {
  return 5;
}

void ProteinRestraintForce::getParameters(vector<Parameter> &parameters) const {
  parameters.push_back(
                       Parameter("-sphereconstant",Value(sphereK, ConstraintValueType::NoConstraints()),3.0, Text("Sphere constant parameter")));
  parameters.push_back(
                       Parameter("-sphereradius",Value(sphereradius, ConstraintValueType::NoConstraints()), 12.0, Text("Sphere radius parameter")));
  parameters.push_back(
                       Parameter("-switchstart", Value(mySwitchon, ConstraintValueType::NotNegative()), 12.0, Text("Switch start = radius")));
  parameters.push_back(
                       Parameter("-switchend", Value(myCutoff, ConstraintValueType::Positive()), 15.0, Text("Switch end=1 > radius")));
  parameters.push_back(
		       Parameter("-atom", Value(myAtom, ConstraintValueType::Positive()), 1, Text("Atom to restrain")));
}

ProteinRestraintForce ProteinRestraintForce::make(const vector<Value> &values){
  return (ProteinRestraintForce(values[0], values[1], values[2], values[3], values[4]));
}
