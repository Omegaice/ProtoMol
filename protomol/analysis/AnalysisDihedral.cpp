#include <protomol/analysis/AnalysisDihedral.h>

#include <protomol/topology/TopologyUtilities.h>
#include <protomol/module/MainModule.h>
#include <protomol/ProtoMolApp.h>

#include <protomol/base/MathUtilities.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/base/SystemUtilities.h>

#include <protomol/io/XYZWriter.h>
#include <protomol/io/CheckpointConfigWriter.h>

#include <sstream>
#include <iostream>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

const string AnalysisDihedral::keyword("AnalyzeDihedral");

AnalysisDihedral::AnalysisDihedral()
	: Analysis(true), mPsiMin(0), mPsiMax(0), mPhiMin(0), mPhiMax(0) {}

AnalysisDihedral::AnalysisDihedral(const std::vector<int> index, const Real psiMin, const Real psiMax, const Real phiMin, const Real phiMax)
	: Analysis(true), mIndex(index), mPsiMin(psiMin), mPsiMax(psiMax), mPhiMin(phiMin), mPhiMax(phiMax) {}

void AnalysisDihedral::doInitialize() {
	int _N = app->topology->atoms.size();

	std::vector<int> residue_id;

	int current_id = app->topology->atoms[0].residue_seq;
	residue_id.push_back(current_id);

	for( int j = 0; j < _N; j++ ) {
		if( app->topology->atoms[j].residue_seq != current_id ) {
			current_id = app->topology->atoms[j].residue_seq;
			residue_id.push_back(current_id);
		}
	}

	const int num_residues = residue_id.size();

	residues_alpha_c.resize(num_residues);
	residues_phi_n.resize(num_residues);
	residues_psi_c.resize(num_residues);

	std::vector<int> atom_residue(_N);

	int last_atom = 0;
	for( int i = 0; i < num_residues; i++ ) {
		residues_alpha_c[i] = residues_phi_n[i] = residues_psi_c[i] = -1;

		int res_idx = 0;
		bool in_block = false;
		for( int j = last_atom; j < _N; j++ ) {
			if( app->topology->atoms[j].residue_seq == residue_id[i] ) {
				in_block = true;//only allow contigous atoms in 1 residue, numbers may be re-used.

				if( app->topology->atoms[j].name.compare("CA") == 0 ) {
					residues_alpha_c[i] = j;//flag valid and mark alpha carbons
				}

				atom_residue[j] = i;

				res_idx++;
			} else {
				last_atom = j;
				if( in_block ) {
					break;		//out if was in block
				}
			}
		}
	}

	const int bonds_size = app->topology->bonds.size();
	for( int i = 0; i < bonds_size; i++ ) {
		const int atom1 = app->topology->bonds[i].atom1;
		const int atom2 = app->topology->bonds[i].atom2;

		if( app->topology->atoms[atom1].name.compare("CA") == 0 ) {
			if( app->topology->atoms[atom2].name.compare("N") == 0 ) {
				residues_phi_n[atom_residue[atom1]] = atom2;
			}
			if( app->topology->atoms[atom2].name.compare("C") == 0 ) {
				residues_psi_c[atom_residue[atom1]] = atom2;
			}
		} else {
			if( app->topology->atoms[atom2].name.compare("CA") == 0 ) {
				if( app->topology->atoms[atom1].name.compare("N") == 0 ) {
					residues_phi_n[atom_residue[atom2]] = app->topology->bonds[i].atom1;
				}
				if( app->topology->atoms[atom1].name.compare("C") == 0 ) {
					residues_psi_c[atom_residue[atom2]] = atom1;
				}
			}
		}
	}
}

void AnalysisDihedral::doIt(long step) {}

const Real AnalysisDihedral::computeDihedral(const int a1, const int a2, const int a3, const int a4) const {
	Vector3D r12 = app->topology->minimalDifference(( app->positions )[a2], ( app->positions )[a1]);
	Vector3D r23 = app->topology->minimalDifference(( app->positions )[a3], ( app->positions )[a2]);
	Vector3D r34 = app->topology->minimalDifference(( app->positions )[a4], ( app->positions )[a3]);

	Vector3D a = r12.cross(r23);
	Vector3D b = r23.cross(r34);
	Vector3D c = r23.cross(a);

	Real cosPhi = a.dot(b) / ( a.norm() * b.norm());
	Real sinPhi = c.dot(b) / ( c.norm() * b.norm());

	return -atan2(sinPhi, cosPhi);
}

void AnalysisDihedral::doRun(long step) {
	std::vector<Real> phiAngle(mIndex.size()), psiAngle(mIndex.size());

	for( int i = 0; i < mIndex.size(); i++ ) {
		const int alpha_c = residues_alpha_c[mIndex[i]];
		const int phi_n = residues_phi_n[mIndex[i]];
		const int psi_c = residues_psi_c[mIndex[i]];

		for( int dihedral = 0; dihedral < app->topology->rb_dihedrals.size(); dihedral++ ) {
			const int a1 = app->topology->rb_dihedrals[dihedral].atom1;
			const int a2 = app->topology->rb_dihedrals[dihedral].atom2;
			const int a3 = app->topology->rb_dihedrals[dihedral].atom3;
			const int a4 = app->topology->rb_dihedrals[dihedral].atom4;

			if( a1 == phi_n && a2 == alpha_c && a3 == psi_c ) {
				phiAngle[i] = rtod(computeDihedral(a1, a2, a3, a4));
			}

			if( a2 == phi_n && a3 == alpha_c && a4 == psi_c ) {
				psiAngle[i] = rtod(computeDihedral(a1, a2, a3, a4));
			}
		}
	}

	bool isFolded = true;
	for( int i = 0; i < mIndex.size(); i++ ) {
		if( ( psiAngle[i] < mPsiMin && psiAngle[i] > mPsiMax ) || ( phiAngle[i] < mPhiMin && phiAngle[i] > mPhiMax ) ) {
			isFolded = false;
			break
		}
	}

	if( isFolded ) {
		report << "[AnalyzeDihedral] Molecule satisfied dihedral angle analysis for residues ";

		for( int j = 0; j < mIndex.size(); j++ ) {
			report << mIndex[j] << "(Phi:" << phiAngle << "°, Psi:" << psiAngle << "°)"
			if( j < mIndex.size() - 1 ) {
				report << ", ";
			}
		}

		report << endr;

		mShouldStop = true;
	}
}

Analysis *AnalysisDihedral::doMake(const vector<Value> &values) const {
	std::istringstream sStream;

	// Parse Indices
	sStream.str(values[1]);
	sStream.clear();

	std::vector<int> index;
	while( !sStream.eof()) {
		int value = -1;
		sStream >> value;
		sStream.ignore(std::numeric_limits<std::streamsize>::max(), ',');

		if( value != -1 ) {
			index.push_back(value);
		}
	}

	// Parse Psi
	Real psiMin, psiMax;

	sStream.str(values[2]);
	sStream >> psiMin;
	sStream.ignore(std::numeric_limits<std::streamsize>::max(), ',');
	sStream >> psiMax;

	// Parse Phi
	Real phiMin, phiMax;

	sStream.str(values[3]);
	sStream >> phiMin;
	sStream.ignore(std::numeric_limits<std::streamsize>::max(), ',');
	sStream >> phiMax;

	return new AnalysisDihedral(index, psiMin, psiMax, phiMin, phiMax);
}

bool AnalysisDihedral::isIdDefined(const Configuration *config) const {
	return config->valid(getId()) && (( *config )[getId()] == true );
}

void AnalysisDihedral::getParameters(vector<Parameter> &parameter) const {
	parameter.push_back( Parameter( getId(), Value( true ), true ) );
	parameter.push_back( Parameter( getId() + "Index", Value("", ConstraintValueType::NoConstraints()), "") );
	parameter.push_back( Parameter( getId() + "PsiRange", Value("", ConstraintValueType::NoConstraints()), "" ) );
	parameter.push_back( Parameter( getId() + "PhiRange", Value("", ConstraintValueType::NoConstraints()), "" ) );
}

bool AnalysisDihedral::adjustWithDefaultParameters(vector<Value> &values, const Configuration *config) const {
	if( !checkParameterTypes(values)) { return false; }

	if( !values[1].valid()) { values[1] = ""; }
	if( !values[2].valid()) { values[2] = ""; }
	if( !values[3].valid()) { values[3] = ""; }

	return checkParameters(values);
}
