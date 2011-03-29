//
//     CTRUNCATE
//     Copyright (C) 2006-2011 Norman Stein, Charles Ballard
//
//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Application.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//

#include "ctruncate_matthews.h"

namespace ctruncate {

	//                       A   R   N   D   C   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
	const char Matthews::_ResidueName1[21] = {'A','R','N','D','C','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'};
	const int Matthews::_Catoms[21] = { 3,  6,  4,  4,  3,  3,  5,  5,  2,  6,  6,  6,  6,  5,  9,  5,  3,  4, 11,  9,  5 };
	const int Matthews::_Hatoms[21] = { 7, 14,  8,  7,  7,  7, 10,  9,  5,  9, 13, 13, 14, 11, 11,  9,  7,  9, 12, 11, 11 };
	const int Matthews::_Natoms[21] = { 1,  4,  2,  1,  1,  1,  2,  1,  1,  3,  1,  1,  2,  1,  1,  1,  1,  1,  2,  1,  1 };
	const int Matthews::_Oatoms[21] = { 2,  2,  3,  4,  2,  2,  3,  4,  2,  2,  2,  2,  2,  2,  2,  2,  3,  3,  2,  3,  2 };
	const int Matthews::_Satoms[21] = { 0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0 };
	const double Matthews::_rbin[15] = {1.199,1.501,1.650,1.801,1.901,2.001,2.201,2.401, 
		2.601,2.801,3.101,3.501,5.001,5.001,5.001};
	const double Matthews::_p0[15] = {0.085,0.312,0.400,0.503,0.597,0.729,1.052,1.781,
		2.852,3.386,3.841,4.281,4.592,1.503,0.257};
	const double Matthews::_vmbar[15] = {2.052,2.102,2.122,2.132,2.140,2.155,2.171,2.182,
		2.191,2.192,2.205,2.211,2.210,2.256,2.324};
	const double Matthews::_wcoeff[15] = {0.213,0.214,0.220,0.225,0.226,0.231,0.236,0.241,
		0.242,0.244,0.244,0.244,0.245,0.446,0.327};
	const double Matthews::_acoeff[15] = {28.38,102.7,187.5,339.3,434.1,540.5,686.2,767.9,
		835.9,856.9,854.0,849.6,846.7,136.6,47.10};
	const double Matthews::_scoeff[15] = {0.953,0.807,0.775,0.702,0.648,0.640,0.635,0.589,
		0.584,0.542,0.500,0.485,0.480,1.180,0.466};
	
	const float Matthews::_densProt = (1.0f/0.74f);
	const float Matthews::_densDNA = (1.0f/0.5f);
	
	Matthews::Matthews(bool protein, bool rna) : _protein(protein), _rna(rna)
	{
		if (_protein ) {
			if (_rna ) {
				_density = 0.25*_densDNA + 0.75*_densProt; //assume 25%/75% DNA/protein, as in Kantardjieff
			} else {
				_density = _densProt;
			}
		} else {
			_density = _densDNA;
		}
	}
	
	int Matthews::operator() (clipper::Cell& cell, clipper::Spacegroup& spacegroup, clipper::SEQfile& file, float reso)
	{
		clipper::MMoleculeSequence seq;
		file.import_molecule_sequence( seq );
		clipper::MPolymerSequence poly = seq[0];
		clipper::String sequence = poly.sequence();
		
		std::vector<int> numatoms(5,0);
		for (int i=0; i<sequence.length(); i++) {
			for (int j=0; j<21; j++) {
				if (sequence[i] == _ResidueName1[j]) {
					numatoms[0] += _Catoms[j];
					numatoms[1] += _Natoms[j];
					numatoms[2] += _Oatoms[j];
					numatoms[3] += _Hatoms[j];
					numatoms[4] += _Satoms[j];
					break;
				}
			}
		}
		_weight = 12.0*numatoms[0] + 14.0*numatoms[1] + 16.0*numatoms[2] +1.0*numatoms[3] +32.0*numatoms[4];
		
		int nmol = 0;
		
		if (_protein ) {
			if (_rna ) {
				nmol = calcComp(cell, spacegroup, reso);
			} else {
				nmol = calcProt(cell, spacegroup, reso);
			}
		} else {
			nmol = calcDNA(cell, spacegroup, reso);
		}
		return nmol;
		
	}
	
	int Matthews::operator() (clipper::Cell& cell, clipper::Spacegroup& spacegroup, int nresidues, float reso)
	{
		_reso = reso;
		// use average weight of C G A T residue (NOT base pair)
		float dna_weight = 325.96*nresidues;
		
		float prot_weight = (12.0*5 + 14.0*1.35 + 16.0*1.5 +1.0*8 +32.0*0.05)*nresidues;
		
		int nmol = 0;
		
		if (_protein ) {
			if (_rna ) {
				_weight = 0.25*dna_weight + 0.75*prot_weight; //assume 25%/75% DNA/protein, as in Kantardjieff
				nmol = calcComp(cell, spacegroup, reso);
			} else {
				_weight = prot_weight;
				nmol = calcProt(cell, spacegroup, reso);
			}
		} else {
			_weight = dna_weight;
			nmol = calcDNA(cell, spacegroup, reso);
		}
		return nmol;
	}
	
	int Matthews::calcProt(clipper::Cell& cell, clipper::Spacegroup& spacegroup, float reso)
	{
		float volume = cell.volume();
		float nsym = spacegroup.num_symops();
		int maxmols1 = std::floor(0.602f*volume*(_density)/(_weight*nsym) )+1;
		int nmol = 0;
		_reso = reso;
		
		float maxprob = 0.0f;
		
		_cmath.resize(maxmols1);
		_solvent.resize(maxmols1);
		_prob.resize(maxmols1);
		_probt.resize(maxmols1);
		
		for (int i = 0 ; i != maxmols1 ; ++i ) {
			_prob[i] = 0.0f;
			_probt[i] = 0.0f;
		}
		
		int imols = 1;
		
		for ( ; imols != maxmols1 ; ++imols ) {
			_cmath[imols] = volume/(_weight*float(imols)*nsym);
			float solvent = (1.0-1.0f/(0.602f*_cmath[imols]*_density))*100.0;
			if (_solvent[imols] < 0.0f ) break;
			_solvent[imols] = solvent;
			
			for (int j = 0; j != 13 ; ++j ) {
				_probt[imols] += vmProb(_cmath[imols],_p0[j],_vmbar[j],_wcoeff[j],_acoeff[j],_scoeff[j]);
			}
			_probt[0] += _probt[imols];
			if ( _probt[imols] >= maxprob) {
				nmol = imols;
				maxprob = _probt[imols];
			}
		}
		
		if ( reso < 99.0f ) {
			maxprob = 0.0f;
			int ireso = 0;
			for ( ; ireso != 13 ; ++ireso ) {
				if ( reso <= _rbin[ireso] ) break;
			}
			for (int i = 1 ; i != imols ; ++i ) {
				_prob[i] = vmProb(_cmath[i],_p0[ireso],_vmbar[ireso],_wcoeff[ireso],_acoeff[ireso],_scoeff[ireso]);
				_prob[0] += _prob[i];
				if ( _prob[i] >= maxprob) {
					nmol = i;
					maxprob = _prob[i];
				}				
			}
			for (int i = 1 ; i != imols ; ++i ) _prob[i]/=_prob[0];
		}
		for (int i = 1 ; i != imols ; ++i ) _probt[i]/=_probt[0];	
		
		return nmol;
	}
	
	
	int Matthews::calcComp(clipper::Cell& cell, clipper::Spacegroup& spacegroup, float reso)
	{
		float volume = cell.volume();
		float nsym = spacegroup.num_symops();
		int maxmols1 = std::floor(0.602f*volume*(_density)/(0.602f*_weight*nsym) )+1;
		int nmol = 0;
		int COMP = 13;
		
		float maxprob = 0.0f;
		
		_cmath.resize(maxmols1);
		_solvent.resize(maxmols1);
		_prob.resize(maxmols1);
		_probt.resize(maxmols1);
		
		for (int i = 1 ; i != maxmols1 ; ++i ) {
			_cmath[i] = volume/(_weight*i*nsym);
			_solvent[i] = (1.0-1.0f/(0.602f*_cmath[i]*_density))*100.0;
			
			_prob[i] = vmProb(_cmath[i],_p0[COMP],_vmbar[COMP],_wcoeff[COMP],_acoeff[COMP],_scoeff[COMP]);
			_prob[0] += _prob[i];
			if ( _prob[i] >= maxprob) {
				nmol = i;
				maxprob = _probt[i];
			}				
		}
		for (int i = 1 ; i != maxmols1 ; ++i ) _prob[i]/=_prob[0];
		
		return nmol;
	}
	
	
	
	int Matthews::calcDNA(clipper::Cell& cell, clipper::Spacegroup& spacegroup, float reso)
	{
		float volume = cell.volume();
		float nsym = spacegroup.num_symops();
		int maxmols1 = std::floor(0.602f*volume*(_density)/(0.602f*_weight*nsym) )+1;
		int nmol = 0;
		int DNA = 14;
		
		float maxprob = 0.0f;
		
		_cmath.resize(maxmols1);
		_solvent.resize(maxmols1);
		_prob.resize(maxmols1);
		_probt.resize(maxmols1);
		
		for (int i = 1 ; i != maxmols1 ; ++i ) {
			_cmath[i] = volume/(_weight*i*nsym);
			_solvent[i] = (1.0-1.0f/(0.602f*_cmath[i]*_density))*100.0;
			_prob[i] = vmProb(_cmath[i],_p0[DNA],_vmbar[DNA],_wcoeff[DNA],_acoeff[DNA],_scoeff[DNA]);
			_prob[0] += _prob[i];
			if ( _prob[i] >= maxprob) {
				nmol = i;
				maxprob = _probt[i];
			}				
		}
		for (int i = 1 ; i != maxmols1 ; ++i ) _prob[i]/=_prob[0];
		
		return nmol;
	}
	
	
	
	float Matthews::matthews(clipper::Cell& cell, clipper::Spacegroup& spacegroup, int nmol)
	{
		/*  Assume 64% of cell volume solvent in dna crystal, with density 1.0/0.5
		 Assume 60% of cell volume solvent in dna/protein crystal..
		 Assume 47% of cell volume solvent in protein crytal, with density 1.0/0.74
		 */
		float volume = cell.volume();
		float nsym = spacegroup.num_symops();
		
		float dna_weight = (1.0-0.64)*_densDNA*(0.602*volume)/(nmol*nsym);
		float prot_weight = (1.0-0.47)*_densProt*(0.602*volume)/(nmol*nsym);
		if (_protein ) {
			if (_rna ) {
				_weight = 0.25*dna_weight + 0.75*prot_weight; //assume 25%/75% DNA/protein, as in Kantardjieff
			} else {
				_weight = prot_weight;
			}
		} else {
			_weight = dna_weight;
		}
		return volume/(_weight*nmol*nsym);
	}
	
	void Matthews::summary(){
		int maxmol1 = _cmath.size();
		printf("$TABLE: Matthews coefficients:\n");
		printf("$GRAPHS");
		if (_protein ) {
			if (_rna ) {
				printf(": 25%% DNA/75%% protein, as in Kantardjieff");
			} else {
				printf(": Protein crystal");
			}
		} else {
			printf(": DNA crystal");
		}
		if ( _reso < 99.0f ) {
			printf(" computed at resolution of %5.3f :A:1,2,3,4\n", _reso,maxmol1);
			printf("$$ Nmol/asym Matthews_Coeff sovlent_frac P(%5.3f) $$\n$$\n",_reso);
			
			for (int i=1; i != maxmol1 ; ++i) {
				printf("  %3d       %6.2f          %6.2f       %6.2f\n", i, _cmath[i], _solvent[i]/100.0f, _prob[i]);
			}
		} else {
			printf(" :A:1,2,3,4 \n", maxmol1);
			printf("$$ Nmol/asym Matthews Coeff solvent_frac P() $$\n$$\n");
			
			for (int i=1; i != maxmol1 ; ++i) {
				printf("  %3d       %6.2f          %6.2f       %6.2f\n", i, _cmath[i], _solvent[i]/100.0f, _probt[i]);
			}
		}
		printf("$$\n\n");
		
	}
	
}
