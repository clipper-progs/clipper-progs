//
//     CTRUNCATE
//     Copyright (C) 2006-2011 Norman Stein, Charles Ballard
//
//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Application.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//

#ifndef __CTRUNCATE_WILSON_H
#define __CTRUNCATE_WILSON_H

#include "clipper/clipper.h"
#include "clipper/clipper-minimol.h"
#include "clipper/clipper-ccp4.h"
#include "ctruncate_utils.h"

namespace ctruncate {
	
	class Scattering {
	public:
		enum MODE { PROTEIN, NUCLEIC, COMPLEX };
        enum AtomType_t { H = 0, C = 5, N = 6, O = 7, P=14, S=15 };

		Scattering( MODE mode=PROTEIN) : _mode(mode) { init(); }
		~Scattering() {}
        
        MODE mode() const { return _mode; }
        MODE mode( MODE mode) { MODE tmp = _mode; _mode = mode; return tmp; }
        
        clipper::ftype numAtoms(const std::string&) const;
        clipper::ftype numAtoms(int) const;
        
        clipper::ftype f(const clipper::ftype invresolsq, clipper::ftype power=2.0) const;
        
        clipper::ftype operator()(const clipper::HKL_data_base& data, const int nres, clipper::ftype power=2.0);
        clipper::ftype operator()(const clipper::HKL_data_base& data, const clipper::MPolymerSequence& poly, clipper::ftype power=2.0);
        clipper::ftype operator()(const clipper::HKL_data_base& data, clipper::ftype solvent, clipper::ftype power=2.0);
        clipper::ftype operator()(const clipper::HKL_data_base& data, clipper::MMDBManager& mmdb, clipper::ftype power=2.0);
        clipper::ftype operator()(const clipper::HKL_data_base& data, std::vector< clipper::Atom>& atoms, clipper::ftype power=2.0);
		clipper::ftype operator()(const int nres, const clipper::Spacegroup& spg, clipper::ftype power=2.0 );
		clipper::ftype operator()(const clipper::MPolymerSequence& poly, const clipper::Spacegroup& spg, clipper::ftype power=2.0);
		clipper::ftype operator()(const clipper::Cell& cell, const clipper::Spacegroup& spg, clipper::ftype solvent=0.5, clipper::ftype power=2.0 );
		clipper::ftype operator()(clipper::MMDBManager& mmdb, const clipper::Spacegroup& spg, clipper::ftype power=2.0);
		clipper::ftype operator()(std::vector<clipper::Atom>& mmdb, const clipper::Spacegroup& spg, clipper::ftype power=2.0 );
		
	private:
		void init();
		static const std::string ProteinAtomNames[5];  //atom names
		static const char ProteinResidueNames[20];  //residue names
		static const int ProteinCatoms[20];         // number of C per residue
		static const int ProteinHatoms[20];         // number of H per residue
		static const int ProteinNatoms[20];         // number of N per residue
		static const int ProteinOatoms[20];         // number of O per residue
		static const int ProteinSatoms[20];         // number of S per residue
		static const clipper::ftype ProteinComp[5];          // average composition
		
		static const std::string NucleicAtomNames[5];
		static const char NucleicResidueNames[5];
		static const int NucleicCatoms[5];
		static const int NucleicHatoms[5];
		static const int NucleicNatoms[5];
		static const int NucleicOatoms[5];
		static const int NucleicPatoms[5];
		static const clipper::ftype NucleicComp[5]; 
		
        MODE _mode;                                   // mode of operation
		std::vector<std::string> atomtype;            // name of scatterer
		std::vector<clipper::ftype> scattering;       // number of electrons per scatter
        std::vector<clipper::ftype> numatoms;         // number of atomtype in unit cell
		
	};
	
	class WilsonB {
	public:
		enum MODE {STRAIGHT, BEST, RNA };
		
		WilsonB ( WilsonB::MODE _mode = BEST ) : mode(_mode), nresidue_supplied(false), sequence_supplied(false), _a(-1.0), _totalscat(300,0.0)
        { scattering.mode((mode == RNA ) ? Scattering::NUCLEIC : Scattering::PROTEIN ); }
		
		clipper::ftype operator()(const clipper::HKL_data_base& data, const clipper::Range<clipper::ftype>* range = NULL, ctruncate::Rings* ice = NULL);
		
		clipper::ftype operator()(const clipper::HKL_data_base& data, int nresidues, const clipper::Range<clipper::ftype>* range = NULL, ctruncate::Rings* ice = NULL);
		
		clipper::ftype operator()(const clipper::HKL_data_base& data, clipper::MPolymerSequence& poly, const clipper::Range<clipper::ftype>* range = NULL, ctruncate::Rings* ice = NULL);
		
		float a() {return _a;}
		float b() {return _b;}
		float B() {return 2.0*_a; }
		float intercept () { return exp(_b); }
		float siga() {return _siga;}
		float sibb() {return _sigb;}
		
		void summary();
		void plot (int nbins = 60);
		void output();
		std::stringstream& xml_output(std::stringstream&);
        
        clipper::ftype f(clipper::ftype);
        clipper::ftype f(clipper::HKL_data_base::HKL_reference_index&);
		
	private:
        ctruncate::Scattering scattering;   //hold information on cell contents
		MODE mode;          //mode
        clipper::ftype _a;  // slope
        clipper::ftype _b;  // intercept
        clipper::ftype _siga; // uncertainty in intercept
        clipper::ftype _sigb; //uncertainty in b-value
		const clipper::HKL_data_base* intensity;
		clipper::Range<clipper::ftype> activeRange; //range of data used in calculation
        
		bool nresidue_supplied, sequence_supplied; //book keeping
		int nresidues; // number of residues
		//std::vector<int> numatoms;  // number of atoms
        std::vector<clipper::ftype> _totalscat;  //number of scattering electrons in cell (scale for BEST)
		float maxres;
		
        static const std::string AtomNames[];
		
		clipper::ftype obs(const clipper::HKL_data_base&, const clipper::HKL_data_base::HKL_reference_index&);              // return observed value as I or f*f
		clipper::ftype sigobs(const clipper::HKL_data_base&, const clipper::HKL_data_base::HKL_reference_index&);           // return suitable sigma
		void wilson_straight(const clipper::HKL_data_base&, clipper::Range<clipper::ftype>& range, ctruncate::Rings& ice);
		void wilson_best(const clipper::HKL_data_base&, clipper::Range<clipper::ftype>& range, ctruncate::Rings& ice);
		void wilson_rna(const clipper::HKL_data_base&, clipper::Range<clipper::ftype>& range, ctruncate::Rings& ice);
        
	};
		
}
#endif
