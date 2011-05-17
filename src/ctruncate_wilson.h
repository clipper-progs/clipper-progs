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
	
	std::vector<float> wilson_plot(clipper::HKL_data<clipper::data32::I_sigI>& isig, float maxres, int nbins, 
								   CCP4Program& prog, clipper::HKL_data<clipper::data32::I_sigI>& norm);
	
	std::vector<float> wilson_plot(clipper::HKL_data<clipper::data32::I_sigI>& isig, int nresidues, float maxres, 
								   int nbins, CCP4Program& prog, clipper::HKL_data<clipper::data32::I_sigI>& norm);
	
	std::vector<float> wilson_plot(clipper::HKL_data<clipper::data32::I_sigI>& isig, clipper::MPolymerSequence& poly, 
								   float maxres, int nbins, CCP4Program& prog, clipper::HKL_data<clipper::data32::I_sigI>& norm);
	
	
	class WilsonB {
	public:
		enum MODE {STRAIGHT, BEST };
		
		WilsonB ( WilsonB::MODE _mode = BEST ) : mode(_mode), nresidue_supplied(false), sequence_supplied(false), _a(-1.0), 
		_bscale(1.0f), _boff(0.0f)
		{ numatoms.resize(5); }
		
		float operator()(clipper::HKL_data<clipper::data32::I_sigI>& isig, ctruncate::Rings* ice = NULL);
		
		float operator()(clipper::HKL_data<clipper::data32::I_sigI>& isig, int nresidues, ctruncate::Rings* ice = NULL);
		
		float operator()(clipper::HKL_data<clipper::data32::I_sigI>& isig, clipper::MPolymerSequence& poly, ctruncate::Rings* ice = NULL);
		
		float a() {return _a;}
		float b() {return _b;}
		float B() {return 2.0*_a; }
		float intercept () { return exp(_b); }
		float siga() {return _siga;}
		float sibb() {return _sigb;}
		
		void summary();
		void plot ();
		void plot(clipper::HKL_data<clipper::data32::I_sigI>& ref, clipper::String& name);
		
	private:
		MODE mode;
		float _a;  // slope
		float _b;  // intercept
		float _siga; // uncertainty in intercept
		float _sigb; //uncertainty in b-value
		float _bscale; // scale best to scattering
		float _boff; //offset on best to scattering
		clipper::HKL_data<clipper::data32::I_sigI>* intensity;
		
		bool nresidue_supplied, sequence_supplied; //book keeping
		int nresidues; // number of residues
		std::vector<int> numatoms;  // number of atoms
		float maxres;
		
		static const std::string AtomNames[5];  //atom names
		static const char ResidueNames[21];  //residue names
		static const int Catoms[21];         // number of C per residue
		static const int Hatoms[21];         // number of H per residue
		static const int Natoms[21];         // number of N per residue
		static const int Oatoms[21];         // number of O per residue
		static const int Satoms[21];         // number of S per residue
		
		void wilson_straight(clipper::HKL_data<clipper::data32::I_sigI>& isig, ctruncate::Rings& ice);
		void wilson_best(clipper::HKL_data<clipper::data32::I_sigI>& isig, ctruncate::Rings& ice);
	};
		
}
#endif
