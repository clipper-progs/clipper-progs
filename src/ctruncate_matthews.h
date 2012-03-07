//
//     CTRUNCATE
//     Copyright (C) 2010-2011 Charles Ballard
//
//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Application.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//

#ifndef __CTRUNCATE_MATTHEWS_H
#define __CTRUNCATE_MATTHEWS_H

#include "clipper/clipper.h"
#include "clipper/clipper-contrib.h"
#include "clipper/clipper-ccp4.h"
#include "clipper/clipper-minimol.h"

namespace ctruncate {
	
	class Matthews
	{
	public:
		Matthews(bool protein = true, bool rna = false);
		
		int operator() (clipper::Cell& cell, clipper::Spacegroup& spacegroup, clipper::SEQfile& file, float reso = 99.0f );
		int operator() (clipper::Cell& cell, clipper::Spacegroup& spacegroup, int nresidues, float reso = 99.0f );
		
		float matthews(clipper::Cell& cell, clipper::Spacegroup& spacegroup, int nmol);
		float matthews(int i) { return _cmath[i]; }
		
		float solvent(int i) { return _solvent[i]; }
		
		float prob(int i) { return _prob[i]; }
		
		void summary();
		
	private:
		// prob calc
		float vmProb(double x, double y0, double xc, double wt, double a, double s) 
		{
			double z = (x-xc)/wt;
			return y0+a*(std::exp(-std::exp(-z)-z*s+1.0));
		}
		int calcProt(clipper::Cell& cell, clipper::Spacegroup& spacegroup, float reso);
		int calcComp(clipper::Cell& cell, clipper::Spacegroup& spacegroup, float reso);
		int calcDNA(clipper::Cell& cell, clipper::Spacegroup& spacegroup, float reso);
		
		/* array of contents based on mmdb single letter codes */
		static const char _ResidueName1[21];
		static const int _Catoms[21];
		static const int _Hatoms[21];
		static const int _Natoms[21];
		static const int _Oatoms[21];
		static const int _Satoms[21];
		
		/* Arrays   :       1: rbin, 2: p0, 3: vmbar, 4: w, 5: a, 6: s
		 columns :        1-13 protein data inclusive to corresponding binmax value
		 14 DNA (all resolutions)
		 15 DNA/Protein complexes(25%/75%) */
		static const double _rbin[15];
		static const double _p0[15];
		static const double _vmbar[15];
		static const double _wcoeff[15];
		static const double _acoeff[15];
		static const double _scoeff[15];
		
		static const float _densProt;
		static const float _densDNA;
		
		bool _protein;
		bool _rna;
		float _weight;
		float _density;
		float _reso;
		std::vector<float> _cmath;
		std::vector<float> _solvent;
		std::vector<float> _prob;
		std::vector<float> _probt;
	};
	
}

#endif