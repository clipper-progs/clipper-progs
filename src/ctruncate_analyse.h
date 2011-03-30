//
//     CTRUNCATE
//     Copyright (C) 2006-2011 Norman Stein, Charles Ballard
//
//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Application.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//

#ifndef __CTRUNCATE_ANALYSE_H
#define __CTRUNCATE_ANALYSE_H

#include "clipper/clipper.h"
#include "clipper/clipper-ccp4.h"

namespace ctruncate {
	
	int cumulative_plot(clipper::HKL_data<clipper::data32::I_sigI>& isig, clipper::ResolutionFn& Sigma);
	
	void yorgo_modis_plot(clipper::HKL_data<clipper::data32::F_sigF>& fsig, float maxres, int nbins, CCP4Program& prog);
	
	//--------------------------------------------------------------
	
	class PattPeak 
	{
	public:
		PattPeak( float maxinvres, int nbins = 20, float temp = -1.0f );
		
		float operator() (clipper::BasisFn_spline& basis_fo, clipper::ResolutionFn& resol_fn);
		void setBins(int nbins) { _nbins = nbins; _patterson.resize(nbins); }
		void setMaxInvRes(float maxinvres) { _maxinvres = maxinvres; }
		void setWidth(float width) { _width = width; }
		
		float p000() { return _P000; }
		float sigma() { return _sigma; }
		float optRes();
		
	private:
		void calcOriginPeak(clipper::BasisFn_spline& basis_fo, clipper::ResolutionFn& resol_fn);
		void fitOriginPeak(clipper::BasisFn_spline& basis_fo, clipper::ResolutionFn& resol_fn);
		
		float _P000;
		float _sigma;
		
		int _nbins; // number of bins
		float _maxinvres;  // maximum 1/res to use
		float _width;
		std::vector<float> _patterson;
	};
	
}

#endif
