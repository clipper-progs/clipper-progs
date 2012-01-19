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
#include "alt_hkl_datatypes.h"

namespace ctruncate {
	
	int cumulative_plot(clipper::HKL_data<clipper::data32::I_sigI>& isig, clipper::ResolutionFn& Sigma);
	
	int cumulative_plot(clipper::HKL_data<clipper::data32::I_sigI>& isig, clipper::HKL_data<clipper::data32::I_sigI>& Sigma);
	
	void yorgo_modis_plot(clipper::HKL_data<clipper::data32::F_sigF>& fsig, float maxres, int nbins, CCP4Program& prog, clipper::U_aniso_orth& uao);
	
	void yorgo_modis_plot(clipper::HKL_data<clipper::data32::I_sigI>& isig, float maxres, int nbins, CCP4Program& prog, clipper::U_aniso_orth& uao);
	
	//--------------------------------------------------------------
	
	class PattPeak 
	{
	public:
		PattPeak( float maxinvres, int nbins = 20, float temp = -1.0f );
		
		float operator() (clipper::BasisFn_spline& basis_fo, clipper::ResolutionFn& resol_fn);
		float operator() (clipper::HKL_data<clipper::data32::I_sigI>& Sigma);
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
	
// vs resolution and half dataset CC
template<class T> class AnomStats
{
public:
    enum TYPE { F, I };
    AnomStats(clipper::HKL_data<clipper::datatypes::J_sigJ_ano<T> >& hkl_data, int nbins=60);
    AnomStats(clipper::HKL_data<clipper::datatypes::G_sigG_ano<T> >& hkl_data, int nbins=60);
protected:
    const T&    obs_pl( const clipper::datatypes::J_sigJ_ano<T>& f ) { return f.I_pl(); }
    const T&    obs_mi( const clipper::datatypes::J_sigJ_ano<T>& f ) { return f.I_mi(); }
    const T& sigobs_pl( const clipper::datatypes::J_sigJ_ano<T>& f ) { return f.sigI_pl(); }
    const T& sigobs_mi( const clipper::datatypes::J_sigJ_ano<T>& f ) { return f.sigI_pl(); }
    const T&    obs_pl( const clipper::datatypes::G_sigG_ano<T>& f ) { return f.f_pl(); }
    const T&    obs_mi( const clipper::datatypes::G_sigG_ano<T>& f ) { return f.f_mi(); }
    const T& sigobs_pl( const clipper::datatypes::G_sigG_ano<T>& f ) { return f.sigf_pl(); }
    const T& sigobs_mi( const clipper::datatypes::G_sigG_ano<T>& f ) { return f.sigf_pl(); }
 
private:
    int _nbins;
};
    
}
#endif
