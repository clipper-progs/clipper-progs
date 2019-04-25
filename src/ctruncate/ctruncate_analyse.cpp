//
//     CTRUNCATE
//     Copyright (C) 2006-2011 Norman Stein, Charles Ballard
//
//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Application.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//

#include "ctruncate_analyse.h"
#include "ctruncate_utils.h"
#include "intensity_target.h"
#include "ctruncate_wilson.h"
#include "cpsf_utils.h"
#include "best.h"
#include "best_rna.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

#define ASSERT assert
#include <assert.h>

namespace ctruncate {
	
	int cumulative_plot(clipper::HKL_data<clipper::data32::I_sigI>& isig, clipper::HKL_data<clipper::data32::I_sigI>& Sigma)
	{
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		
		// construct cumulative distribution function for intensity (using Z rather than E)
		clipper::Range<double> intensity_range_centric;
		clipper::Range<double> intensity_range_acentric;
		// changed from jsig to isig
		
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			if ( !isig[ih].missing() ) {
				if ( ih.hkl_class().centric() ) intensity_range_centric.include( (isig[ih].I()/ih.hkl_class().epsilon() )/Sigma[ih].I() );
				else intensity_range_acentric.include( (isig[ih].I()/ih.hkl_class().epsilon() ) /Sigma[ih].I() );
			}
		}
		//printf("C2: %20.5f %20.5f\n",intensity_range_centric.max(),intensity_range_centric.min());
		//printf("A2: %20.5f %20.5f\n",intensity_range_acentric.max(),intensity_range_acentric.min());
		
		int ncentric = 0;
		clipper::Generic_ordinal intensity_ord_c;
		clipper::Generic_ordinal intensity_ord_a;
		intensity_ord_c.init( intensity_range_centric );
		intensity_ord_a.init( intensity_range_acentric );
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			if ( !isig[ih].missing() ) {
				if ( ih.hkl_class().centric() ) {
					intensity_ord_c.accumulate( ( isig[ih].I()/ih.hkl_class().epsilon() )/Sigma[ih].I()  );
					ncentric++;
				}
				else intensity_ord_a.accumulate( ( isig[ih].I()/ih.hkl_class().epsilon() )/Sigma[ih].I() );
			}
		}
		intensity_ord_c.prep_ordinal();
		intensity_ord_a.prep_ordinal();
		
		// theoretical values for cumulative intensity distribution
		double acen[51] = {0.0,
			0.0392106, 0.0768837, 0.1130796, 0.1478562, 0.1812692, 0.2133721, 0.2442163, 0.2738510, 0.3023237, 0.3296800,
			0.3559636, 0.3812166, 0.4054795, 0.4287909, 0.4511884, 0.4727076, 0.4933830, 0.5132477, 0.5323336, 0.5506710,
			0.5682895, 0.5852171, 0.6014810, 0.6171071, 0.6321206, 0.6465453, 0.6604045, 0.6737202, 0.6865138, 0.6988058,
			0.7106158, 0.7219627, 0.7328647, 0.7433392, 0.7534030, 0.7630722, 0.7723623, 0.7812881, 0.7898639, 0.7981035,
			0.8060200, 0.8136260, 0.8209339, 0.8279551, 0.8347011, 0.8411826, 0.8474099, 0.8533930, 0.8591416, 0.8646647};
		
		double cen[51] = {0.0,
			0.1585194, 0.2227026, 0.2709655, 0.3108435, 0.3452792, 0.3757939, 0.4032988, 0.4283924, 0.4514938, 0.4729107,
			0.4928775, 0.5115777, 0.5291583, 0.5457398, 0.5614220, 0.5762892, 0.5904133, 0.6038561, 0.6166715, 0.6289066,
			0.6406032, 0.6517983, 0.6625250, 0.6728131, 0.6826895, 0.6921785, 0.7013024, 0.7100815, 0.7185345, 0.7266783,
			0.7345289, 0.7421010, 0.7494079, 0.7564625, 0.7632764, 0.7698607, 0.7762255, 0.7823805, 0.7883348, 0.7940968,
			0.7996745, 0.8050755, 0.8103070, 0.8153755, 0.8202875, 0.8250491, 0.8296659, 0.8341433, 0.8384867, 0.8427008};
		
		double acen_twin[51] = {0.0,
			0.0030343, 0.0115132, 0.0245815, 0.0414833, 0.0615519, 0.0842006, 0.1089139, 0.1352404, 0.1627861, 0.1912079, 
			0.2202081, 0.2495299, 0.2789524, 0.3082868, 0.3373727, 0.3660750, 0.3942806, 0.4218963, 0.4488460, 0.4750691, 
			0.5005177, 0.5251562, 0.5489585, 0.5719077, 0.5939942, 0.6152149, 0.6355726, 0.6550744, 0.6737317, 0.6915590, 
			0.7085736, 0.7247951, 0.7402450, 0.7549459, 0.7689218, 0.7821971, 0.7947971, 0.8067470, 0.8180725, 0.8287987, 
			0.8389511, 0.8485543, 0.8576328, 0.8662106, 0.8743109, 0.8819565, 0.8891694, 0.8959710, 0.9023818, 0.9084218}; 
		
		
		printf("$TABLE: Cumulative intensity distribution:\n");
		printf("$GRAPHS");
		printf(": Cumulative intensity distribution (Acentric and centric):N:1,2,3,4,5,6:\n$$");
		printf(" Z Acent_theor Acent_twin Acent_obser Cent_theor Cent_obser $$\n$$\n");
		double x = 0.0;
		double deltax=0.04;
		for (int i=0; i<=50; i++) {
			if (ncentric) printf("%10.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n", x, acen[i], acen_twin[i], intensity_ord_a.ordinal(x), cen[i], intensity_ord_c.ordinal(x));
			else printf("%10.5f %8.5f %8.5f %8.5f %8.5f -\n", x, acen[i], acen_twin[i], intensity_ord_a.ordinal(x), cen[i]);
			x += deltax;
		}
		printf("$$\n\n");
		
		// count entries where acentric distribution lower than expected - sign of twinning
		x = 0.08;
		deltax = 0.12;
		int ntw = 0;
		for (int i=0;i<4; i++) {
			if ( (acen[3*i+2] - intensity_ord_a.ordinal(x))/acen[3*i+2] > 0.4 ) ntw ++;
			x += deltax;
		}
		return ntw;
	}
	
	int cumulative_plot(clipper::HKL_data<clipper::data32::I_sigI>& isig, clipper::ResolutionFn& Sigma)
	{
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		
		// construct cumulative distribution function for intensity (using Z rather than E)
		clipper::Range<double> intensity_range_centric;
		clipper::Range<double> intensity_range_acentric;
		// changed from jsig to isig
		
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			if ( !isig[ih].missing() ) {
				if ( ih.hkl_class().centric() ) intensity_range_centric.include( (isig[ih].I()/ih.hkl_class().epsilon() )/Sigma.f(ih) );
				else intensity_range_acentric.include( (isig[ih].I()/ih.hkl_class().epsilon() ) /Sigma.f(ih) );
			}
		}
		//printf("C2: %20.5f %20.5f\n",intensity_range_centric.max(),intensity_range_centric.min());
		//printf("A2: %20.5f %20.5f\n",intensity_range_acentric.max(),intensity_range_acentric.min());
		
		int ncentric = 0;
		clipper::Generic_ordinal intensity_ord_c;
		clipper::Generic_ordinal intensity_ord_a;
		intensity_ord_c.init( intensity_range_centric );
		intensity_ord_a.init( intensity_range_acentric );
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			if ( !isig[ih].missing() ) {
				if ( ih.hkl_class().centric() ) {
					intensity_ord_c.accumulate( ( isig[ih].I()/ih.hkl_class().epsilon() )/Sigma.f(ih)  );
					ncentric++;
				}
				else intensity_ord_a.accumulate( ( isig[ih].I()/ih.hkl_class().epsilon() )/Sigma.f(ih) );
			}
		}
		intensity_ord_c.prep_ordinal();
		intensity_ord_a.prep_ordinal();
		
		// theoretical values for cumulative intensity distribution
		double acen[51] = {0.0,
			0.0392106, 0.0768837, 0.1130796, 0.1478562, 0.1812692, 0.2133721, 0.2442163, 0.2738510, 0.3023237, 0.3296800,
			0.3559636, 0.3812166, 0.4054795, 0.4287909, 0.4511884, 0.4727076, 0.4933830, 0.5132477, 0.5323336, 0.5506710,
			0.5682895, 0.5852171, 0.6014810, 0.6171071, 0.6321206, 0.6465453, 0.6604045, 0.6737202, 0.6865138, 0.6988058,
			0.7106158, 0.7219627, 0.7328647, 0.7433392, 0.7534030, 0.7630722, 0.7723623, 0.7812881, 0.7898639, 0.7981035,
			0.8060200, 0.8136260, 0.8209339, 0.8279551, 0.8347011, 0.8411826, 0.8474099, 0.8533930, 0.8591416, 0.8646647};
		
		double cen[51] = {0.0,
			0.1585194, 0.2227026, 0.2709655, 0.3108435, 0.3452792, 0.3757939, 0.4032988, 0.4283924, 0.4514938, 0.4729107,
			0.4928775, 0.5115777, 0.5291583, 0.5457398, 0.5614220, 0.5762892, 0.5904133, 0.6038561, 0.6166715, 0.6289066,
			0.6406032, 0.6517983, 0.6625250, 0.6728131, 0.6826895, 0.6921785, 0.7013024, 0.7100815, 0.7185345, 0.7266783,
			0.7345289, 0.7421010, 0.7494079, 0.7564625, 0.7632764, 0.7698607, 0.7762255, 0.7823805, 0.7883348, 0.7940968,
			0.7996745, 0.8050755, 0.8103070, 0.8153755, 0.8202875, 0.8250491, 0.8296659, 0.8341433, 0.8384867, 0.8427008};
		
		double acen_twin[51] = {0.0,
			0.0030343, 0.0115132, 0.0245815, 0.0414833, 0.0615519, 0.0842006, 0.1089139, 0.1352404, 0.1627861, 0.1912079, 
			0.2202081, 0.2495299, 0.2789524, 0.3082868, 0.3373727, 0.3660750, 0.3942806, 0.4218963, 0.4488460, 0.4750691, 
			0.5005177, 0.5251562, 0.5489585, 0.5719077, 0.5939942, 0.6152149, 0.6355726, 0.6550744, 0.6737317, 0.6915590, 
			0.7085736, 0.7247951, 0.7402450, 0.7549459, 0.7689218, 0.7821971, 0.7947971, 0.8067470, 0.8180725, 0.8287987, 
			0.8389511, 0.8485543, 0.8576328, 0.8662106, 0.8743109, 0.8819565, 0.8891694, 0.8959710, 0.9023818, 0.9084218}; 
		
		
		printf("$TABLE: Cumulative intensity distribution:\n");
		printf("$GRAPHS");
		printf(": Cumulative intensity distribution (Acentric and centric):N:1,2,3,4,5,6:\n$$");
		printf(" Z Acent_theor Acent_twin Acent_obser Cent_theor Cent_obser $$\n$$\n");
		double x = 0.0;
		double deltax=0.04;
		for (int i=0; i<=50; i++) {
			if (ncentric) printf("%10.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n", x, acen[i], acen_twin[i], intensity_ord_a.ordinal(x), cen[i], intensity_ord_c.ordinal(x));
			else printf("%10.5f %8.5f %8.5f %8.5f %8.5f -\n", x, acen[i], acen_twin[i], intensity_ord_a.ordinal(x), cen[i]);
			x += deltax;
		}
		printf("$$\n\n");
		
		// count entries where acentric distribution lower than expected - sign of twinning
		x = 0.08;
		deltax = 0.12;
		int ntw = 0;
		for (int i=0;i<4; i++) {
			if ( (acen[3*i+2] - intensity_ord_a.ordinal(x))/acen[3*i+2] > 0.4 ) ntw ++;
			x += deltax;
		}
		return ntw;
	}
	
	
	//--------------------------------------------------------------
	
	PattPeak::PattPeak(float maxinvres, int nbins, float temp ) : _maxinvres(maxinvres), _nbins(nbins), _patterson(nbins,0.0f)
	{
		float coef = 1.5f;
		float dmax = (sqrt(1.0f/maxinvres)/3.0f)*2.0f*1.5f;
		float btemp = ( temp > 0.0f ) ? temp : 0.0f ;
		float dmax1 =  std::sqrt(btemp/(clipper::Util::twopi()*clipper::Util::twopi()*2.0f) )*2.0f*1.5f ;
		_width = ( dmax1 > dmax ) ? dmax1 : dmax ;
	}
	
	float PattPeak::operator()(clipper::BasisFn_spline& basis_fo, clipper::ResolutionFn& resol_fn)
	{
		calcOriginPeak(basis_fo, resol_fn);
		fitOriginPeak(basis_fo, resol_fn);
		
		return optRes();
	}
	
	float PattPeak::operator()(clipper::HKL_data<clipper::data32::I_sigI>& Sigma)
	{
		int nprm = 60;
		const clipper::HKL_info& hklinf = Sigma.hkl_info();
		std::vector<double> params( nprm, 1.0 );
		clipper::BasisFn_spline basis_fo( Sigma, nprm, 2.0 );
		TargetFn_meanInth<clipper::data32::I_sigI> target_fo( Sigma, 1);
		clipper::ResolutionFn resol_fn( hklinf, basis_fo, target_fo, params );
		
		calcOriginPeak(basis_fo, resol_fn);
		fitOriginPeak(basis_fo, resol_fn);
		
		return optRes();
	}
	
	float PattPeak::optRes()
	{
		float width_res = 0.715*1.0f/_maxinvres;
		float width_patt = 2.0f*_sigma;
		
		return std::sqrt((width_patt*width_patt+width_res*width_res)/2.0f);
	}
	
	// calculate the patterson origin peak in 1-d using fitted data
	
	void PattPeak::calcOriginPeak(clipper::BasisFn_spline& basis_fo, clipper::ResolutionFn& resol_fn) 
	/* -----------------------------------------------------------
	 
	 <I(s)> = average intensity , dS = DETRH
	 
	 P(r) = 4pi * Int <I> * (sin(2pisr)/(2pisr)) * s^2 * ds
	 
	 -----------------------------------------------------------*/
	{
		float widthd = _width/float(_nbins);
		float widthr = _maxinvres/float(_nbins);
		
		for (int id = 0 ; id != _nbins ; ++id ) {
			float d = widthd*(float(id)+0.5f);
			
			for ( int ir=0; ir!=_nbins; ++ir ) {
				float res = widthr*(float(ir)+0.5f); 
				float rsq = res*res;
				float intensity = basis_fo.f_s( rsq, resol_fn.params() );
				float sr = clipper::Util::twopi()*res*d;
				
				_patterson[id] += 2.0f*intensity * std::sin(sr)*rsq*widthr/(res*d);			}
		}
		return;
	}
	
	void PattPeak::fitOriginPeak(clipper::BasisFn_spline& basis_fo, clipper::ResolutionFn& resol_fn) 
	/* fit gaussain to OriginPeak
	 */
	{
		std::vector<clipper::ftype> weights(_nbins,1.0f);
		std::vector<clipper::ftype> x(_nbins);
		std::vector<clipper::ftype> y(_nbins);
		
		for( int i = 0 ; i != _nbins ; ++i) {
			float dist = (float(i)+0.5f)*_width/float(_nbins);
			x[i] = 0.25*dist*dist;
			y[i] = std::log(1.0f/_patterson[i]);
		}
		
        clipper::ftype a, b, siga, sigb;
		
		straight_line_fit(x,y,weights,_nbins,a,b,siga,sigb);
		
		_P000 = std::exp(-b);
		
		float btemp = 0.25*a;
		
		_sigma = std::sqrt(1.0f/std::abs(2.0f*btemp) );
		
		return;	
	}
	
	
	//******Completeness***************************************************************************
	
	template<class D> void Completeness<D>::operator() (clipper::HKL_data<D>& isig, clipper::Resolution reso )
	{
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		
		_t = type(isig[isig.first()]); //set our type
		
		//if resolution not set use observed
		_reso = reso;
        if (_reso.is_null() ) _reso = isig.hkl_info().resolution();
		
		// calculate completeness
		std::vector<clipper::ftype> sumov(_nbins,0.0);
		std::vector<clipper::ftype> summeas(_nbins,0.0);
		
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			// bin number different in C because arrays start at zero
			float mult = 1.0/ih.hkl_class().epsilonc();
			int bin = int( double(_nbins) * ih.invresolsq() / _reso.invresolsq_limit() - 0.001);
			//if (bin >= _nbins || bin < 0) printf("Warning: (completeness) illegal bin number %d\n", bin);
			if ( bin < _nbins && bin >= 0 ) {
				sumov[bin] += mult;
				if ( !isig[ih].missing() ||
				!clipper::Util::is_nan(obs(isig[ih])) || !clipper::Util::is_nan(sigobs(isig[ih])) ) {
					_compsig[bin] += mult;
					clipper::ftype Ival=obs(isig[ih]);
					clipper::ftype Isig=sigobs(isig[ih]);
					clipper::ftype isigi = Ival/Isig;
					if (isigi >= 1.0f ) _compsig1[bin] += mult;
					if (isigi >= 2.0f ) _compsig2[bin] += mult;
					if (isigi >= 3.0f ) _compsig3[bin] += mult;
					_standard[bin] += (isigi > 0.01 ) ? mult/(isigi) : mult/100.0 ;
                                        summeas[bin] += mult;
				}
			}
		}
		for (int i=0; i !=_nbins; ++i) {
			if (sumov[i] > 0.0) {
				_compsig[i] /= sumov[i];
				_compsig1[i] /= sumov[i];
				_compsig2[i] /= sumov[i];
				_compsig3[i] /= sumov[i];
				_standard[i] /= summeas[i];
			}
		}
	}
		
	template<class D> void Completeness<D>::plot() 
	{
		if (_t == I ) {
			printf("\n$TABLE: Intensity Completeness:\n");
			printf("$GRAPHS");
			printf(": Completeness v resolution:N:1,2,3,4,5,6:\n");
			printf("$$ 1/resol^2 completeness sig1 sig2 sig3 standard$$\n$$\n");
		} else {
			printf("\n$TABLE: Structure factor completeness:\n");
			printf("$GRAPHS");
			printf(": Completeness v resolution:N:1,2,3,4,5,6:\n");
			printf("$$ 1/resol^2 completeness sig1 sig2 sig3 standard$$\n$$\n");
		}
		
		for(int i=0;i!=_nbins;++i){
			double res = _reso.invresolsq_limit()*(double(i)+0.5)/double(_nbins);
			printf("%10.6f %8.4f %8.4f %8.4f %8.4f %8.4f\n",res, _compsig[i],_compsig1[i],_compsig2[i],_compsig3[i]
				   ,_standard[i]);
		}
		printf("$$\n\n");
	}
	
	
	//-------class tNCS-------------------------------------------------
	//-------tNCS peak search-------------------------------------------
	
	/*template <class T, template<class> class D> const std::vector<clipper::Symop>& tNCS::operator() (clipper::HKL_data<D<T> >& hkldata, clipper::Resolution r)
	{
		_base = &hkldata;
		_reso = r;
		
		clipper::HKL_info hklinf(hkldata.hkl_info());
		clipper::Cell cell = hklinf.cell();
		clipper::Spacegroup spgr = hklinf.spacegroup();
		
		// check for pseudo translation (taken from cpatterson)
		// get Patterson spacegroup
		clipper::HKL_info hklp;
		clipper::Spacegroup
		pspgr( clipper::Spgr_descr( spgr.generator_ops().patterson_ops() ) );
		hklp.init( pspgr, cell, reso, true );
		
		// make patterson coeffs
		clipper::HKL_data<clipper::datatypes::F_phi<T> > fphi( hklp );
		for ( clipper::HKL_data_base::HKL_reference_index  ih = fphi.first(); !ih.last(); ih.next() ) {
			D<T> i = hkldata[ih.hkl()];
			if ( !i.missing() ) {
				fphi[ih].f() = I(i);
				fphi[ih].phi() = 0.0 ;
			}
		}
		
		// make grid if necessary
		clipper::Grid_sampling grid( pspgr, cell, reso );
		
		// make xmap
		clipper::Xmap<float> patterson( pspgr, cell, grid );
		patterson.fft_from( fphi );
		
		//peak search in patterson
		PeakSearch pksch;                      // peak search object
		
		int npeak = 10;
		
		const std::vector<int>& ppks = pksch( patterson );
		
		if (ppks.size() == 0 || ppks.size() == 1 ) return peaks;
		
		clipper::ftype rho0 = pksch.zero();
		clipper::ftype top_peak = patterson.get_data( ppks[0] ) - rho0;
		int i = 0;
		clipper::ftype pval(0.0);
		
		for (int i = 1 ; i != ppks.size() ; ++i ) {
			clipper::ftype next_peak = patterson.get_data( ppks[i] ) - rho0;
			clipper::Coord_frac c0 = patterson.coord_of( ppks[i] ).coord_frac(grid);
			clipper::ftype ratio = next_peak/top_peak;
			clipper::ftype dist2 = std::sqrt(c0.lengthsq(cell) );
			// look for peaks > 20% of origin peak and at least 14A distant from origin
			// precentage estimate is Zwartz CCP4 Newsletter 42
			if (dist2 > 14.0 ) {
				const clipper::ftype aval = 0.0679;
				const clipper::ftype bval = 3.56;
				pval = (1.0 - std::exp(-std::pow(ratio/(aval*(T(1.0)-ratio)),-bval)) )*100.0;
				if (pval < 1.0) {
					//clipper::Rtop<clipper::ftype> tmp;
					_peaks.push_back(clipper::Symop(clipper::RTop_frac(clipper::Mat33<clipper::ftype>::identity(), c0) ) );
					_peak_prob.push_back(pval);
					_peak_height.push_back(ratio);
				}
			}	
		} 
		
		return _peaks;
	} */
	
	//-------tNCS ouput to stdout----------------------------------------
	void tNCS::output() const
	{
		printf("\n\nTRANSLATIONAL NCS:\n\n");
		if ( _peaks.size() && _peak_prob[0] < 1.0 ) {
			clipper::Vec3<clipper::ftype> c0 = _peaks[0].trn();
			printf("Translational NCS has been detected at (%6.3f, %6.3f, %6.3f).\n  The probability based on peak ratio is %5.2f%% that this is by chance (with resolution limited to %5.2f A). \n", c0[0],c0[1],c0[2],_peak_prob[0],_reso.limit() );
			printf("This will have a major impact on the twinning estimates and effectiveness of the truncate procedure\n\n");
			std::cout << "Peak#       Location       Ratio  Q-score" << std::endl;
			for (int i = 0; i != _peaks.size() ; ++i ) {
				printf("%2d    (%5.2f,%5.2f,%5.2f) %5.2f  %5.2f\n",i+1,(_peaks[i].trn() )[0],(_peaks[i].trn() )[1],(_peaks[i].trn() )[2],_peak_height[i],_peak_prob[i] );
			}
		} else {
			printf("No translational NCS detected (with resolution limited to %5.2f A)\n\n", _reso.limit() );
            if (_peaks.size() ) {
                clipper::Vec3<clipper::ftype> c0 = _peaks[0].trn();
                printf("The highest peak located at (%6.3f, %6.3f, %6.3f).\n  The probability based on peak ratio is %5.2f%% that this is by chance (with resolution limited to %5.2f A). \n", c0[0],c0[1],c0[2],_peak_prob[0],_reso.limit() );
			}
        }
        printf("The analysis uses the peak heights in the patterson map that are further than 14 A (approx. 4 Ca-Ca) from the origin.  The presence of a large off origin peak (above 20%%) and/or a very low Q-score, below 1.0, is a string indicator of the presence of tNCS.  An intermidiate Q-score, between 5.0 and 1.0, may indicate weak tNCS or be the result of cross vector of a large scatterer such as a cluster or heavy metal.\n\n");
        printf("Reference: P. Zwarts CCP4 Newsletter 42\n\n");
    }
	
	//-------tNCS xml to stringstream-----------------------------------
	std::stringstream& tNCS::xml_output(std::stringstream& ss) const {
		ss.precision(3);
		ss << "<translationalNCS>" << std::endl;
		ss << "  <detected>" << ((hasNCS()) ? "Yes" : "No") << "</detected>" << std::endl;
		if ( _peaks.size() && _peak_prob[0] < 1.0 ) {
			clipper::Vec3<clipper::ftype> c0 = _peaks[0].trn();
			ss << "  <Comment id=\"tncs\">" << std::endl;
			ss << "Translational NCS has been detected at (" << std::setw(6) << c0[0] << "," << c0[1] << "," <<c0[2] <<")";
			ss << " .The probability that this is by chance is " << std::setw(5) << std::setprecision(2) << _peak_prob[0] << std::endl;
			ss << "This will have a major impact on the twinning estimates and effectiveness of the truncate procedure." << std::endl;
			ss << "</Comment>";
			ss << "  <PeakList>" << std::endl;
			ss.precision(3);
			for (int i = 0; i != _peaks.size() ; ++i ) {
				ss << "  <Peak>" << std::endl;
				ss << "    <number>" << i+1 << "</number>" << std::endl;
				ss << "    <coordinate type=\"map\" separator=\" \">" << std::setw(5) << (_peaks[i].trn() )[0] << " " << (_peaks[i].trn() )[1] << " " << (_peaks[i].trn() )[2] << "</coordinate>" << std::endl;
				ss << "    <height>" << std::setw(5) << _peak_height[i] << "</height>" << std::endl;
				ss << "    <Q-score>" << std::setw(5) << _peak_prob[i] << "</Q-score>" << std::endl;
				ss << "  </Peak>" << std::endl;
			}
			ss << "  </PeakList>" << std::endl;
		} else {
			ss << "  <Comment id=\"tncs\">" << std::endl;
			ss << "No translational NCS detected." << std::endl;
            if (_peaks.size() ) {
                clipper::Vec3<clipper::ftype> c0 = _peaks[0].trn();
                ss << "The largest peak is at at (" << std::setw(6) << c0[0] << "," << c0[1] << "," <<c0[2] <<")";
                ss << " .The probability that this is by chance is " << std::setw(5) << std::setprecision(2) << _peak_prob[0] << std::endl;
            }
			ss << "  </Comment>" << std::endl;
		}
        ss << "  <Comment id=\"tncsExplain\">" << std::endl;
        ss << " The analysis uses the peak heights in the patterson map that are further than 14 A (approx. 4 Ca-Ca) from the origin.  The presence of a large off origin peak (above 20%) and/or a very low Q-score, below 1.0, is a string indicator of the presence of tNCS.  An intermidiate Q-score, between 5.0 and 1.0, may indicate weak tNCS or be the result of cross vector of a large scatterer such as a cluster or heavy metal." << std::endl;
        ss << "  </Comment>" << std::endl;
        ss << "<Reference id=\"tncsExplain\">";
        ss << "Reference: P. Zwarts CCP4 Newsletter 42";
        ss << "</Reference>" << std::endl;
		ss << "</translationalNCS>" << std::endl;
		return ss;
	}
	
    
	//---------AnomStats--------------------------------------
    // run on scaled on merged data, however the best stats will be from the unmerged set.
    // get from aimless
    AnomStats::AnomStats(const clipper::HKL_data_base& hkldata)
    {
		_base = const_cast<clipper::HKL_data_base *>(&hkldata);
		
		_binner  = ResolStats_base(hkldata);
		_meas    = AnomStats_measurability(_binner);
		_bij     = AnomStats_bijveot(_binner);
		_signoise= AnomStats_signoise(_binner);

        int nbins = _binner.size();
		
		clipper::ftype MEAS_STAT = 0.05;
		clipper::ftype SIG_STAT = (is_intensity() ) ? 1.3 : 1.3;
		clipper::ftype BIJ_STAT = (is_intensity() ) ? 0.010 : 0.006;
		
        //assume values decrease monatomically.  Cut at measurability 5%
        //Dauter Acta D62 (2006) 867
        //Zwart Acta D61 (2005) 1437
        for (int i1 = 0 ; i1 != nbins ; ++i1 ) {
            if ( _meas[i1] > MEAS_STAT ) _meas_range.include(_binner[i1]);
        }
        
        //assume values decrease monatomically.  Cut at DeltaAnom at 1.3
        //Dauter Acta D62 (2006) 867
        //Schneider Acta D58 (2002) 1772
        for (int i1 = 0 ; i1 != nbins ; ++i1 ) {
            if ( _signoise[i1] > SIG_STAT ) _signoise_range.include(_binner[i1]);
        }
        
        //assume values decrease monatomically.  Cut at deltaF/F 0.6%
        //Zwart Acta D61 (2005) 1437
        //Wang Methods Enzymol 115 (1985) 90
        for ( int i1 = 0 ; i1 != nbins ; ++i1 ) {
            if ( _bij[i1] > BIJ_STAT ) _bij_range.include(_binner[i1]);
        }
        
	}
	
	//void AnomStats::output()
	void AnomStats::output()
	{
        std::cout << "Estimated limits of anomalous signal:" << std::endl << std::endl;
        
        if (is_intensity() )  {       
			std::cout << "  Bijvoet ratio (<deltaI>/<I>) > 1.0%              : ";
			if (_bij_range.min() < _bij_range.max() ) 
				std::cout << 1.0/std::sqrt(_bij_range.min() ) << " - " << 1.0/std::sqrt(_bij_range.max() ) << " A " << std::endl;
			else 
				std::cout << "NaN - NaN" << std::endl;
			std::cout << "  Anomalous signal to noise (<deltaI/sigdI>) > 1.3 : ";
			if (_signoise_range.min() < _signoise_range.max() ) 
				std::cout << 1.0/std::sqrt(_signoise_range.min() ) << " - " << 1.0/std::sqrt(_signoise_range.max() ) << " A " << std::endl;
			else 
				std::cout << "NaN - NaN" << std::endl;
			std::cout << "  Measurability limit (Nanon/Nov) > 5%             : ";
			if (_meas_range.min() < _meas_range.max() ) 
				std::cout << 1.0/std::sqrt(_meas_range.min() ) << " - " << 1.0/std::sqrt(_meas_range.max() ) << " A " << std::endl;
			else 
				std::cout << "NaN - NaN" << std::endl;
		} else {
			std::cout << "  Bijvoet ratio (<deltaF>/<F>) > 0.6%              : ";
			if (_bij_range.min() < _bij_range.max() ) 
				std::cout << 1.0/std::sqrt(_bij_range.min() ) << " - " << 1.0/std::sqrt(_bij_range.max() ) << " A " << std::endl;
			else 
				std::cout << "NaN - NaN" << std::endl;
			std::cout << "  Anomalous signal to noise (<deltaF/sigdF>) > 1.3 : " ;
			if (_signoise_range.min() < _signoise_range.max() ) 
				std::cout << 1.0/std::sqrt(_signoise_range.min() ) << " - " << 1.0/std::sqrt(_signoise_range.max() ) << " A " << std::endl;
			else 
				std::cout << "NaN - NaN" << std::endl;
			std::cout << "  Measurability limit (Nanon/Nov) > 5%             : ";
			if (_meas_range.min() < _meas_range.max() ) 
				std::cout << 1.0/std::sqrt(_meas_range.min() ) << " - " << 1.0/std::sqrt(_meas_range.max() ) << " A " << std::endl;
			else 
				std::cout << "NaN - NaN" << std::endl;
		}
        std::cout << std::endl;
        if ( _signoise_range.min() > _signoise_range.max() || _meas_range.min() > _meas_range.max() ) {
            std::cout << "  Warning: NO anomalous signal." << std::endl;
        } else if ( ( _signoise_range.min() < 0.01 && _signoise_range.max() > 0.11 ) && // low resolution lower than 10A, high beyond 3 A
            ( _meas_range.min() < 0.01 && _meas_range.max() > 0.11 ) ) {
            std::cout << " Significant anomalous signal to high resolution." << std::endl;
        } else if ( ( _signoise_range.min() < 0.01 && _signoise_range.max() > 0.04 ) && // low resolution lower than 10A, high beyond 5 A
                   ( _meas_range.min() < 0.01 && _meas_range.max() > 0.04 ) ) {
            std::cout << " Significant anomalous signal to medium resolution." << std::endl;
        } else if ( ( _signoise_range.min() < 0.01 && _signoise_range.max() > 0.04 ) || // low resolution lower than 10A, high beyond 5 A
                   ( _meas_range.min() < 0.01 && _meas_range.max() > 0.04 ) ) {
            std::cout << " Inconsistent results, there may be some anomalous signal." << std::endl;
        } else if ( _signoise_range.min() > 0.04 || _meas_range.min() > 0.04 )  {
            std::cout << " Warning: NO anomalous signal. Anomalous signal is probably artefact of noise at high resolution." << std::endl;
        } else {
            std::cout << " Some anomalous signal at intermediate resolution.  Low resolution signal missing." << std::endl;
        }

		std::cout << std::endl;
        std::cout << "  These calculations are performed using scaled and merged data.  More accurate estimates of the limit of the anomalous signal\n  can be obtained using scaled and unmerged data in the half dataset correlation calculation of aimless. " << std::endl << std::endl;
	//}
		std::cout << "  The measurability is defined as the fraction of the anomalous differences for which the signal to noise for deltaI,\n  I(+) and I(-) are all > 3.";
        std::cout << " The resolution limits are derived from I/sigI of 1.2 for the signal to noise, and measurability above 5%." << std::endl;
		std::cout << "  For well processed data having more than 5% of the measured data satisfying this criteria is a good\n  indicator, particularly when combined with";
		std::cout << " a significate mean signal to noise of the anomalous difference in a\n  resolution shell." << std::endl; 
	
	//void AnomStats::output_loggraph()
	//{
		if (is_intensity() )  {
			printf("\n$TABLE: Intensity anomalous analysis:\n");
			printf("$GRAPHS");
			printf(": Mn(dI/sigdI) v resolution:N:1,2:\n");
			printf(": Mn(dI/I) v resolution:N:1,3:\n");
			printf(": Mesurability v resolution:N:1,4:\n");
			printf("$$ 1/resol^2 Mn(dI/sigdI)) Mn(dI/I) measurability$$\n$$\n");			
		} else {
			printf("\n$TABLE: Structure factor anomalous analysis:\n");
			printf("$GRAPHS");
			printf(": Mn(dF/sigdF) v resolution:N:1,2:\n");
			printf(": Mn(dF/F) v resolution:N:1,3:\n");
			printf(": Mesurability v resolution:N:1,4:\n");
			printf("$$ 1/resol^2 Mn(dF/sigdF)) Mn(dF/F) measurability$$\n$$\n");
		}			
		int nbins = _binner.size();
		for(int i=0;i!=nbins;++i){
			printf("%10.6f %12.4e %12.4e %12.4e\n",_binner[i],_signoise[i],_bij[i],_meas[i]);
		}
		printf("$$\n\n");
    }
	
	std::stringstream& AnomStats::xml_output(std::stringstream& ss)
	{
		std::string is = (is_intensity() ) ? "intensity" : "amplitude" ;
		ss << "<AnomStatistics type=\"" << is << "\">" << std::endl;
		ss << "  <ResolutionRange id=\"Signal to Noise Ratio\" unit=\"Angstrom\" >" << std::endl;
		if (_signoise_range.min() < _signoise_range.max() ) 
			ss << std::fixed << std::setw(5) << std::setprecision(2) << "    <min>" << 1.0/std::sqrt(_signoise_range.min() ) << "</min>\n    <max>"
			<< 1.0/std::sqrt(_signoise_range.max() ) << "</max>" << std::endl;
		else 
			ss <<"    <min> NaN</min>\n    <max> NaN</max>" << std::endl;
		ss << "  </ResolutionRange>" << std::endl;
		ss << "  <ResolutionRange id=\"Measurability Limit\" unit=\"Angstrom\" >" << std::endl;
		if (_meas_range.min() < _meas_range.max() ) 
			ss << std::fixed << std::setw(5) << std::setprecision(2) << "    <min>" << 1.0/std::sqrt(_meas_range.min() ) << "</min>\n    <max>"
			<< 1.0/std::sqrt(_meas_range.max() ) << "</max>" << std::endl;
		else 
			ss <<"    <min> NaN</min>\n    <max> NaN</max>" << std::endl;
		ss << "  </ResolutionRange>" << std::endl;
		ss << "  <Comment id='AnomDescription'>" << std::endl;
		ss << "  These calculations are performed using scaled and merged data.  More accurate estimates of the limit of the anomalous signal\n  can be obtained using scaled and unmerged data in the half dataset correlation calculation of aimless. " << std::endl << std::endl;
		ss << " The measurability is defined as the fraction of the anomalous differences for which the signal to noise for deltaI,\n  I(+) and I(-) are all > 3.";
        ss << " The resolution limits are derived from I/sigI of 1.2 for the signal to noise, and measurability above 5%";
		ss << " For well processed data having more than 5% of the measured data satisfying this criteria is a good\n   indicator, particularly when combined with";
		ss << " a significate mean signal to noise of the anomalous difference in a\n  resolution shell." << std::endl;
        ss << "  </Comment>" << std::endl;
        ss << "  <Comment id='AnomResult'>" << std::endl;
        if ( _signoise_range.min() > _signoise_range.max() || _meas_range.min() > _meas_range.max() ) {
            ss << "  Warning: NO anomalous signal." << std::endl;
        } else if ( ( _signoise_range.min() < 0.01 && _signoise_range.max() > 0.11 ) && // low resolution lower than 10A, high beyond 3 A
            ( _meas_range.min() < 0.01 && _meas_range.max() > 0.11 ) ) {
            ss << "  Significant anomalous signal to high resolution." << std::endl;
        } else if ( ( _signoise_range.min() < 0.01 && _signoise_range.max() > 0.04 ) && // low resolution lower than 10A, high beyond 5 A
            ( _meas_range.min() < 0.01 && _meas_range.max() > 0.04 ) ) {
            ss << "  Significant anomalous signal to medium resolution." << std::endl;
        } else if ( ( _signoise_range.min() < 0.01 && _signoise_range.max() > 0.04 ) || // low resolution lower than 10A, high beyond 5 A
                   ( _meas_range.min() < 0.01 && _meas_range.max() > 0.04 ) ) {
            ss << "  Inconsistent results, there may be some anomalous signal." << std::endl;
        } else if ( _signoise_range.min() > 0.04 || _meas_range.min() > 0.04 )  {
            ss << "  Warning: NO anomalous signal. Anomalous signal is probably artefact of noise at high resolution." << std::endl;
        } else {
            ss << "  Some anomalous signal at intermediate resolution.  Low resolution signal missing." << std::endl;
        }
		ss << "  </Comment>" << std::endl;
		ss << "</AnomStatistics>" << std::endl;
		/*
		 <CCP4Table groupID="graph" id"I-anom" title="Intensity anomalous analysis">
		 <plot>
		 <title>Resolution plots for anomalous scattering indicators</title>
		 <yrange min="0" max="None"/>
		 <plotline xcol="  1" ycol="  4" >
		 <symbolsize>  0</symbolsize>
		 <linestyle>-</linestyle>
		 <colour>red</colour>
		 </plotline>
		 <plotline xcol="  1" ycol="  3" >
		 <symbolsize>  0</symbolsize>
		 <linestyle>-</linestyle>
		 <colour>green</colour>
		 </plotline>
		 <plotline xcol="  1" ycol="  2" >
		 <linestyle>-</linestyle>
		 <colour>black</colour>
		 </plotline>
		 </plot>
		 <headers separator=" ">
		 1/resol^2 Mn(dF/sigdF)) Mn(dF/F) measurability
		 </headers>
		 <data>
		 ... ... ... ...
		 </data>
		 </CCP4Table>
*/
		int nbins = _binner.size();
		
		ss << "<CCP4Table groupID=\"graph\" id=\"anomalous " << is << " plot\" title=\""<< is << " anomalous analysis plots\">" << std::endl;
		ss << "<plot>" << std::endl;
		ss << "<title>Resolution plots for anomalous scattering indicators</title>" << std::endl;
		ss << "<xscale>oneoversqrt</xscale>" << std::endl;
		ss << "<yrange min=\"0\" max=\"None\"/>" << std::endl;
		ss << "<plotline xcol=\"1\" ycol=\"4\" >" <<  std::endl;
		ss << "<symbolsize>  0</symbolsize>" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>red</colour>" << std::endl;
		ss << "</plotline>" << std::endl;
		ss << "<plotline xcol=\"1\" ycol=\"3\" >" << std::endl;
		ss << "<symbolsize>  0</symbolsize>" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>green</colour>" << std::endl;
		ss << "</plotline>" << std::endl;
		ss << "<plotline xcol=\"  1\" ycol=\"  2\" >" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>black</colour>" << std::endl;
		ss << "</plotline>" << std::endl;
		ss << "</plot>" << std::endl;
		if (this->is_intensity()) ss << "<headers separator=\" \">\n 1/resol^2 Mn(dI/sigdI)) Mn(dI/I) measurability \n </headers>" << std::endl;
			else ss << "<headers separator=\" \">\n 1/resol^2 Mn(dF/sigdF)) Mn(dF/F) measurability \n </headers>" << std::endl;
		ss << "<data>" << std::endl;;
		for(int i=0;i!=nbins;++i){
			ss << std::setw(6) << _binner[i] << " " << _signoise[i] << " " << _bij[i] << " " << _meas[i] << std::endl;
		}
		ss << "</data>" << std::endl;
		ss << "</CCP4Table>" << std::endl;
		return ss;
	}
	
	//******ResolStats_base************************************************************************
	void ResolStats_base::init(const clipper::HKL_data_base& hkldata, const clipper::Range<clipper::ftype>& range, const bool missing, const int nreflns)
	{
		_base = const_cast<clipper::HKL_data_base *>(&hkldata);
        _range.include(range.max() );
        _range.include(range.min() );
        _s_ord.init( range, 1000 );
		if (missing) {
			//_s_ord.init( hkldata, 1.0 );;
            for (clipper::HKL_data_base::HKL_reference_index ih = hkldata.first_data(); !ih.last(); hkldata.next_data(ih) )
            if ( range.contains(ih.invresolsq()) ) _s_ord.accumulate( ih.invresolsq() );
		} else {
			//_s_ord.init( hkldata.hkl_info(), 1.0);
            for (clipper::HKL_data_base::HKL_reference_index ih = hkldata.hkl_info().first(); !ih.last(); ih.next() )
                if ( range.contains(ih.invresolsq()) ) _s_ord.accumulate( ih.invresolsq() );
		}
        _s_ord.prep_ordinal();
		
		int nbins(0);
        int Nreflections(0);
        if (missing) {
            for (clipper::HKL_data_base::HKL_reference_index ih = hkldata.first_data(); !ih.last(); hkldata.next_data(ih) )
                if ( range.contains(ih.invresolsq()) ) ++Nreflections;
		} else {
            for (clipper::HKL_data_base::HKL_reference_index ih = hkldata.hkl_info().first(); !ih.last(); ih.next() )
                if ( range.contains(ih.invresolsq()) ) ++Nreflections;
            
		}
		//int Nreflections(hkldata.num_obs() );
		//const int nreflns(500);
		{
			//for ( clipper::HKL_data_base::HKL_reference_index ih = hkldata.first(); !ih.last(); ih.next() ) {
			//		 if (!hkldata.missing(ih.index() ) )++Nreflections;
			//}
			if ( nbins == 0 && nreflns != 0 ) {
				nbins = std::max( Nreflections/nreflns , 1);
				//} else if ( nreflns == 0 && nprm2 != 0 ) {
				//nprm = nbins;
			} else {
				//nprm2 = std::max( Nreflections/nreflns , nprm2);
				double np1(nbins+0.499);
				double np2(Nreflections/nreflns);
				double np(std::sqrt(np1*np1*np2*np2/(np1*np1+np2*np2) ) );
				nbins = std::max( int(np), 1 );
			}
		}
		_b_reso.resize(nbins,0.0);
		_b_contains.resize(nbins,0.0);
        _b_range.resize(nbins);
				
        if (missing) {
            for ( clipper::HKL_data_base::HKL_reference_index ih = hkldata.first_data(); !ih.last(); hkldata.next_data(ih) ) {
                clipper::ftype s = ih.invresolsq();
                if ( range.contains(s) ) {
                    clipper::ftype eps = (this->is_intensity() ) ? 1.0/ih.hkl_class().epsilonc() : 1.0/std::sqrt(ih.hkl_class().epsilonc());
                    int bin = clipper::Util::bound( 0,clipper::Util::intf( clipper::ftype(nbins) * _s_ord.ordinal( s ) ), nbins-1 );
                    _b_reso[bin] += eps*s;
                    _b_contains[bin] += eps;
                    _b_range[bin].include(s);
                }
            }
        } else {
            for (clipper::HKL_data_base::HKL_reference_index  ih = hkldata.hkl_info().first(); !ih.last(); ih.next() ) {
                clipper::ftype s = ih.invresolsq();
                if ( range.contains(s) ) {
                    clipper::ftype eps = (this->is_intensity() ) ? 1.0/ih.hkl_class().epsilonc() : 1.0/std::sqrt(ih.hkl_class().epsilonc());
                    int bin = clipper::Util::bound( 0,clipper::Util::intf( clipper::ftype(nbins) * _s_ord.ordinal( s ) ), nbins-1 );
                    _b_reso[bin] += eps*s;
                    _b_contains[bin] += eps;
                    _b_range[bin].include(s);
                }
            }

        }
        
		for (int i=0 ; i != nbins ; ++i) _b_reso[i] /= _b_contains[i];

		return;
	}
	
	void ResolStats_base::init(const ResolStats_base& base)
	{
        _range.include(base._range.max() );
        _range.include(base._range.min() );
		_s_ord = base._s_ord;
		_b_reso = base._b_reso;
        _b_contains = base._b_contains;
        _b_range = base._b_range;
		_base = base._base;
	}
	
	int ResolStats_base::operator()(const clipper::ftype s) const
	{
		int nbins = this->size();
		return clipper::Util::bound( 0,clipper::Util::intf( clipper::ftype(nbins) * _s_ord.ordinal( s ) ), nbins-1 );
	}
	
	clipper::ftype ResolStats_base::operator[](const int index) const
	{
		int nbins = this->size();
		int bin = clipper::Util::bound( 0, index, nbins-1 );
		return _b_reso[bin];
	}
	
	ResolStats_base& ResolStats_base::operator=(const ResolStats_base& orig)
	{
		init(orig);
		return *this;
	}
	
	//******AnomStats_measurability***************************************************************
	
	AnomStats_measurability::AnomStats_measurability(const clipper::HKL_data_base& hkldata) 
	{
		init(hkldata);
		calc(hkldata);
	}
	
	AnomStats_measurability::AnomStats_measurability(const ResolStats_base& base) 
	{
		init(base);
		calc(*(this->parent() ) );
	}
		
	void AnomStats_measurability::calc(const clipper::HKL_data_base& hkldata)
	{
		clipper::ftype MEAS_LIMIT;
		MEAS_LIMIT = (this->is_intensity() ) ? 3.0 : 1.5;
		int nbins = this->size();
				
		_meas.resize(nbins,0.0);
		
		std::vector<clipper::ftype> sumov(nbins,0.0);
		
		clipper::xtype working[hkldata.data_size()];
		
		for ( clipper::HKL_data_base::HKL_reference_index ih = hkldata.first_data(); !ih.last(); hkldata.next_data(ih) ) {
			if (!hkldata.missing(ih.index() ) ) {
				clipper::ftype eps = (this->is_intensity() ) ? 1.0/ih.hkl_class().epsilonc() : 1.0/std::sqrt(ih.hkl_class().epsilonc());
				clipper::ftype s = ih.invresolsq();
				hkldata.data_export(ih.hkl(),working);
				clipper::ftype Ip(working[0]);
				clipper::ftype sp(working[1]);
				clipper::ftype Im(working[2]);
				clipper::ftype sm(working[3]);
				clipper::ftype dI(0.0), ds(1.0);
				if (!clipper::Util::is_nan(Ip) && !clipper::Util::is_nan(Im) ) {
					dI = std::fabs(Ip - Im);
					ds = std::sqrt(sm*sm+sp*sp);
				}
				int bin = this->operator()(s);
				if ( dI/ds >= MEAS_LIMIT && Ip/sp >= MEAS_LIMIT && Im/sm >= MEAS_LIMIT ) {
					_meas[bin] += eps;
				}
				sumov[bin] += eps;
			}
		}
		
		for (int i=0 ; i != nbins ; ++i) _meas[i] /= sumov[i];
	}
	
	clipper::ftype AnomStats_measurability::operator[](const int index) const
	{
		int nbins = this->size();
		int bin = clipper::Util::bound( 0, index, nbins-1 );
		return _meas[bin];
	}	
	
	AnomStats_measurability& AnomStats_measurability::operator=(const AnomStats_measurability& orig)
	{
		init(orig);
		_meas = orig._meas;
		return *this;
	}
	
	//******AnomStats_bijveot*********************************************************************
	
	AnomStats_bijveot::AnomStats_bijveot(const clipper::HKL_data_base& hkldata) 
	{
		init(hkldata);
		calc(hkldata);
	}
	
	AnomStats_bijveot::AnomStats_bijveot(const ResolStats_base& base) 
	{
		init(base);
		calc(*(this->parent() ) );
	}
	
	void AnomStats_bijveot::calc(const clipper::HKL_data_base& hkldata)
	{
		int nbins = this->size();
		
		_meandI.resize(nbins,0.0);
		
		std::vector<clipper::ftype> meanI(nbins,0.0);
		
		clipper::xtype working[hkldata.data_size()];
				
		for ( clipper::HKL_data_base::HKL_reference_index ih = hkldata.first_data(); !ih.last(); hkldata.next_data(ih) ) {
			if (!hkldata.missing(ih.index() ) ) {
				clipper::ftype eps = (this->is_intensity() ) ? 1.0/ih.hkl_class().epsilonc() : 1.0/std::sqrt(ih.hkl_class().epsilonc());
				clipper::ftype s = ih.invresolsq();
				hkldata.data_export(ih.hkl(),working);
				clipper::ftype Ip(working[0]);
				clipper::ftype Im(working[2]);
				int bin = this->operator()(s);
				if (!clipper::Util::is_nan(Ip) && !clipper::Util::is_nan(Im) ) {
					_meandI[bin] += eps*std::fabs(Ip - Im);
					meanI[bin] += eps*0.5*(Ip+Im);
				} else if (!clipper::Util::is_nan(Ip) ) {
					meanI[bin] += eps*Ip;
				} else if (!clipper::Util::is_nan(Im) ) {
					meanI[bin] += eps*Im;
				}
			}
		}
		
		for (int i=0 ; i != nbins ; ++i) _meandI[i] /= meanI[i];
	}
	
	clipper::ftype AnomStats_bijveot::operator[](const int index) const
	{
		int nbins = this->size();
		int bin = clipper::Util::bound( 0, index, nbins-1 );
		return _meandI[bin];
	}	
	
	AnomStats_bijveot& AnomStats_bijveot::operator=(const AnomStats_bijveot& orig)
	{
		init(orig);
		_meandI = orig._meandI;
		return *this;
	}
	
	//******AnomStats_signoise*******************************************************************
	AnomStats_signoise::AnomStats_signoise(const clipper::HKL_data_base& hkldata) 
	{
		init(hkldata);
		calc(hkldata);
	}
	
	AnomStats_signoise::AnomStats_signoise(const ResolStats_base& base) 
	{
		init(base);
		calc(*(this->parent() ) );
	}
	
			
	clipper::ftype AnomStats_signoise::operator[](const int index) const
	{
		int nbins = this->size();
		int bin = clipper::Util::bound( 0, index, nbins-1 );
		return _meandI[bin];
	}	
	
	void AnomStats_signoise::calc(const clipper::HKL_data_base& hkldata)
	{
		int nbins = this->size();
		
		_meandI.resize(nbins,0.0);
		
		std::vector<clipper::ftype> sumov(nbins,0.0);
		
		clipper::xtype working[hkldata.data_size()];
		
		for ( clipper::HKL_data_base::HKL_reference_index ih = hkldata.first_data(); !ih.last(); hkldata.next_data(ih) ) {
			if (!hkldata.missing(ih.index() ) ) {
				clipper::ftype eps = (this->is_intensity() ) ? 1.0/ih.hkl_class().epsilonc() : 1.0/std::sqrt(ih.hkl_class().epsilonc());
				clipper::ftype s = ih.invresolsq();
				hkldata.data_export(ih.hkl(),working);
				clipper::ftype Ip(working[0]);
				clipper::ftype sp(working[1]);
				clipper::ftype Im(working[2]);
				clipper::ftype sm(working[3]);
				int bin = this->operator()(s);
				if (!clipper::Util::is_nan(Ip) && !clipper::Util::is_nan(Im) ) {
					clipper::ftype dI = std::fabs(Ip - Im);
					clipper::ftype ds = std::sqrt(sp*sp+sm*sm);
					_meandI[bin] += eps*dI/ds;
				}
				sumov[bin] += eps;
			}
		}
		
		for (int i=0 ; i != nbins ; ++i) _meandI[i] /= sumov[i];
	}
	
	AnomStats_signoise& AnomStats_signoise::operator=(const AnomStats_signoise& orig)
	{
		init(orig);
		_meandI = orig._meandI;
		return *this;
	}
	
	//******ResoCorrels***************************************************************************
	
	template<class D1, class D2> void ResoCorrel<D1,D2>::operator() (clipper::HKL_data<D1>& isig1, clipper::HKL_data<D2>& isig2, 
																	 clipper::Range<clipper::ftype> reso )
	{
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		
		_t = type(isig1[isig1.first()]); //set our type
		
		//if resolution not set use observed
		_reso = reso;		
		if (reso.range() < 0.0 ) {
			_reso = isig2.hkl_info().invresolsq_range();
		}
		
		clipper::ftype num_symops = isig1.hkl_info().spacegroup().num_symops();
		
		// calculate completeness
		clipper::ftype eof(0.0), ecf(0.0), eoof(0.0), eccf(0.0), eocf(0.0), nf(0.0);
		std::vector<clipper::ftype> v_eof(_nbins,0.0);
		std::vector<clipper::ftype> v_ecf(_nbins,0.0);
		std::vector<clipper::ftype> v_eoof(_nbins,0.0);
		std::vector<clipper::ftype> v_eccf(_nbins,0.0);
		std::vector<clipper::ftype> v_eocf(_nbins,0.0);
		std::vector<clipper::ftype> v_nf(_nbins,0.0);
		
		for ( HRI ih = isig1.first(); !ih.last(); ih.next() ) {
			// bin number different in C because arrays start at zero
			if ( _reso.contains(ih.invresolsq() )) {
				clipper::ftype mult = num_symops/ih.hkl_class().epsilonc();
				
				int bin = int( double(_nbins) * (ih.invresolsq()- _reso.min()) / _reso.range());
				if ( bin < _nbins && bin >= 0 ) {
  					if ( ih.hkl() != clipper::HKL::zero() &&
						(!clipper::Util::is_nan(obs(isig1[ih])) || !clipper::Util::is_nan(sigobs(isig1[ih]) ) ) &&
						(!clipper::Util::is_nan(obs(isig2[ih.hkl()])) || !clipper::Util::is_nan(sigobs(isig2[ih.hkl()]) ) ) ) {
						clipper::ftype ff = obs(isig1[ih]);
						clipper::ftype foof = obs(isig2[ih.hkl()]);
						nf += mult;
						eof += foof*mult;
						eoof += foof * foof*mult;
						ecf += ff*mult;
						eocf += foof* ff*mult;
						eccf += ff * ff*mult;
						v_nf[bin] += mult;
						v_eof[bin] += foof*mult;
						v_eoof[bin] += foof * foof*mult;
						v_ecf[bin] += ff*mult;
						v_eocf[bin] += foof* ff*mult;
						v_eccf[bin] += ff * ff*mult;
					}
                }
			}
		}
 		for (int i=0; i !=_nbins; ++i) {
			_comp[i] = ( v_nf[i] * v_eocf[i] - v_eof[i] * v_ecf[i] )
					/ sqrt ( ( v_nf[i] * v_eccf[i] - v_ecf[i] * v_ecf[i] ) * ( v_nf[i] * v_eoof[i] - v_eof[i] * v_eof[i] ) );
		}
		_cc = ( nf * eocf - eof * ecf ) / sqrt ( ( nf * eccf - ecf * ecf ) * ( nf * eoof - eof * eof ) );
	}
	
	template<class D1,class D2> void ResoCorrel<D1,D2>::plot() 
	{
		printf("\n$TABLE: Resolution dependent correlation coefficient (CC=%5.4f):\n",_cc);
		printf("$GRAPHS");
		printf(": CC v resolution:N:1,2:\n");
		printf("$$ 1/resol^2 CC$$\n$$\n");
		for(int i=0;i!=_nbins;++i){
			double res = _reso.range()*(double(i)+0.5)/double(_nbins)+_reso.min();
			printf("%10.6f %8.4f \n",res, ( clipper::Util::is_nan(_comp[i]) ) ? 0.0 : _comp[i]);
		}
		printf("$$\n\n");
	}
	
	//----Rings analysis----------------------------------------------
	
	bool IceRings_analyse::present() {	
		bool icer = false;
		for ( int i = 0; i != _rings->Nrings(); ++i) if (_rings->Reject(i) ) icer = true;
		
		return icer;
	}
	
	std::string IceRings_analyse::output() {
		
		clipper::ftype maxres = _data->hkl_info().resolution().invresolsq_limit();
		
		std::stringstream ss;
		
		if ( _rings->MeanSSqr(0) <= maxres && _rings->MeanSigI(0) > 0.0f ) {
			ss << "ICE RING SUMMARY:\n\n";
			ss << " reso  ice_ring  mean_I mean_Sigma Estimated_I   Ratio Zscore Completeness Ave_Completeness\n";
			for ( int i = 0; i != _rings->Nrings(); ++i) {
				float reso = _rings->MeanSSqr(i);
				if ( reso <= maxres && _rings->MeanSigI(i) > 0.0f ) {
					float imean = _rings->MeanI(i);
					float sigImean = _rings->MeanSigI(i);
					//float expectedI = exp( log(basis_ice.f_s( reso, Sigma.params()) ) + param_gauss[1]*reso);
					ss << std::fixed << std::setw(5) << std::setprecision(2) << 1.0f/std::sqrt(reso) << " "
					<< ( (_rings->Reject(i) ) ? "  yes  " : "  no   " )
					<< std::setw(10) << imean << " " 
					<< std::setw(10) << sigImean << " " 
					<< std::setw(11) << _ideal_rings.MeanI(i) << "  "
                    << std::setw(6) << imean/_ideal_rings.MeanI(i) << " "
					<< std::setw(6) << (imean-_ideal_rings.MeanI(i))/sigImean << "    "
					<< std::setw(8) << _rings->Comp(i) << " "
					<< std::setw(8) << _comp[i] << std::endl;
					//printf("%6.2f %-10.2f %-10.2f %-10.2f %-6.2f %-6.2f %-6.2f\n",1.0f/std::sqrt(reso),imean,sigImean,_ideal_rings.MeanI(i),
					//	   (imean-_ideal_rings.MeanI(i))/sigImean,_rings->Comp(i),_comp[bin] );
				}
			}
            if ( _wB == NULL ) {
                ss << std::endl << "The ice rings table shows data Z-scores and completeness for ice ring sensitive resolutions in comparison with neighbouring ice-ring insensitive resolutions. Large z-scores and low completeness at give a strong hint to the presence of ice rings. It may be required to exclude these resolution ranges. " << std::endl;
            } else {
                ss << std::endl << " The ice rings table shows data Z-scores and completeness for ice ring sensitive resolutions in comparison with the Wilson B-factor fit. Large z-scores and low completeness at give a strong hint to the presence of ice rings. It may be required to exclude these resolution ranges." << std::endl;
            }
        }
		return ss.str();
	}
	
	std::stringstream& IceRings_analyse::xml_output(std::stringstream& ss) {
        clipper::ftype maxres = _data->hkl_info().resolution().invresolsq_limit();
		ss << "<IceRingsAnalysis>" << std::endl;
		for ( int i = 0; i != _rings->Nrings(); ++i) {
			float reso = _rings->MeanSSqr(i);
			float imean = _rings->MeanI(i);
			float sigImean = _rings->MeanSigI(i);
            if ( reso <= maxres && sigImean > 0.0f ) {
                ss << "  <Ring>" << std::endl;
                ss << "    <Number>" << i+1 << "</Number>" << std::endl;
                ss << "    <Resolution>" << std::fixed << std::setw(5) << std::setprecision(2) << 1.0f/std::sqrt(reso) << "</Resolution>" << std::endl;
                ss << "    <Reject>" << ( (_rings->Reject(i) ) ? "  yes  " : "  no   " ) << "</Reject>" << std::endl;
                ss << "    <Imean>" << std::setw(10) << imean << "</Imean>" << std::endl;
                ss << "    <SigImean>" << sigImean << "</SigImean>" << std::endl;
                ss << "    <Ratio>" << imean/_ideal_rings.MeanI(i) << "</Ratio>" << std::endl;
                ss << "    <Z-score>" << (imean-_ideal_rings.MeanI(i))/sigImean << "</Z-score>" << std::endl;
                ss << "    <Completeness>" << _rings->Comp(i) << "</Completeness>" << std::endl;
                ss << "    <ExpectCompleteness>" << _comp[i] << "</ExpectCompleteness>" << std::endl;
                ss << "  </Ring>" << std::endl;
		}
        }
        ss << "<Comment id='IceRingsAnalysis'>" << std::endl;
        if ( _wB == NULL ) {
            ss << " The ice rings table shows data Z-scores and completeness for ice ring sensitive resolutions in comparison with neighbouring ice-ring insensitive resolutions. Large z-scores and low completeness at give a strong hint to the presence of ice rings. It may be required to exclude these resolution ranges." << std::endl;
        } else {
            ss << " The ice rings table shows data Z-scores and completeness for ice ring sensitive resolutions in comparison with the Wilson B-factor fit. Large z-scores and low completeness at give a strong hint to the presence of ice rings. It may be required to exclude these resolution ranges." << std::endl;
        }
        ss << "</Comment>" << std::endl;
		ss << "</IceRingsAnalysis>" << std::endl;
		return ss;
	}
    
    //----Outlier Rings analysis----------------------------------------------
	
	bool OutlierRings_analyse::present() {
		bool icer = false;
		for ( int i = 0; i != _rings->Nrings(); ++i) if (_rings->Reject(i) ) icer = true;
		
		return icer;
	}
    
    float OutlierRings_analyse::percentage() {
        int nb(0);
        for ( int i = 0; i != _rings->Nrings(); ++i) if (_rings->Reject(i) ) ++nb;
        
        return float(nb)/float(_rings->Nrings() );
    }
	
	std::string OutlierRings_analyse::output() {
		
		clipper::ftype maxres = _data->hkl_info().resolution().invresolsq_limit();
		
		std::stringstream ss;
		
        bool done(false);
        for ( int i = 0; i != _rings->Nrings(); ++i) {
            if (_rings->MeanSSqr(i) <= maxres && _rings->MeanSigI(i) > 0.0f) done = true;
        }
		if ( done ) {
			ss << "OUTLIER RING SUMMARY:\n\n";
            if ( present() ) {
                ss << "Outliers total " << std::fixed << std::setw(5) << std::setprecision(1) << percentage()*100 << "% of the bins." << std::endl << std::endl;
                ss << " reso    mean_I mean_Sigma Estimated_I  Ratio Zscore Completeness Ave_Completeness\n";
                for ( int i = 0; i != _rings->Nrings(); ++i) {
                    if ( _rings->Reject(i) && _rings->MeanSigI(i) > 0.0f ) {
                        float reso = _rings->MeanSSqr(i);
                        float imean = _rings->MeanI(i);
                        float sigImean = _rings->MeanSigI(i);
                        //float expectedI = exp( log(basis_ice.f_s( reso, Sigma.params()) ) + param_gauss[1]*reso);
                        ss << std::fixed << std::setw(5) << std::setprecision(2) << 1.0f/std::sqrt(reso)
                        << std::setw(10) << imean << " "
                        << std::setw(10) << sigImean << " "
                        << std::setw(11) << _ideal_rings.MeanI(i) << " "
                        << std::setw(6) << imean/_ideal_rings.MeanI(i) << " "
                        << std::setw(6) << (imean-_ideal_rings.MeanI(i))/sigImean << " "
                        << std::setw(8) << _rings->Comp(i) << " "
                        << std::setw(8) << _comp[i] << std::endl;
                    }
                }
            } else {
                ss << "No problem resolution rings found." << std::endl;
            }
            if ( _wB == NULL ) {
                ss << std::endl << "The outlier rings table shows data Z-scores and completeness for problem resolution bins in comparison with neighbouring resolutions. Large z-scores and low completeness at give a strong hint to the presence of problems. It may be required to exclude these resolution ranges. " << std::endl;
            } else {
                ss << std::endl << " The outlier rings table shows data Z-scores and completeness for problem resolution bins in comparison with the Wilson B-factor fit. Large z-scores and low completeness at give a strong hint to the presence of problems. It may be required to exclude these resolution ranges." << std::endl;
            }
        }
        return ss.str();
    }
	
	std::stringstream& OutlierRings_analyse::xml_output(std::stringstream& ss) {
        clipper::ftype maxres = _data->hkl_info().resolution().invresolsq_limit();
		ss << "<OutlierRingsAnalysis>" << std::endl;
        ss << "  <present>" << ( ( present() ) ? "  yes  " : "  no   " ) << "</present>" << std::endl;
        ss << "  <percentage>" << std::fixed << std::setw(5) << std::setprecision(1) <<  percentage()*100 << "</percentage>" << std::endl;
		for ( int i = 0; i != _rings->Nrings(); ++i) {
			float reso = _rings->MeanSSqr(i);
			float imean = _rings->MeanI(i);
			float sigImean = _rings->MeanSigI(i);
            if ( _rings->Reject(i) && sigImean > 0.0f ) {
                ss << "  <Ring>" << std::endl;
                ss << "    <Number>" << i+1 << "</Number>" << std::endl;
                ss << "    <Resolution>" << std::fixed << std::setw(5) << std::setprecision(2) << 1.0f/std::sqrt(reso) << "</Resolution>" << std::endl;
                ss << "    <Reject>" << ( (_rings->Reject(i) ) ? "  yes  " : "  no   " ) << "</Reject>" << std::endl;
                ss << "    <Imean>" << std::setw(10) << imean << "</Imean>" << std::endl;
                ss << "    <SigImean>" << sigImean << "</SigImean>" << std::endl;
                ss << "    <Z-score>" << (imean-_ideal_rings.MeanI(i))/sigImean << "</Z-score>" << std::endl;
                ss << "    <Completeness>" << _rings->Comp(i) << "</Completeness>" << std::endl;
                ss << "    <ExpectCompleteness>" << _comp[i] << "</ExpectCompleteness>" << std::endl;
                ss << "  </Ring>" << std::endl;
            }
        }
        ss << "<Comment id='OutlierRingsAnalysis'>" << std::endl;
        if ( _wB == NULL ) {
            ss << " The outlier rings table shows data Z-scores and completeness for problem resolution bins in comparison with neighbouring resolutions. Large z-scores and low completeness at give a strong hint to the presence of problems. It may be required to exclude these resolution ranges." << std::endl;
        } else {
            ss << " The outlier rings table shows data Z-scores and completeness for problem resolution bins in comparison with the Wilson B-factor fit. Large z-scores and low completeness at give a strong hint to the presence of problems. It may be required to exclude these resolution ranges." << std::endl;
        }
        ss << "</Comment>" << std::endl;
		ss << "</OutlierRingsAnalysis>" << std::endl;
		return ss;
	}

	
	//******HKLStats_completeness***************************************************************
	
	HKLStats_completeness::HKLStats_completeness(const clipper::HKL_data_base& hkldata,clipper::ftype val) 
	{
		_val=val;
		init(hkldata);
		calc(hkldata);
	}
	
	HKLStats_completeness::HKLStats_completeness(const ResolStats_base& base, clipper::ftype val)
	{
		_val=val;
		init(base);
		calc(*(this->parent() ) );
	}
	
	void HKLStats_completeness::calc(const clipper::HKL_data_base& hkldata)
	{
		int nbins = this->size();
		
		_completeness.resize(nbins,0.0);
		
		std::vector<clipper::ftype> sumov(nbins,0.0);
		
		clipper::xtype working[hkldata.data_size()];
		for ( clipper::HKL_data_base::HKL_reference_index ih = hkldata.first(); !ih.last(); ih.next() ) {
			clipper::ftype eps = (this->is_intensity() ) ? 1.0/ih.hkl_class().epsilonc() : 1.0/std::sqrt(ih.hkl_class().epsilonc());
			int bin = this->operator()(ih.invresolsq() );
			if (!hkldata.missing(ih.index() ) ) {
				clipper::ftype I(0.0);
				clipper::ftype sig(1.0);
				clipper::ftype s = ih.invresolsq();
				hkldata.data_export(ih.hkl(),working);
				if (this->is_anomalous() ) {
					clipper::ftype Ip(working[0]);
					clipper::ftype sp(working[1]);
					clipper::ftype Im(working[2]);
					clipper::ftype sm(working[3]);
					if (!clipper::Util::is_nan(Ip) && !clipper::Util::is_nan(Im) && sp > 0.0 && sm > 0.0 ) {
						I = eps*0.5*(Ip+Im);
						sig = eps*0.5*sqrt(sm*sm+sp*sp);
					} else if (!clipper::Util::is_nan(Ip) && sp > 0.0 ) {
						I = eps*Ip;
						sig = eps*sp;
					} else if (!clipper::Util::is_nan(Im) && sm > 0.0 ) {
						I = eps*Im;
						sig = eps*sm;
					}
				} else {
					clipper::ftype Ip(working[0]);
					clipper::ftype sp(working[1]);
					if (!clipper::Util::is_nan(Ip) && sp > 0.0 ) {
						I = eps*Ip;
						sig = eps*sp;
					}
				}
				
				if ( I/sig >= _val ) {
					_completeness[bin] += eps;
				}
			}
			sumov[bin] += eps;
		}
		
		for (int i=0 ; i != nbins ; ++i) _completeness[i] /= sumov[i];
	} 
	
	clipper::ftype HKLStats_completeness::operator[](const int index) const
	{
		int nbins = this->size();
		int bin = clipper::Util::bound( 0, index, nbins-1 );
		return _completeness[bin];
	}	
	
	HKLStats_completeness& HKLStats_completeness::operator=(const HKLStats_completeness& orig)
	{
		init(orig);
		_val = orig._val;
		_completeness = orig._completeness;
		return *this;
	}
	
	//******HKLStats_Rstandard***************************************************************
	
	HKLStats_Rstandard::HKLStats_Rstandard(const clipper::HKL_data_base& hkldata,clipper::ftype val) 
	{
		init(hkldata);
		calc(hkldata);
	}
	
	HKLStats_Rstandard::HKLStats_Rstandard(const ResolStats_base& base, clipper::ftype val) 
	{
		init(base);
		calc(*(this->parent() ) );
	}
	
	void HKLStats_Rstandard::calc(const clipper::HKL_data_base& hkldata)
	{
		int nbins = this->size();
		
		_Rstandard.resize(nbins,0.0);
		
		std::vector<clipper::ftype> sumov(nbins,0.0);
		
		clipper::xtype working[hkldata.data_size()];
		for ( clipper::HKL_data_base::HKL_reference_index ih = hkldata.first(); !ih.last(); ih.next() ) {
			if (!hkldata.missing(ih.index() ) ) {
				clipper::ftype I;
				clipper::ftype sig;
				clipper::ftype eps = (this->is_intensity() ) ? 1.0/ih.hkl_class().epsilonc() : 1.0/std::sqrt(ih.hkl_class().epsilonc());
				clipper::ftype s = ih.invresolsq();
				hkldata.data_export(ih.hkl(),working);
				if (this->is_anomalous() ) {
					clipper::ftype Ip(working[0]);
					clipper::ftype sp(working[1]);
					clipper::ftype Im(working[2]);
					clipper::ftype sm(working[3]);
					if (!clipper::Util::is_nan(Ip) && !clipper::Util::is_nan(Im) && sp > 0.0 && sm > 0.0 ) {
						I = eps*0.5*(Ip+Im);
						sig = eps*0.5*sqrt(sm*sm+sp*sp);
					} else if (!clipper::Util::is_nan(Ip) && sp > 0.0 ) {
						I = eps*Ip;
						sig = eps*sp;
					} else if (!clipper::Util::is_nan(Im) && sm > 0.0 ) {
						I = eps*Im;
						sig = eps*sm;
					}
				} else {
					clipper::ftype Ip(working[0]);
					clipper::ftype sp(working[1]);
					if (!clipper::Util::is_nan(Ip) && sp > 0.0 ) {
						I = eps*Ip;
						sig = eps*sp;
					}
				}
				
				int bin = this->operator()(ih.invresolsq() );
				if ( is_intensity() ) {
					_Rstandard[bin] += ((I>0.0) ? sig/(2.0*sqrt(I) ) : 0.0 );
				} else {
					_Rstandard[bin] += sig;
				}
				sumov[bin] += I;
			}
		}
		
		for (int i=0 ; i != nbins ; ++i) _Rstandard[i] /= sumov[i];
	}
	
	clipper::ftype HKLStats_Rstandard::operator[](const int index) const
	{
		int nbins = this->size();
		int bin = clipper::Util::bound( 0, index, nbins-1 );
		return _Rstandard[bin];
	}	
	
	HKLStats_Rstandard& HKLStats_Rstandard::operator=(const HKLStats_Rstandard& orig)
	{
		init(orig);
		_Rstandard = orig._Rstandard;
		return *this;
	}
	
	//---------HKLAnalysis--------------------------------------
    // run on scaled on merged data, however the best stats will be from the unmerged set.
    // get from aimless
	
	clipper::ftype HKLAnalysis::ACCEPTABLE = 0.85;
	
    /*template<class T, template<class> class D> void HKLAnalysis::operator()(const D<T>& hkldata)
    {
		_base = const_cast<clipper::HKL_data_base *>(&hkldata);
		
		_binner  = ResolStats_base(hkldata);
		_completeness[0]    = HKL_completeness(_binner);
		_completeness[1]    = HKL_completeness(_binner,1.0);
		_completeness[2]    = HKL_completeness(_binner,1.5);
		_completeness[3]    = HKL_completeness(_binner,2.0);
		_completeness[4]    = HKL_completeness(_binner,2.5);
		_completeness[5]    = HKL_completeness(_binner,3.0);
		
		for (int ii=0; ii != _completeness.size() ; ++ii) (_completeness[ii])(hkldata);
		
        int nbins = _binner.size();
		
		int NBINS = _binner.size();
		for (int ii=0; ii != _completeness.size() ; ++ii) {
			if ((_completeness[ii])[0]) >= ACCEPTABLE) _activerange.include(hkldata.hkl_info().invresolsq_range().min());
			if ((_completeness[ii])[NBINS-1]) >= ACCEPTABLE) _activerange.include(hkldata.hkl_info().invresolsq_range().max());
			for ( int i=1 ; i != NBINS ; ++i) {
				int i1=i-1;
				if ( (_completeness[ii])[i] >= ACCEPTABLE && (_completeness[ii])[i1]) < ACCEPTABLE ) {
					_activerange[ii].include(0.5*(_binner[i]+_binner[i1]) );
				} else if ( (_completeness[ii])[i] < ACCEPTABLE && (_completeness[ii])[i1]) >= ACCEPTABLE ) {
					_activerange[ii].include(0.5*(_binner[i]+_binner[i1]) );
				}
			}
		}
		
		return;
	}*/
	/*
				
				// try on Istandard
				if (active_range.max() == -999999999 && active_range.min() == 999999999 ) {
					printf("         Attempt resolution range estimate using Istandard < 1.0 (Ideally would use 0.2) \n");
					int i = 0;
					clipper::Range<double> range(reso_range.min(),reso_range.max() );
					for ( ; i != NBINS-1 ; ++i) {
						if ( compt.standard(compt.bin2invresolsq(i)) < 1.0 && compt.standard(compt.bin2invresolsq(i+1)) < 1.0 ) break;
					}
					if ( i != (NBINS-1) ) {
						int j = NBINS-1;
						for ( ; j != 1 ; --j) {
							if ( compt.standard(compt.bin2invresolsq(j)) < 1.0 && compt.standard(compt.bin2invresolsq(j-1)) < 1.0 ) break;
						}
						if (j != 0 )
							if (i != 0) {
								float d = (compt.bin2invresolsq(i)+compt.bin2invresolsq(i-1))/2.0;
								active_range.include(d);
							} else {
								active_range.include(reso_range.min() );
							}
						if (j != NBINS-1 ) {
							float d = (compt.bin2invresolsq(j)+compt.bin2invresolsq(j+1))/2.0;
							active_range.include(d);
						} else {
							active_range.include(range.max() );
						}
					} else {
						active_range = range;
					}
					if ( active_range.max() != -999999999 && active_range.min() != 999999999 ) {
						printf("         Resolution Range of this data is %7.3fA to %7.3fA\n",1.0/std::sqrt(active_range.min() ), 1.0/std::sqrt(active_range.max() ) );
					}
				}
				
	*/	
	const clipper::Range<clipper::ftype>& HKLAnalysis::active_range() {
		if (_active.min() < _active.max() ) return _active;
        
        const clipper::ftype ALIM(0.2), BLIM(0.35), RR(4.0), RLIM(3.0);
        
        clipper::ftype rmax = _data->hkl_info().resolution().limit();
        
        std::vector<clipper::ftype> per(_activerange.size() );
        
        for (int i = 0; i != _activerange.size() ; ++i) per[i] = float(_binner(_activerange[i].max() )-_binner(_activerange[i].min() )+1 )/_binner.size();
        
		int ii(0), i3(0);
        for (int i=0; i != _activerange.size(); ++i) {
            if ( int(_completeness[i].IoversigI() ) == 3 ) i3 = ii = i;
        }
    
        // first i/sigI > 3 limits
        clipper::ftype rmin = 1.0/std::sqrt((_activerange[i3]).min() );
        clipper::ftype amax = 1.0/std::sqrt((_activerange[i3]).max() );
        clipper::ftype rr = std::fabs(rmin-amax);
        if ( ( rr > RR && (amax < RLIM || amax < rmax+0.5 ) && per[i3] > ALIM ) || per[i3] > BLIM )  {
            _active = _activerange[i3];
            --ii;
        }
        
        // second R standard > 0.1
        if (ii == i3) {
            clipper::Range<clipper::ftype> rrange;
            int NBINS = _binner.size();
            clipper::ftype acceptable = 0.1;
            if (_Rstandard[0] <= acceptable && _Rstandard[0] > 0.0) rrange.include(_data->hkl_info().invresolsq_range().min());
            if (_Rstandard[NBINS-1] <= acceptable && _Rstandard[NBINS-1] > 0.0) rrange.include(_data->hkl_info().invresolsq_range().max());
            for ( int i=1 ; i != NBINS ; ++i) {
                int i1=i-1;
                if ( (_Rstandard[i] <= acceptable && _Rstandard[i] > 0.0 ) && _Rstandard[i1] > acceptable ) {
                    rrange.include(_binner[i1]+(acceptable-_Rstandard[i1])/(_Rstandard[i]-_Rstandard[i1])*(_binner[i]-_binner[i1]) );
                } else if ( _Rstandard[i] > acceptable && (_Rstandard[i1] <= acceptable && _Rstandard[i1] > 0.0 ) ) {
                    rrange.include(_binner[i1]+(acceptable-_Rstandard[i1])/(_Rstandard[i]-_Rstandard[i1])*(_binner[i]-_binner[i1]) );
                }
            }
            rmin = 1.0/std::sqrt(rrange.min() );
            amax = 1.0/std::sqrt(rrange.max() );
            rr = std::fabs(rmin-amax);
            clipper::ftype per = float(_binner(rrange.max() )-_binner(rrange.min() ) +1 )/_binner.size();
            if ( ( rr > RR && (amax < RLIM || amax < rmax+0.5 ) && per > ALIM ) || per > BLIM )  {
                _active = rrange;
                --ii;
            }
        }
        
        //next try lower sigma levels
        if (ii == i3) {
            for (; ii != 0 ; --ii) {
                rmin = 1.0/std::sqrt((_activerange[ii]).min() );
                amax = 1.0/std::sqrt((_activerange[ii]).max() );
                rr = std::fabs(rmin-amax);
                if ( ( rr > RR && (amax < RLIM || amax < rmax+0.5 ) && per[ii] > ALIM ) || per[ii] > BLIM )  {
                    _active = _activerange[ii];
                    break;
                }
            }
        }
        
        // reduce acceptable limit in steps down to 35% for I/sigI > 3
        if (ii == 0) {
            clipper::Range<clipper::ftype> rrange;
            int NBINS = _binner.size();
            ii = 10;
            for (; ii != 0 ; --ii) { // try down to 30% completeness
                float acceptable = float(ii)*0.05 + 0.3;
                if ((_completeness[i3])[0] >= acceptable) rrange.include(_data->hkl_info().invresolsq_range().min());
                if ((_completeness[i3])[NBINS-1] >= acceptable) rrange.include(_data->hkl_info().invresolsq_range().max());
                for ( int i=1 ; i != NBINS ; ++i) {
                    int i1=i-1;
                    if ( (_completeness[i3])[i] >= acceptable && (_completeness[i3])[i1] < acceptable ) {
                        rrange.include(_binner[i1]+(acceptable-(_completeness[i3])[i1])/((_completeness[i3])[i]-(_completeness[i3])[i1])*(_binner[i]-_binner[i1]) );
                    } else if ( (_completeness[i3])[i] < acceptable && (_completeness[i3])[i1] >= acceptable ) {
                        rrange.include(_binner[i1]+(acceptable-(_completeness[i3])[i1])/((_completeness[i3])[i]-(_completeness[i3])[i1])*(_binner[i]-_binner[i1]) );
                    }
                }
                rmin = 1.0/std::sqrt(rrange.min() );
                amax = 1.0/std::sqrt(rrange.max() );
                rr = std::fabs(rmin-amax);
                clipper::ftype per = float(_binner(rrange.max() )-_binner(rrange.min() ) +1 )/_binner.size();
                if ( ( rr > RR && (amax < RLIM || amax < rmax+0.5 ) && per > ALIM ) || per > BLIM )  {
                    _active = rrange;
                    break;
                }
            }
            
            // finally just use SHELX style limits
            if (ii == 0) _active = clipper::Range<clipper::ftype>(_data->hkl_info().invresolsq_range().min(),1.0/std::pow(rmax+0.5,2) );
        }
		return _active;	
	}
	
	//void AnomStats::output()
	void HKLAnalysis::output()
	{
        active_range();
        int i3(0), ia(0);
        for (int ii=0; ii != _activerange.size(); ++ii) {
            if ( int(_completeness[ii].IoversigI() ) == 3 ) i3 = ii;
        }
        for (int ii=0; ii != _activerange.size(); ++ii) {
            if ( _active.max() == (_activerange[ii]).max() && _active.min() == (_activerange[ii]).min() ) ia = ii;
        }
		if (this->is_intensity() ) {
			printf("\nCOMPLETENESS ANALYSIS (using intensities):\n");
			printf("\nThe following uses I/sigI Completeness levels, in particular targeting completeness above 85%%.  The Completeness with I/sigma above 3 indicates a strong signal (A better estimate is available using CC1/2 in aimless).\n");
            std::cout << std::endl;
            std::cout << "   I/sigI>N            range(A)      %refln" << std::endl;
        } else {
            printf("\nCOMPLETENESS ANALYSIS (using amplitudes):\n");
			printf("\nThe following uses F/sigF completeness levels.  A better estimate is available using CC1/2.\n");
            std::cout << std::endl;
			std::cout << "   F/sigF>N            range(A)      %refln" << std::endl;
        }
        for (int ii = _activerange.size()-1; ii != 0 ; --ii) {
            clipper::ftype rmax = _data->hkl_info().resolution().limit();
            clipper::ftype rmin = 1.0/std::sqrt((_activerange[ii]).min() );
            clipper::ftype amax = 1.0/std::sqrt((_activerange[ii]).max() );
            std::cout << "    " << std::fixed << std::setprecision(1) << std::setw(4) << float(_completeness[ii].IoversigI()) <<  "        ";
            if ((_activerange[ii]).min() < (_activerange[ii]).max() )
                std::cout << "    " << std::fixed << std::setprecision(2) << std::setw(5) << rmin << " - " << std::setw(5) << amax << "    " << std::fixed << std::setprecision(1) << std::setw(5) << 100.*float(_binner(_activerange[ii].max() )-_binner(_activerange[ii].min() )+1 )/_binner.size() << "  " << (( ii == ia  )  ? "***" : "") << std::endl;
            else
                std::cout << "      NaN -   NaN      0.0" << std::endl;
        }
        std::cout << "     N/A            ";
        if ((_activerange[0]).min() < (_activerange[0]).max() )
            std::cout << std::fixed << std::setw(5) << std::setprecision(2) << 1.0/std::sqrt((_activerange[0]).min() ) << " - "  << std::setw(5) << 1.0/std::sqrt((_activerange[0]).max() ) << std::endl;
        else
            std::cout << "      NaN -   NaN      0.0" << std::endl;
        std::cout << std::endl;
        if ( (_activerange[i3]).max() == -999999999 && (_activerange[i3]).min() == 999999999 ) {
            if (this->is_intensity() ) printf("WARNING: The resolution range with I/sigI > 3");
            else printf("WARNING: The resolution range with F/sigF > 3");
            printf(" and completeness above %4.2f could not be\n",ACCEPTABLE);
            printf("determined.  The completeness of this data is poor. In order not to discard too much data\nthe resolution range for analysing the statistics has been relaxed to, %6.2f - %6.2fA:\nnote the statistics output will be less accurate.\n\n",1.0/std::sqrt(_active.min() ), 1.0/std::sqrt(_active.max() ));
        } else if  (i3 != ia ) {
            if (this->is_intensity() ) printf("WARNING: The resolution range with I/sigI > 3");
            else printf("WARNING: The resolution range with F/sigF > 3");
            printf(" with completeness above\n0.85 is small, %6.2fA to %6.2fA.  This strong data corresponds to\napproximately %3d%% of the reflections in the file.",1.0/std::sqrt((_activerange[i3]).min() ), 1.0/std::sqrt((_activerange[i3]).max() ),int(100*float(_binner(_activerange[i3].max() )-_binner(_activerange[i3].min() )+1 )/_binner.size()) );
            printf("  In order not to\ndiscard too much data the resolution range for analysing the statistics\nhas been relaxed to, %6.2f - %6.2fA (%3d%%): note the statistics output will be less accurate.\n\n",1.0/std::sqrt(_active.min() ), 1.0/std::sqrt(_active.max() ), int(100*float(_binner(_active.max() )-_binner(_active.min() )+1 )/_binner.size()) );
            if (_activerange[i3].min() > 0.0011 ) printf("Strong data does not extend to low resolutions.\n"); // low resolution higher than 30A
            if (_activerange[i3].max() < 0.11 ) printf("Strong data does not extend to high resolutions.\n"); // high resolution lower than 3A
            if ( ( _activerange[0].max() - _activerange[2].max() ) > 0.1 ) printf("Processed data contains a lot of weak data at high resolution.  Moments plots and Wilson plot should be consulted\n"); // large amount of weak data beyond limit.
            std::cout << std::endl;
            
            
        } else {
            if (this->is_intensity() ) printf("The resolution range with I/sigI > 3");
            else printf("The resolution range with F/sigF > 3");
            printf(" with completeness above %4.2f, the estimated strong data resolution range ",ACCEPTABLE);
            printf("of this data, is %6.2fA to %6.2fA.\n  This corresponds to approximately %3d%% of the reflections in the file.\n\n",1.0/std::sqrt(_active.min() ), 1.0/std::sqrt(_active.max() ), int(100*float(_binner(_activerange[i3].max() )-_binner(_activerange[i3].min() )+1 )/_binner.size()));
        }
    
        if (is_intensity() )  {
			printf("\n$TABLE: Intensity Completeness analysis:\n");
			printf("$GRAPHS");
			printf(": Completeness & (I/sigI)>N v resolution:N:1,2,3,4,5,6,7,8:\n");
			printf(": Completeness & Rstandard v resolution:N:1,2,9:\n");
			printf("$$ 1/resol^2 Completeness (I/s>15) (I/s>10) (I/s>5) (I/s>3) (I/s>2) (I/s>1) Rstandard$$\n$$\n");
		} else {
			printf("\n$TABLE: Structure factor Completeness analysis:\n");
			printf("$GRAPHS");
			printf(": Completeness & %%(F/sigF)>N v resolution:N:1,2,3,4,5:\n");
			printf("$$ 1/resol^2 Compl (I/s>15) (I/s>10) (I/s>5) (F/s>3) (F/s>2) (F/s>1) Rstandard$$\n$$\n");
		}			
		int nbins = _binner.size();
		for(int i=0;i!=nbins;++i){
			printf("    %5.4f   %5.3f        %5.3f    %5.3f    %5.3f   %5.3f   %5.3f   %5.3f   %6.3f\n",_binner[i],_completeness[0][i],_completeness[6][i],_completeness[5][i],_completeness[4][i],_completeness[3][i],_completeness[2][i],_completeness[1][i],_Rstandard[i]);
		}
		printf("$$\n\n");
		
        std::cout << "The completeness at various resolution limit plots gives the completeness after applying a I/sigI cutoff.  The profiles give an indication of the quality of the data.  The Rstandard plot (<sigF>/<F>) gives an alternative indicator.  Strongly recorded resolution bins would typically have values below 0.1." << std::endl;
        
        if (is_intensity() )  printf("\nLow Resolution Intensity Completeness analysis:\n");
        else printf("\nLow Structure factor Completeness analysis:\n");
        printf("   1/resol^2    Range         Completeness \n");
		nbins = _lb.size();
		for(int i=0;i!=nbins;++i){
			printf("%10.4f   %5.3f-%5.3f   %5.3f [%4.1f:%4.1f]\n",_lb[i],(_lcompleteness.binRange(i)).min(),(_lcompleteness.binRange(i)).max(),_lcompleteness[i],_lcompleteness[i]*_lb.binContains(i),_lb.binContains(i)  );
		}
		printf("\n\n");
        
        std::cout << "Low completeness at low resolution can lead to map distortions and other difficulties.  This often arises through experimental effects such as incorrectly alligned crystal, poorly positioned backstop, or over exposure." << std::endl << std::endl << std::endl;
		//completeness with direction
		
		//ice rings
		std::cout << _ira.output() << std::endl;
		
		//wilson plot
		_wilsonB.output();
        
        //outlier rings
		std::cout << _ora.output() << std::endl;
    }
	
	std::stringstream& HKLAnalysis::xml_output(std::stringstream& ss)
	{
        int i3(0),ia(0);
        active_range();
        for (int ii=0; ii != _activerange.size(); ++ii) {
            if ( int(_completeness[ii].IoversigI() ) == 3 ) i3 = ii;
        }
        for (int ii=0; ii != _activerange.size(); ++ii) {
            if ( _active.max() == (_activerange[ii]).max() && _active.min() == (_activerange[ii]).min() ) ia = ii;
        }
		std::string is = (is_intensity() ) ? "intensity" : "amplitude" ;
		ss << "<DataStatistics type=\"" << is << "\">" << std::endl;
		ss << "  <ResolutionRange id=\"Completeness\" unit=\"Angstrom\" >" << std::endl;
		if ( (_activerange[0]).min() < (_activerange[0]).max() ) 
			ss << std::fixed << std::setprecision(2) << "    <min>" << 1.0/std::sqrt((_activerange[0]).min() ) << "</min>\n    <max>"
			<< 1.0/std::sqrt((_activerange[0]).max() ) << "</max>\n    <percentage>100</percentage>" << std::endl;
		else 
			ss <<"    <min> NaN</min>\n    <max> NaN</max>\n    <percentage> NaN</percentage> " << std::endl;
		ss << "  </ResolutionRange>" << std::endl;
        //
        for (int i = _completeness.size()-1; i != 0; --i) {
            ss << "  <ResolutionRange id=\"Completeness" << int(_completeness[i].IoversigI() ) << "\" unit=\"Angstrom\"";
            if (this->is_intensity() ) ss << " title=\"Completeness I/sigmaI &gt; " << _completeness[i].IoversigI() << "\"";
            else ss << " title=\"Completeness F/sigmaF &gt; " << _completeness[i].IoversigI() << "\"";
            ss << ">" << std::endl;
            if ( (_activerange[i]).min() < (_activerange[i]).max() )
                ss << std::fixed << std::setprecision(2) << "    <min>" << 1.0/std::sqrt((_activerange[i]).min() ) << "</min>\n    <max>"
                << std::fixed << std::setprecision(2) << 1.0/std::sqrt((_activerange[i]).max() ) << "</max>\n         <percentage>" <<
                int(100*float(_binner(_activerange[i].max() )-_binner(_activerange[i].min() )+1 )/_binner.size() ) << "</percentage>" << std::endl;
            else
                ss <<"    <min> NaN</min>\n    <max> NaN</max>\n     <percentage> NaN</percentage>" << std::endl;
            ss << "  </ResolutionRange>" << std::endl;
        }
        ss << " <Comment id=\"CompletenessReso\">" << std::endl;
        ss << "The Completeness max and min resolution indicte the resolution extremes where obs/sig(obs) &gt; N accounts for more than 85 % of observations." << std::endl;
        ss << "  </Comment>" << std::endl;
		if (this->is_intensity() ) {
			ss << "  <Comment id=\"Completeness3\" >" << std::endl;
			ss << "The Completeness with I/sigma above 3 indicates a strong signal, " << std::endl;
            if ( (_activerange[i3]).min() < (_activerange[i3]).max() ) {
                ss << "here it is from " << std::fixed << std::setprecision(2) << 1.0/std::sqrt(_activerange[i3].min() ) << " to " << 1.0/std::sqrt(_activerange[i3].max() ) << "A." << std::endl;
            } else {
                ss << " here it could not be determined.  This is very poor data."  << std::endl;
            }
			ss << "(The CC1/2 on unmerged intensities gives, perhaps, the best estimate of the useable data.)" << std::endl;
			ss << "  </Comment>" << std::endl;
		}
		else {
			ss << "  <Comment id=\"Completeness3\" >" << std::endl;
			ss << "The Completeness with F/sigma above 3 indicates very strong data, " << std::endl;
            if ( (_activerange[i3]).min() < (_activerange[i3]).max() ) {
                ss << "here it is from " << std::fixed << std::setprecision(2) << 1.0/std::sqrt(_activerange[i3].min() ) << " to " << 1.0/std::sqrt(_activerange[i3].max() ) << "A (approximately " << int(100*float(_binner(_activerange[i3].max() )-_binner(_activerange[i3].min() )+1) )/_binner.size() << " percent of the reflections)." << std::endl;
            } else {
                ss << " here it is not defined."  << std::endl;
            }
			ss << "The CC1/2 on unmerged intensities gives, perhaps, the best estimate of the useable data." << std::endl;
			ss << "  </Comment>" << std::endl;
		}
        ss << "<Comment id=\"ActiveRange\">" << std::endl;
        if ( (_activerange[i3]).max() == -999999999 && (_activerange[i3]).min() == 999999999 )
            ss << "&lt;li&gt; In order to not discard all the data the resolution range for analysing the statistics has been relaxed to " << std::fixed << std::setprecision(2) << 1.0/std::sqrt(_active.min() ) << " to " << 1.0/std::sqrt(_active.max() ) << "A. &lt;/li&gt;" << std::endl;
        else if (ia == i3 )
            ss << "&lt;li&gt; This will be used in the analysis of the data. &lt;/li&gt;" << std::endl;
        else
            ss << "&lt;li&gt; In order not to discard too much data the resolution range for analysing the statistics has been relaxed to " << std::fixed << std::setprecision(2) << 1.0/std::sqrt(_active.min() ) << " to " << 1.0/std::sqrt(_active.max() ) << "A. &lt;/li&gt;" << std::endl;
        if (_activerange[i3].min() > 0.0011 ) ss << "&lt;li&gt; Strong data does not extend to low resolutions. &lt;/li&gt;" << std::endl; // low resolution higher than 30A
        if (_activerange[i3].max() < 0.11 ) ss << "&lt;li&gt; Strong data does not extend to high resolutions. &lt;/li&gt;" << std::endl; // high resolution lower than 3A
        if ( ( _activerange[0].max() - _activerange[2].max() ) > 0.1 ) ss << "&lt;li&gt; Processed data contains a lot of weak data at high resolution.  Moments plots and Wilson plot should be consulted. &lt;/li&gt;" << std::endl; // large amount of weak data beyond limit.
        ss << "</Comment>" << std::endl;
        ss << "<Comment id=\"CompletenessQuality\">" << std::endl;
        if ( (_activerange[i3]).max() == -999999999 && (_activerange[i3]).min() == 999999999 )
            ss << "Completeness test shows very poor data.";
        else if (ia == i3 && _activerange[i3].min() < 0.0011 && _activerange[i3].max() > 0.11 && ( _activerange[0].max() - _activerange[2].max() ) < 0.1 )
            ss << "Completeness test shows good data." << std::endl;
        else
            ss << "Completeness test shows some issues." << std::endl;
        ss << "</Comment>" << std::endl;
		ss << "</DataStatistics>" << std::endl;
		/*
		 <CCP4Table groupID="graph" id"I-anom" title="Intensity anomalous analysis">
		 <plot>
		 <title>Resolution plots for anomalous scattering indicators</title>
		 <yrange min="0" max="None"/>
		 <plotline xcol="  1" ycol="  4" >
		 <symbolsize>  0</symbolsize>
		 <linestyle>-</linestyle>
		 <colour>red</colour>
		 </plotline>
		 <plotline xcol="  1" ycol="  3" >
		 <symbolsize>  0</symbolsize>
		 <linestyle>-</linestyle>
		 <colour>green</colour>
		 </plotline>
		 <plotline xcol="  1" ycol="  2" >
		 <linestyle>-</linestyle>
		 <colour>black</colour>
		 </plotline>
		 </plot>
		 <headers separator=" ">
		 1/resol^2 Mn(dF/sigdF)) Mn(dF/F) measurability
		 </headers>
		 <data>
		 ... ... ... ...
		 </data>
		 </CCP4Table>
		 */
		int nbins = _binner.size();
		
		ss << "<CCP4Table groupID=\"graph\" id=\"completeness\" title=\""<< is << " completeness vs resolution\">" << std::endl;
		ss << "<plot>" << std::endl;
		if (this->is_intensity() ) ss << "<title>Resolution plots for completeness at various I/sigmaI levels</title>" << std::endl;
		else ss << "<title>Resolution plots for completeness at various F/sigmaF levels</title>" << std::endl;
		ss << "<xscale>oneoversqrt</xscale>" << std::endl;
        ss << "<yrange min=\"0\" max=\"None\"/>" << std::endl;
		ss << "<plotline xcol=\"1\" ycol=\"8\" >" <<  std::endl;
		ss << "<symbolsize>  0</symbolsize>" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>yellow</colour>" << std::endl;
		ss << "</plotline>" << std::endl;
        ss << "<yrange min=\"0\" max=\"None\"/>" << std::endl;
		ss << "<plotline xcol=\"1\" ycol=\"7\" >" <<  std::endl;
		ss << "<symbolsize>  0</symbolsize>" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>blue</colour>" << std::endl;
		ss << "</plotline>" << std::endl;
        ss << "<yrange min=\"0\" max=\"None\"/>" << std::endl;
		ss << "<plotline xcol=\"1\" ycol=\"6\" >" <<  std::endl;
		ss << "<symbolsize>  0</symbolsize>" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>pink</colour>" << std::endl;
		ss << "</plotline>" << std::endl;
		ss << "<yrange min=\"0\" max=\"None\"/>" << std::endl;
		ss << "<plotline xcol=\"1\" ycol=\"5\" >" <<  std::endl;
		ss << "<symbolsize>  0</symbolsize>" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>brown</colour>" << std::endl;
		ss << "</plotline>" << std::endl;		
		ss << "<plotline xcol=\"1\" ycol=\"4\" >" <<  std::endl;
		ss << "<symbolsize>  0</symbolsize>" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>red</colour>" << std::endl;
		ss << "</plotline>" << std::endl;
		ss << "<plotline xcol=\"1\" ycol=\"3\" >" << std::endl;
		ss << "<symbolsize>  0</symbolsize>" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>green</colour>" << std::endl;
		ss << "</plotline>" << std::endl;
		ss << "<plotline xcol=\"  1\" ycol=\"  2\" >" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>black</colour>" << std::endl;
		ss << "</plotline>" << std::endl;
		ss << "</plot>" << std::endl;
		ss << "<plot>" << std::endl;
		ss << "<title>Resolution plots for completeness and Rstandard</title>" << std::endl;
		ss << "<xscale>oneoversqrt</xscale>" << std::endl;
		ss << "<yrange min=\"0\" max=\"None\"/>" << std::endl;
		ss << "<plotline xcol=\"1\" ycol=\"9\" >" <<  std::endl;
		ss << "<symbolsize>  0</symbolsize>" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>red</colour>" << std::endl;
		ss << "</plotline>" << std::endl;
		ss << "<plotline xcol=\"  1\" ycol=\"  2\" >" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>black</colour>" << std::endl;
		ss << "</plotline>" << std::endl;
		ss << "</plot>" << std::endl;
        ss << "<headers separator=\" \">\n 1/resol^2 completeness";
        for (int ii = _completeness.size()-1; ii != 0; --ii) {
            if (this->is_intensity() ) ss << " I/sigI&gt;" << int(_completeness[ii].IoversigI());
            else ss << " F/sigF&gt;" << int(_completeness[ii].IoversigI());
        }
        ss << " Rstandard\n </headers>" << std::endl;
		ss << "<data>" << std::endl;;
		for(int i=0;i!=nbins;++i){
			ss << std::fixed << std::setw(10) << std::setprecision(4) << _binner[i] << " " << std::setw(5) << std::setprecision(3) << _completeness[0][i] << " ";
            for (int ii = _completeness.size()-1; ii != 0; --ii)
                ss << " " << std::fixed << std::setw(5) << std::setprecision(3) << _completeness[ii][i] ;
			ss << " " << std::fixed << std::setw(5) << std::setprecision(3) << _Rstandard[i] << std::endl;
		}
		ss << "</data>" << std::endl;
        ss << "</CCP4Table>" << std::endl;
        
        ss << "<DataStatistics>" << std::endl;
        ss << "<Comment id=\"Completeness\">" << std::endl;
        ss << "The completeness at various resolution limit plots gives the completeness after applying a I/sigI cutoff.  The profiles give an indication of the quality of the data.  The Rstandard plot (&lt;sigF&gt;/ &lt;F&gt;) gives an alternative indicator.  Strongly recorded resolution bins would typically have values below 0.1" << std::endl;
        ss << "</Comment>" << std::endl;
        nbins = _lb.size();
        ss << "<LowResoCompleteness>" << std::endl;
        for(int i=0; i!=nbins;++i) {
            ss << "<bin id=\"LowResoCompleteness\">" << i+1 << "</bin>" << std::endl;
            ss << "  <ResolutionRange id=\"LowResoCompleteness\" unit=\"Angstrom\" >" << std::endl;
            ss << std::fixed << std::setprecision(2) << "    <min>" << 1.0/std::sqrt((_lb.binRange(i)).min() ) << "</min>\n    <max>"
            << std::fixed << std::setprecision(2) <<  1.0/std::sqrt((_lb.binRange(i)).max() ) << "</max>" << std::endl;
            ss << "  </ResolutionRange>" << std::endl;
            ss << "  <Completeness id=\"LowResoCompleteness\">" << _lcompleteness[i] << "</Completeness>" << std::endl;
            ss << "  <nTot>" << std::fixed << std::setprecision(1) << _lb.binContains(i) << "</nTot>" << std::endl;
            ss << "  <nObs>" << std::fixed << std::setprecision(1) << _lb.binContains(i)*_lcompleteness[i] <<  "</nObs>" << std::endl;
        }
        ss << "</LowResoCompleteness>" << std::endl;
        ss << " <Comment id=\"CompletenessLow\">" << std::endl;
        ss << "Low completeness at low resolution can lead to map distortions and other difficulties.  This often arises through experimental effects such as incorrectly alligned crystal, poorly positioned backstop, or over exposure." << std::endl;
        ss << "</Comment>" << std::endl;
        ss << "</DataStatistics>" << std::endl;
		
		_ira.xml_output(ss);
		
		_wilsonB.xml_output(ss);
		
        _ora.xml_output(ss);
        
		return ss;
	}
	

    //---------Instantiate templates------------------------------------
    
	template class Completeness<clipper::datatypes::F_sigF<clipper::ftype32> >;
	template class Completeness<clipper::datatypes::I_sigI<clipper::ftype32> >;
	
	//template class tNCS<clipper::ftype32>;
	//template class tNCS<clipper::ftype64>;
		
	template class ResoCorrel<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::datatypes::I_sigI<clipper::ftype32> >;
	template class ResoCorrel<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::datatypes::F_sigF<clipper::ftype32> >;
	
	//template bool IceRings_analyse::operator()<clipper::data32::I_sigI>(clipper::HKL_data<clipper::data32::I_sigI>&,ctruncate::Rings&);
}
