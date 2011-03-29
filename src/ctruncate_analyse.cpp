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

namespace ctruncate {
	
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
		std::vector<float> weights(_nbins,1.0f);
		std::vector<float> x(_nbins);
		std::vector<float> y(_nbins);
		
		for( int i = 0 ; i != _nbins ; ++i) {
			float dist = (float(i)+0.5f)*_width/float(_nbins);
			x[i] = 0.25*dist*dist;
			y[i] = std::log(1.0f/_patterson[i]);
		}
		
		float a, b, siga, sigb;
		
		straight_line_fit(x,y,weights,_nbins,a,b,siga,sigb);
		
		_P000 = std::exp(-a);
		
		float btemp = 0.25*b;
		
		_sigma = std::sqrt(1.0f/std::abs(2.0f*btemp) );
		
		return;	
	}
	
	
}