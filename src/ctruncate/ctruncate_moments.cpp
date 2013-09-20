//
//     CTRUNCATE
//     Copyright (C) 2006-2011 Norman Stein, Charles Ballard
//
//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Application.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//

#include "ctruncate_moments.h"
#include "intensity_target.h"
#include <cstdio>

namespace ctruncate {
		

	//-------Moments calculation -------------------------------------------------------
	
	template <class D> Moments<D>::Moments(const clipper::HKL_data<D>& data) {
		hkl_data = &data;
		calc_range = data.hkl_info().invresolsq_range();
		std::string types(this->type());
		mode = (!types.compare("F_sigF") ) ? true : false;
		pn = ( mode ) ? 2.0 : 1.0;
		for (int i=0; i !=5; ++i ) acache[i]=ccache[i]=0.0;
		norm(data);
	}
	
	template <class D> Moments<D>::Moments(const clipper::HKL_data<D>& data, const clipper::Range<clipper::ftype>& range) {
		hkl_data = &data;
		calc_range = range;
		std::string types(this->type());
		mode = (!types.compare("F_sigF") ) ? true : false;
		pn = ( mode ) ? 2.0 : 1.0;
		for (int i=0; i !=5; ++i ) acache[i]=ccache[i]=0.0;
		norm(data);
	}
	
	template <class D> Moments<D>::~Moments() {
	}
	
	template <class D> void Moments<D>::norm(const clipper::HKL_data<D>& data) {
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		na=nc=0;
		clipper::HKL_data<D> fa(data.hkl_info() ), fc(data.hkl_info() );
		for ( HRI ih = data.first(); !ih.last(); ih.next() ) {
			if ( !data[ih].missing() ) {
				if ( !ih.hkl_class().centric() ) {
					na += ih.hkl_class().epsilonc();
					fa[ih] = D(obs(data[ih]),sigobs(data[ih]) );
				} else {
					nc += ih.hkl_class().epsilonc();
					fc[ih] = D(obs(data[ih]),sigobs(data[ih]) );
				}
			}
		}
		
		abins = std::min(1000,na/100+1);
		cbins = std::min(1000,nc/100+1);
		
			std::vector<clipper::ftype> params_init( abins, 1.0 );
			clipper::BasisFn_binner acentric_norm(fa, abins, 1.0);
			
			TargetFn_meanInth<D >  atarget_fn( fa, pn );
		clipper::ResolutionFn af( data.hkl_info() , acentric_norm, atarget_fn, params_init );
			
		clipper::BasisFn_binner centric_norm(fc, cbins, 1.0);
			
		TargetFn_meanInth<D >  ctarget_fn( fc, pn );
		clipper::ResolutionFn cf( data.hkl_info() , centric_norm, ctarget_fn, params_init );
		
		
		normalised.init(data );
		for ( HRI ih = data.first(); !ih.last(); ih.next() ) {
			if ( !data[ih].missing() ) {
				if ( !ih.hkl_class().centric() ) {
					//na += ih.hkl_class().epsilonc();
					normalised[ih] = D(obs(data[ih])/std::pow(acentric_norm.f_s(ih.invresolsq(), af.params() ),1.0/pn),sigobs(data[ih]) );
				} else {
					//nc += ih.hkl_class().epsilonc();
					normalised[ih] = D(obs(data[ih])/std::pow(centric_norm.f_s(ih.invresolsq(), cf.params() ),1.0/pn),sigobs(data[ih]) );
				}
			}
		}
	}
	
	template <class D> clipper::ftype Moments<D>::calc(clipper::ftype power, bool centric) {
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		clipper::HKL_data<D> const& data = *hkl_data;
		
		clipper::ftype val(0.0);
		int nval(0);
		for ( HRI ih = data.first(); !ih.last(); ih.next() ) {
			if (calc_range.contains(ih.invresolsq() ) ) {
				if ( !data[ih].missing() ) {
					if ( ih.hkl_class().centric() == centric ) {
						val += ih.hkl_class().epsilonc()*std::pow(double(obs(normalised[ih]) ),double(power) );
						nval += ih.hkl_class().epsilonc();
					}
				}
			}
		}
		return val/clipper::ftype(nval);
	}
	
	template <class D> clipper::ftype Moments<D>::calc_var(bool centric) {
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		const clipper::HKL_data<D>& data = *hkl_data;
		
		clipper::ftype val(0.0);
		int nval(0);
		clipper::ftype mean = calc(2.0,centric);
		for ( HRI ih = data.first(); !ih.last(); ih.next() ) {
			if (calc_range.contains(ih.invresolsq() ) ) {
				if ( !data[ih].missing() ) {
					if ( ih.hkl_class().centric() == centric) {
						val += ih.hkl_class().epsilonc()*std::fabs(std::pow(double(obs(normalised[ih]) ),2.0) - mean);
						nval += ih.hkl_class().epsilonc();
					}
				}
			}
		}
		return val/clipper::ftype(nval);
	}
	
	template <class D> clipper::ftype Moments<D>::acentric_first() {
		if ( acache[0] != 0.0 ) return acache[0];
		return calc(1.0, false);
	}
	
	template <class D> clipper::ftype Moments<D>::acentric_second() {
		if ( acache[1] != 0.0 ) return acache[1];
		return calc(2.0, false);
	}
	
	template <class D> clipper::ftype Moments<D>::acentric_third() {
		if ( acache[2] != 0.0 ) return acache[2];
		return calc(3.0, false);
	}
	
	template <class D> clipper::ftype Moments<D>::acentric_fourth() {
		if ( acache[3] != 0.0 ) return acache[3];
		return calc(4.0, false);
	}
	
	template <class D> clipper::ftype Moments<D>::acentric_variance() {
		if ( acache[4] != 0.0 ) return acache[4];
		return calc_var(false);
	}
	
	template <class D> clipper::ftype Moments<D>::centric_first() {
		if ( ccache[0] != 0.0 ) return ccache[0];
		return calc(1.0, true);
	}
	
	template <class D> clipper::ftype Moments<D>::centric_second() {
		if ( ccache[1] != 0.0 ) return ccache[1];
		return calc(2.0, true);
	}
	
	template <class D> clipper::ftype Moments<D>::centric_third() {
		if ( ccache[2] != 0.0 ) return ccache[2];
		return calc(3.0, true);
	}
	
	template <class D> clipper::ftype Moments<D>::centric_fourth() {
		if ( ccache[3] != 0.0 ) return ccache[3];
		return calc(4.0, true);
	}
	
	template <class D> clipper::ftype Moments<D>::centric_variance() {
		if ( ccache[4] != 0.0 ) return ccache[4];
		return calc_var(true);
	}
		
	template <class D> clipper::ftype Moments<D>::theo_untwinned_centric_first() {
		if ( mode ) return E_THEO_UNTWINNED_CENTRIC_FIRST;
		else return Z_THEO_UNTWINNED_CENTRIC_FIRST;
	}
	
	template <class D> clipper::ftype Moments<D>::theo_untwinned_centric_second() {
		if ( mode ) return E_THEO_UNTWINNED_CENTRIC_SECOND;
		else return E_THEO_UNTWINNED_CENTRIC_FOURTH;
	}
	
	template <class D> clipper::ftype Moments<D>::theo_untwinned_centric_third() {
		if ( mode ) return E_THEO_UNTWINNED_CENTRIC_THIRD;
		else return Z_THEO_UNTWINNED_CENTRIC_THIRD;
	}
	
	template <class D> clipper::ftype Moments<D>::theo_untwinned_centric_fourth() {
		if ( mode ) return E_THEO_UNTWINNED_CENTRIC_FOURTH;
		else return Z_THEO_UNTWINNED_CENTRIC_FOURTH;
	}
		
	template <class D> clipper::ftype Moments<D>::theo_untwinned_centric_variance() {
		if ( mode ) return E_THEO_UNTWINNED_CENTRIC_VAR;
		else return Z_THEO_UNTWINNED_ACENTRIC_VAR;
	}
	
	template <class D> clipper::ftype Moments<D>::theo_perfect_centric_first() {
		if ( mode ) return E_THEO_PERFECT_CENTRIC_FIRST;
		else return Z_THEO_PERFECT_CENTRIC_FIRST;
	}
	
	template <class D> clipper::ftype Moments<D>::theo_perfect_centric_second() {
		if ( mode ) return E_THEO_PERFECT_CENTRIC_SECOND;
		else return E_THEO_PERFECT_CENTRIC_FOURTH;
	}
	
	template <class D> clipper::ftype Moments<D>::theo_perfect_centric_third() {
		if ( mode ) return E_THEO_PERFECT_CENTRIC_THIRD;
		else return Z_THEO_PERFECT_CENTRIC_THIRD;
	}
	
	template <class D> clipper::ftype Moments<D>::theo_perfect_centric_fourth() {
		if ( mode ) return E_THEO_PERFECT_CENTRIC_FOURTH;
		else return Z_THEO_PERFECT_CENTRIC_FOURTH;
	}
	
	template <class D> clipper::ftype Moments<D>::theo_untwinned_acentric_first() {
		return theo_perfect_centric_first();
	}
	
	template <class D> clipper::ftype Moments<D>::theo_untwinned_acentric_second() {
		return theo_perfect_centric_second();
	}
	
	template <class D> clipper::ftype Moments<D>::theo_untwinned_acentric_third() {
		return theo_perfect_centric_third();
	}
	
	template <class D> clipper::ftype Moments<D>::theo_untwinned_acentric_fourth() {
		return theo_perfect_centric_fourth();
	}
	
	template <class D> clipper::ftype Moments<D>::theo_untwinned_acentric_variance() {
		if ( mode ) return E_THEO_UNTWINNED_ACENTRIC_VAR;
		else return Z_THEO_UNTWINNED_ACENTRIC_VAR;
	}
	
	template <class D> clipper::ftype Moments<D>::theo_perfect_acentric_first() {
		if ( mode ) return E_THEO_PERFECT_ACENTRIC_FIRST;
		else return Z_THEO_PERFECT_ACENTRIC_FIRST;
	}
	
	template <class D> clipper::ftype Moments<D>::theo_perfect_acentric_second() {
		if ( mode ) return E_THEO_PERFECT_ACENTRIC_SECOND;
		else return E_THEO_PERFECT_ACENTRIC_FOURTH;
	}
	
	template <class D> clipper::ftype Moments<D>::theo_perfect_acentric_third() {
		if ( mode ) return E_THEO_PERFECT_ACENTRIC_THIRD;
		else return Z_THEO_PERFECT_ACENTRIC_THIRD;
	}
	
	template <class D> clipper::ftype Moments<D>::theo_perfect_acentric_fourth() {
		if ( mode ) return E_THEO_PERFECT_ACENTRIC_FOURTH;
		else return Z_THEO_PERFECT_ACENTRIC_FOURTH;
	}

	//-------Moments estimate twinning fraction -------------------------------------------------------
	
	/*! estimate twinning fraction from <I^2>/<I>^2 estimation
	 \param (void)
	 \return fraction (clipper::ftype) */
	template <class D> clipper::ftype Moments<D>::fraction() {
		clipper::ftype val = (mode) ? acentric_fourth() : acentric_second();
		if ( val > E_THEO_PERFECT_CENTRIC_FOURTH ) {
			return 0.0;
		} else if ( val < E_THEO_PERFECT_ACENTRIC_FOURTH) {
			return 0.5;
		} else {
			return  0.5 - sqrt(0.5*val - 0.75);
		}
	}
	
	//-------Moments summary-------------------------------------------------------
	
	/*! output summary for data
	 \param void
	 \return void */
	template <class D> void Moments<D>::summary() {
		if (mode ) {
			printf("\n\nMEAN ACENTRIC MOMENTS OF E:\n\n");
			printf("<E> = %6.3f (Expected value = %6.3f, Perfect Twin = %6.3f)\n", acentric_first(), theo_untwinned_acentric_first(), theo_perfect_acentric_first() );
			printf("<|E^2-1|> = %6.3f (Expected = %6.3f)\n", acentric_variance(), theo_untwinned_acentric_variance() );
			printf("<F^3>/<F>^3 = %6.3f (Expected value = %6.3f, Perfect Twin = %6.3f)\n", acentric_third(), theo_untwinned_acentric_third(), theo_perfect_acentric_third() );
			printf("<I^2>/<I>^2 = %6.3f (Expected value = %6.3f, Perfect Twin = %6.3f)\n", acentric_fourth(), theo_untwinned_acentric_fourth(), theo_perfect_acentric_fourth());
			printf("(equivalent to twin fraction of %6.3f)\n",this->fraction() );
		} else {
			printf("\n\nMEAN ACENTRIC MOMENTS OF Z:\n\n");
			printf("<I^2>/<I>^2 = %6.3f (Expected = %6.3f, Perfect Twin = %6.3f)\n", acentric_second(), theo_untwinned_acentric_second(), theo_perfect_acentric_second() );
			printf("<I^3>/<I>^3 = %6.3f (Expected value = %6.3f, Perfect Twin = %6.3f)\n", acentric_third(), theo_untwinned_acentric_third(), theo_perfect_acentric_third() );
			printf("<I^4>/<I>^4 = %6.3f (Expected value = %6.3f, Perfect Twin = %6.3f)\n", acentric_fourth(), theo_untwinned_acentric_fourth(), theo_perfect_acentric_fourth());
			printf("(equivalent to twin fraction of %6.3f)\n",this->fraction() );
		}
		printf("\n\n");
	}
	
	
		//-------Moments loggraph-------------------------------------------------------
		
		/*! output loggraph for data
		 \param void
		 \return type void */
		template <class D> void Moments<D>::loggraph() {
			typedef clipper::HKL_data_base::HKL_reference_index HRI;
				
			int ncentric(0);
			int nacentric(0);
			
			for ( HRI ih = normalised.first(); !ih.last(); ih.next() ) {
				if ( !normalised[ih].missing() )
					if (ih.hkl_class().centric()) ncentric += 1;
			}
			
			int nbins(30);
			int ncbins = std::min(nbins, ncentric/10);
				
			std::vector<int> Na(nbins,0),  Nc(ncbins,0);
			std::vector<double> I1a(nbins,0.0), I1c(ncbins,0.0);
			std::vector<double> I2a(nbins,0.0), I2c(ncbins,0.0);
			std::vector<double> I3a(nbins,0.0), I3c(ncbins,0.0);
			std::vector<double> I4a(nbins,0.0), I4c(ncbins,0.0);
			std::vector<double> I5a(nbins,0.0), I5c(ncbins,0.0);
				
			clipper::Range<clipper::ftype> range=normalised.hkl_info().invresolsq_range();
			clipper::ftype maxres = range.max();
			
			clipper::ftype meana = calc(2.0,false);
			clipper::ftype meanc = calc(2.0,true);
			
			for ( HRI ih = normalised.first(); !ih.last(); ih.next() ) {
				if ( !normalised[ih].missing() ) {
					int n =ih.hkl_class().epsilonc();
					double val = obs(normalised[ih]);
					int bin = int( nbins * ( ih.invresolsq() / double(maxres) ) - 0.5  );
					int cbin = int( ncbins * ( ih.invresolsq() / double(maxres) ) - 0.5  );
						if (!ih.hkl_class().centric()) {
							Na[bin]+=n;
							if (val > 0.0) {
								I1a[bin] += n*val;
								I2a[bin] += n*val*val;
								I3a[bin] += n*val*val*val;
								I4a[bin] += n*val*val*val*val;
								I5a[bin] += n*std::fabs(val*val-meana);
							}
						}
						else if (ncentric != 0) {
							Nc[cbin]+=n;
							if (val > 0.0) {
								I1c[cbin] += n*val;
								I2c[cbin] += n*val*val;
								I3c[cbin] += n*val*val*val;
								I4c[cbin] += n*val*val*val*val;
								I5c[cbin] += n*std::fabs(val*val-meanc);
							}
						}
					}
				}
			
			for ( int bin=0 ; bin != nbins ; bin++ ) {
				if ( Na[bin] != 0 ) {
				I1a[bin] /= Na[bin];
				I2a[bin] /= Na[bin];
				I3a[bin] /= Na[bin];
				I4a[bin] /= Na[bin];
				I5a[bin] /= Na[bin];	
				}
			}
			
			for ( int bin=0 ; bin != ncbins ; bin++ ) {
				if ( Nc[bin] != 0 ) {
				I1c[bin] /= Nc[bin];
				I2c[bin] /= Nc[bin];
				I3c[bin] /= Nc[bin];
				I4c[bin] /= Nc[bin];
				I5c[bin] /= Nc[bin];				
				} }
				
			if (mode) {
				printf("$TABLE: Acentric moments of E:\n");
				printf("$GRAPHS");
				printf(": 4th moment of E %5.3f (Expected value = 2, Perfect Twin = 1.5):0|%5.3fx0|5:1,5:\n", acentric_fourth(), maxres);
				printf(": 1st & 3rd moments of E (Expected values = 0.886, 1.329, Perfect twin = 0.94, 1.175):0|%5.3fx0|2:1,2,4:\n", maxres);
				printf("$$ 1/resol^2 <E> <|E**2-1|> <E**3> <E**4> $$\n$$\n");
				
				for (int i=0; i<nbins; i++) {
					double res = maxres * (double(i) + 0.5)/double(nbins);
					double n = double(Na[i]);
					printf("%10.6f %10.6f %10.6f %10.6f %10.6f\n", res, I1a[i], I5a[i], I3a[i], I4a[i]);
				}
				printf("$$\n\n");
				
				if (ncentric != 0) {
					printf("$TABLE: Centric moments of E:\n");
					printf("$GRAPHS");
					printf(": 4th moment of E %5.3f (Expected = 3, Perfect Twin = 2):0|%5.3fx0|5:1,5:\n", centric_fourth(), maxres);
					printf(": 1st & 3rd moments of E (Expected = 0.798, 1.596, Perfect Twin = 0.886, 1.329):0|%5.3fx0|4:1,2,4:\n", maxres);
					printf("$$ 1/resol^2 <E> <|E**2-1|> <E**3> <E**4> $$\n$$\n");
					
					for (int i=0; i<ncbins; i++) {
						double res = maxres * (double(i) + 0.5)/double(ncbins);
						double n = double(Nc[i]);
						printf("%10.6f %10.6f %10.6f %10.6f %10.6f\n", res, I1c[i], I5c[i], I3c[i], I4c[i]);
					}
					printf("$$\n\n");
				}
			} else {
				printf("$TABLE: Acentric moments of I:\n");
				printf("$GRAPHS");
				printf(": 2nd moment of I %5.3f (Expected value = 2, Perfect Twin = 1.5):0|%5.3fx0|5:1,2:\n", acentric_second(), maxres);
				printf(": 3rd & 4th moments of I (Expected values = 6, 24, Perfect twin = 3, 7.5):0|%5.3fx0|2:1,3,4:\n", maxres);
				printf("$$ 1/resol^2 <I**2> <I**3> <I**4> $$\n$$\n");
				
				for (int i=0; i<nbins; i++) {
					double res = maxres * (double(i) + 0.5)/double(nbins);
					double n = double(Na[i]);
					printf("%10.6f %10.6f %10.6f %10.6f\n", res, I2a[i], I3a[i], I4a[i]);
				}
				printf("$$\n\n");
				
				if (ncentric != 0) {
					printf("$TABLE: Centric moments of I:\n");
					printf("$GRAPHS");
					printf(": 2nd moment of I %5.3f (Expected = 3, Perfect Twin = 2):0|%5.3fx0|5:1,2:\n", centric_second(), maxres);
					printf(": 3rd & 4th moments of I (Expected = 15, 105, Perfect Twin = 6, 24):0|%5.3fx0|4:1,3,4:\n", maxres);
					printf("$$ 1/resol^2 <I**2> <I**3> <I**4> $$\n$$\n");
					
					for (int i=0; i<ncbins; i++) {
						double res = maxres * (double(i) + 0.5)/double(ncbins);
						double n = double(Nc[i]);
						printf("%10.6f %10.6f %10.6f %10.6f\n", res, I2c[i], I3c[i], I4c[i]);
					}
					printf("$$\n\n");
				}
			}
			return;
		}
	
	template class Moments<clipper::datatypes::F_sigF<float> >;
	template class Moments<clipper::datatypes::F_sigF<double> >;
	template class Moments<clipper::datatypes::I_sigI<float> >;
	template class Moments<clipper::datatypes::I_sigI<double> >;
	
}

/*

statistic	Untwinned	Perfect		Untwinned	Perfect		Partial		Partial
			centric		centric		acentric	acentric	centric		acentric
<|H|>														2(1-2a)/pi	(1-2a)/2
<|H^2|>														(1-2a)^2/2	(1-2a)^2/3
<|L|>		2/pi		1/2			1/2			3/8
<|L^2|>		1/2			1/3			1/3			1/5
<|E^2-1|>	0.968		too low		0.736		too low
<E>			0.798		0.886		0.886		0.940
<E^3>		1.596		1.329		1.329		1.175
<E^4>,<Z^2>	3.0			2.0			2.0			1.5
<Z^3>		15.0		6.0			6.0			3.0
<Z^4>		105.0		24.0		24.0		7.5
*/
