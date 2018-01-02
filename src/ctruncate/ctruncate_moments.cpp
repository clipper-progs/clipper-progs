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
#include <iomanip>

namespace ctruncate {
    
    
	//-------Moments calculation -------------------------------------------------------
	
	
	/*template <class T, template<class> class D> Moments::Moments(clipper::HKL_data<D<T> >& data) {
     hkl_data = &data;
     calc_range = data.hkl_info().invresolsq_range();
     std::string types(this->type());
     mode = (!types.compare("F_sigF") ) ? true : false;
     pn = ( mode ) ? 2.0 : 1.0;
     for (int i=0; i !=5; ++i ) acache[i]=ccache[i]=0.0;
     norm(data);
     }*/
	
	/*template <class T, template<class> class D> Moments::Moments(clipper::HKL_data<D<T> >& data, clipper::Range<clipper::ftype>& range) {
     hkl_data = &data;
     calc_range = range;
     std::string types(this->type());
     mode = (!types.compare("F_sigF") ) ? true : false;
     pn = ( mode ) ? 2.0 : 1.0;
     for (int i=0; i !=5; ++i ) acache[i]=ccache[i]=0.0;
     norm(data);
     }*/
	
	Moments::~Moments() {
	}
	
	
	/*template <class T, template<class> class D> void Moments::operator()(const clipper::HKL_data<D<T> >& data, const clipper::Range<clipper::ftype>& range) {
     hkl_data = &data;
     calc_range = range;
     std::string types(this->type());
     mode = (!types.compare("F_sigF") ) ? true : false;
     pn = ( mode ) ? 2.0 : 1.0;
     for (int i=0; i !=5; ++i ) acache[i]=ccache[i]=0.0;
     norm(data);
     }*/
	
	/*template <class T, template<class> class D> void Moments::norm(const clipper::HKL_data<D<T> >& data) {
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
     
     TargetFn_meanInth<D<T> >  atarget_fn( fa, pn );
     clipper::ResolutionFn af( data.hkl_info() , acentric_norm, atarget_fn, params_init );
     
     clipper::BasisFn_binner centric_norm(fc, cbins, 1.0);
     
     TargetFn_meanInth<D<T> >  ctarget_fn( fc, pn );
     clipper::ResolutionFn cf( data.hkl_info() , centric_norm, ctarget_fn, params_init );
     
     
     normalised.init(data );
     for ( HRI ih = data.first(); !ih.last(); ih.next() ) {
     if ( !data[ih].missing() ) {
     if ( !ih.hkl_class().centric() ) {
     //na += ih.hkl_class().epsilonc();
     normalised[ih] = (obs(data[ih])/std::pow(acentric_norm.f_s(ih.invresolsq(), af.params() ),1.0/pn),sigobs(data[ih]) );
     } else {
     //nc += ih.hkl_class().epsilonc();
     normalised[ih] = (obs(data[ih])/std::pow(centric_norm.f_s(ih.invresolsq(), cf.params() ),1.0/pn),sigobs(data[ih]) );
     }
     }
     }
     } */
	
	clipper::ftype Moments::calc(clipper::ftype power, bool centric) {
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		
		clipper::ftype val(0.0),nval(0.0);
		bool is(is_intensity() );
		for ( HRI ih = _base->first(); !ih.last(); ih.next() ) {
			if (calc_range.contains(ih.invresolsq() ) ) {
				if ( !_base->missing(ih.index() ) ) {
					if ( ih.hkl_class().centric() == centric ) {
						clipper::ftype eps = ( is ) ? 1.0/ih.hkl_class().epsilon() : 1.0/std::sqrt(ih.hkl_class().epsilon());
						val += std::pow(double(eps*obs(_normalised[ih]) ),double(power) );
						nval += 1.0;
					}
				}
			}
		}
		return val/clipper::ftype(nval);
	}
	
	clipper::ftype Moments::calc_var(bool centric) {
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		clipper::ftype val(0.0),nval(0.0);
		clipper::ftype mean = calc(2.0,centric);
		bool is(is_intensity() );
		for ( HRI ih = _base->first(); !ih.last(); ih.next() ) {
			if (calc_range.contains(ih.invresolsq() ) ) {
				if ( !_base->missing(ih.index() ) ) {
					if ( ih.hkl_class().centric() == centric) {
						clipper::ftype eps = ( is ) ? 1.0/ih.hkl_class().epsilon() : 1.0/std::sqrt(ih.hkl_class().epsilon());
						val += std::pow(eps*obs(_normalised[ih]),2);
						nval += 1.0;
					}
				}
			}
		}
		return val/clipper::ftype(nval)-mean;
	}
	
	clipper::ftype Moments::acentric_first() {
		if ( acache[0] != 0.0 ) return acache[0];
		return calc(1.0, false);
	}
	
	clipper::ftype Moments::acentric_second() {
		if ( acache[1] != 0.0 ) return acache[1];
		return calc(2.0, false);
	}
	
	clipper::ftype Moments::acentric_third() {
		if ( acache[2] != 0.0 ) return acache[2];
		return calc(3.0, false);
	}
	
	clipper::ftype Moments::acentric_fourth() {
		if ( acache[3] != 0.0 ) return acache[3];
		return calc(4.0, false);
	}
	
	clipper::ftype Moments::acentric_variance() {
		if ( acache[4] != 0.0 ) return acache[4];
		return calc_var(false);
	}
	
	clipper::ftype Moments::centric_first() {
		if ( ccache[0] != 0.0 ) return ccache[0];
		return calc(1.0, true);
	}
	
	clipper::ftype Moments::centric_second() {
		if ( ccache[1] != 0.0 ) return ccache[1];
		return calc(2.0, true);
	}
	
	clipper::ftype Moments::centric_third() {
		if ( ccache[2] != 0.0 ) return ccache[2];
		return calc(3.0, true);
	}
	
	clipper::ftype Moments::centric_fourth() {
		if ( ccache[3] != 0.0 ) return ccache[3];
		return calc(4.0, true);
	}
	
	clipper::ftype Moments::centric_variance() {
		if ( ccache[4] != 0.0 ) return ccache[4];
		return calc_var(true);
	}
	
	clipper::ftype Moments::theo_untwinned_centric_first() {
		if ( is_intensity() ) return Z_THEO_UNTWINNED_CENTRIC_FIRST;
		else return E_THEO_UNTWINNED_CENTRIC_FIRST;
	}
	
	clipper::ftype Moments::theo_untwinned_centric_second() {
		if ( is_intensity() ) return E_THEO_UNTWINNED_CENTRIC_FOURTH;
		else return E_THEO_UNTWINNED_CENTRIC_SECOND;
	}
	
	clipper::ftype Moments::theo_untwinned_centric_third() {
		if ( is_intensity() ) return Z_THEO_UNTWINNED_CENTRIC_THIRD;
		else return E_THEO_UNTWINNED_CENTRIC_THIRD;
	}
	
	clipper::ftype Moments::theo_untwinned_centric_fourth() {
		if ( is_intensity() ) return Z_THEO_UNTWINNED_CENTRIC_FOURTH;
		else return E_THEO_UNTWINNED_CENTRIC_FOURTH;
	}
	
	clipper::ftype Moments::theo_untwinned_centric_variance() {
		if ( is_intensity() ) return 0.0;
		else return E_THEO_UNTWINNED_CENTRIC_VAR;
	}
	
	clipper::ftype Moments::theo_perfect_centric_first() {
		if ( is_intensity() ) return Z_THEO_PERFECT_CENTRIC_FIRST;
		else return E_THEO_PERFECT_CENTRIC_FIRST;
	}
	
	clipper::ftype Moments::theo_perfect_centric_second() {
		if ( is_intensity() ) return E_THEO_PERFECT_CENTRIC_FOURTH;
		else return E_THEO_PERFECT_CENTRIC_SECOND;
	}
	
	clipper::ftype Moments::theo_perfect_centric_third() {
		if ( is_intensity() ) return Z_THEO_PERFECT_CENTRIC_THIRD;
		else return E_THEO_PERFECT_CENTRIC_THIRD;
	}
	
	clipper::ftype Moments::theo_perfect_centric_fourth() {
		if ( is_intensity() ) return Z_THEO_PERFECT_CENTRIC_FOURTH;
		else return E_THEO_PERFECT_CENTRIC_FOURTH;
	}
	
    clipper::ftype Moments::theo_perfect_centric_variance() {
		if ( is_intensity() ) return 0.0;
		else return E_THEO_UNTWINNED_ACENTRIC_VAR;
	}
    
	clipper::ftype Moments::theo_untwinned_acentric_first() {
		return theo_perfect_centric_first();
	}
	
	clipper::ftype Moments::theo_untwinned_acentric_second() {
		return theo_perfect_centric_second();
	}
	
	clipper::ftype Moments::theo_untwinned_acentric_third() {
		return theo_perfect_centric_third();
	}
	
	clipper::ftype Moments::theo_untwinned_acentric_fourth() {
		return theo_perfect_centric_fourth();
	}
	
	clipper::ftype Moments::theo_untwinned_acentric_variance() {
		if ( is_intensity() ) return 0.0;
		else return E_THEO_UNTWINNED_ACENTRIC_VAR;
	}
    
    clipper::ftype Moments::theo_perfect_acentric_variance() {
		if ( is_intensity() ) return 0.0;
		else return E_THEO_PERFECT_ACENTRIC_VAR;
	}
	
	clipper::ftype Moments::theo_perfect_acentric_first() {
		if ( is_intensity() ) return Z_THEO_PERFECT_ACENTRIC_FIRST;
		else return E_THEO_PERFECT_ACENTRIC_FIRST;
	}
	
	clipper::ftype Moments::theo_perfect_acentric_second() {
		if ( is_intensity() ) return E_THEO_PERFECT_ACENTRIC_FOURTH;
		else return E_THEO_PERFECT_ACENTRIC_SECOND;
	}
	
	clipper::ftype Moments::theo_perfect_acentric_third() {
		if ( is_intensity() ) return Z_THEO_PERFECT_ACENTRIC_THIRD;
		else return E_THEO_PERFECT_ACENTRIC_THIRD;
	}
	
	clipper::ftype Moments::theo_perfect_acentric_fourth() {
		if ( is_intensity() ) return Z_THEO_PERFECT_ACENTRIC_FOURTH;
		else return E_THEO_PERFECT_ACENTRIC_FOURTH;
	}
	
	//-------Moments estimate twinning fraction -------------------------------------------------------
	
	/*! estimate twinning fraction from <I^2>/<I>^2 estimation
	 \param (void)
	 \return fraction (clipper::ftype) */
	clipper::ftype Moments::fraction() {
		clipper::ftype val = (is_intensity() ) ? acentric_second() : acentric_fourth();
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
	void Moments::summary() {
		if (is_intensity() ) {
			printf("\n\nMEAN ACENTRIC Moments OF Z:\n\n");
			printf("<I^2>/<I>^2 = %6.3f (Expected = %6.3f, Perfect Twin = %6.3f)\n", acentric_second(), theo_untwinned_acentric_second(), theo_perfect_acentric_second() );
			printf("<I^3>/<I>^3 = %6.3f (Expected value = %6.3f, Perfect Twin = %6.3f)\n", acentric_third(), theo_untwinned_acentric_third(), theo_perfect_acentric_third() );
			printf("<I^4>/<I>^4 = %6.3f (Expected value = %6.3f, Perfect Twin = %6.3f)\n", acentric_fourth(), theo_untwinned_acentric_fourth(), theo_perfect_acentric_fourth());
			printf("(equivalent to twin fraction of %6.3f)\n",this->fraction() );
		} else {
			printf("\n\nMEAN ACENTRIC Moments OF E:\n\n");
			printf("<E> = %6.3f (Expected value = %6.3f, Perfect Twin = %6.3f)\n", acentric_first(), theo_untwinned_acentric_first(), theo_perfect_acentric_first() );
			printf("<|E^2-1|> = %6.3f (Expected = %6.3f)\n", acentric_variance(), theo_untwinned_acentric_variance() );
			printf("<F^3>/<F>^3 = %6.3f (Expected value = %6.3f, Perfect Twin = %6.3f)\n", acentric_third(), theo_untwinned_acentric_third(), theo_perfect_acentric_third() );
			printf("<I^2>/<I>^2 = %6.3f (Expected value = %6.3f, Perfect Twin = %6.3f)\n", acentric_fourth(), theo_untwinned_acentric_fourth(), theo_perfect_acentric_fourth());
			printf("(equivalent to twin fraction of %6.3f)\n",this->fraction() );
		}
		printf("\n\n");
	}
	
	void Moments::output() {
        summary();
        
        std::cout << std::endl << "The acentric moments plot can also give an indication of the quality of the high resolution data.  Large deviations often indicate bias in the data." << std::endl << std::endl;
        
        loggraph();
    }
	//-------Moments loggraph-------------------------------------------------------
	
	/*! output loggraph for data
	 \param void
	 \return type void */
	void Moments::loggraph() {
		if (is_intensity() ) {
			if (_a_reso.size() != 1) {
				printf("$TABLE: Acentric Moments of I:\n");
				printf("$GRAPHS");
				printf(": 2nd moment of I %5.3f (Expected value = 2, Perfect Twin = 1.5):0|%5.3fx0|5:1,2:\n", acentric_second(), _a_reso[_a_reso.size()-1]);
				printf(": 3rd & 4th Moments of I (Expected values = 6, 24, Perfect twin = 3, 7.5):0|%5.3fx0|36:1,3,4:\n", _a_reso[_a_reso.size()-1]);
				printf("$$ 1/resol^2   <I**2>     <I**3>     <I**4> $$\n$$\n");
				
				for (int i=0; i !=_a_reso.size(); ++i) {
					printf("%10.6f %10.3f %10.3f %10.3f\n", _a_reso[i], I2a[i], I3a[i], I4a[i]);
				}
				printf("$$\n\n");
			}
			if (_c_reso.size() != 1) {
				printf("$TABLE: Centric Moments of I:\n");
				printf("$GRAPHS");
				printf(": 2nd moment of I %5.3f (Expected = 3, Perfect Twin = 2):0|%5.3fx0|5:1,2:\n", centric_second(), _c_reso[_c_reso.size()-1]);
				printf(": 3rd & 4th Moments of I (Expected = 15, 105, Perfect Twin = 6, 24):0|%5.3fx0|150:1,3,4:\n", _c_reso[_c_reso.size()-1]);
				printf("$$ 1/resol^2   <I**2>     <I**3>     <I**4> $$\n$$\n");
				
				for (int i=0; i!=_c_reso.size(); i++) {
					printf("%10.6f %10.3f %10.3f %10.3f\n", _c_reso[i], I2c[i], I3c[i], I4c[i]);
				}
				printf("$$\n\n");
			}
		} else {
			if (_a_reso.size() != 1) {
				printf("$TABLE: Acentric Moments of E:\n");
				printf("$GRAPHS");
				printf(": 4th moment of E %5.3f (Expected value = 2, Perfect Twin = 1.5):0|%5.3fx0|5:1,5:\n", acentric_fourth(), _a_reso[_a_reso.size()-1]);
				printf(": 1st & 3rd Moments of E (Expected values = 0.886, 1.329, Perfect twin = 0.94, 1.175):0|%5.3fx0|2:1,2,4:\n", _a_reso[_a_reso.size()-1]);
				printf("$$ 1/resol^2   <E>     <|E**2-1|>   <E**3>     <E**4> $$\n$$\n");
				
				for (int i=0; i !=_a_reso.size(); ++i) {
					printf("%10.6f %810.3f %10.3f %10.3f %10.3f\n", _a_reso[i], I1a[i], I5a[i], I3a[i], I4a[i]);
				}
				printf("$$\n\n");
			}
			if (_c_reso.size() != 1) {
				printf("$TABLE: Centric Moments of E:\n");
				printf("$GRAPHS");
				printf(": 4th moment of E %5.3f (Expected = 3, Perfect Twin = 2):0|%5.3fx0|5:1,5:\n", centric_fourth(), _c_reso[_c_reso.size()-1]);
				printf(": 1st & 3rd Moments of E (Expected = 0.798, 1.596, Perfect Twin = 0.886, 1.329):0|%5.3fx0|4:1,2,4:\n", _c_reso[_c_reso.size()-1]);
				printf("$$ 1/resol^2   <E>     <|E**2-1|>   <E**3>     <E**4> $$\n$$\n");
				
				for (int i=0; i != _c_reso.size(); ++i) {
					printf("%10.6f %10.3f %10.3f %10.3f %10.3f\n", _c_reso[i], I1c[i], I5c[i], I3c[i], I4c[i]);
				}
				printf("$$\n\n");
			}
		}
		return;
	}
	
	std::stringstream& Moments::xml_output(std::stringstream& ss) {
        if (_a_reso.size() != 1) {
            ss << "  <Moments id=\"acentric\" type=\""<< ((is_intensity() ) ? "intensity" : "amplitudes") << "\">" << std::endl;
            if (is_intensity() ) {
                ss << "    <Moment id=\"&lt;I^2&gt;/&lt;I&gt;^2\">" << std::endl;
                ss << "       <value id=\"&lt;I^2&gt;/&lt;I&gt;^2\">" << std::fixed << std::setw(6) << std::setprecision(3) << acentric_second() << "</value>" << std::endl;
                ss << "       <untwinned id=\"&lt;I^2&gt;/&lt;I&gt;^2\">" << std::fixed << std::setw(6) << std::setprecision(3) << theo_untwinned_acentric_second() << " </untwinned>" << std::endl;
                ss << "       <twinned id=\"&lt;I^2&gt;/&lt;I&gt;^2\">" << std::fixed << std::setw(6) << std::setprecision(3) << theo_perfect_acentric_second() << " </twinned>" << std::endl;
                ss << "    </Moment>" << std::endl;
                ss << "    <Moment id=\"&lt;I^3&gt;/&lt;I&gt;^3\">" << std::endl;
                ss << "       <value id=\"&lt;I^3&gt;/&lt;I&gt;^3\">" << std::fixed << std::setw(6) << std::setprecision(3) << acentric_third() << "</value>" << std::endl;
                ss << "       <untwinned id=\"&lt;I^3&gt;/&lt;I&gt;^3\">" << std::fixed << std::setw(6) << std::setprecision(3) << theo_untwinned_acentric_third() << " </untwinned>" << std::endl;
                ss << "       <twinned id=\"&lt;I^3&gt;/&lt;I&gt;^3\">" << std::fixed << std::setw(6) << std::setprecision(3) << theo_perfect_acentric_third() << " </twinned>" << std::endl;
                ss << "    </Moment>" << std::endl;
                ss << "    <Moment id=\"&lt;I^4&gt;/&lt;I&gt;^4\">" << std::endl;
                ss << "       <value id=\"&lt;I^4&gt;/&lt;I&gt;^4\">" << std::fixed << std::setw(6) << std::setprecision(3) << acentric_fourth() << "</value>" << std::endl;
                ss << "       <untwinned id=\"&lt;I^4&gt;/&lt;I&gt;^4\">" << std::fixed << std::setw(6) << std::setprecision(3) << theo_untwinned_acentric_fourth() << " </untwinned>" << std::endl;
                ss << "       <twinned id=\"&lt;I^4&gt;/&lt;I&gt;^4\">" << std::fixed << std::setw(6) << std::setprecision(3) << theo_perfect_acentric_fourth() << " </twinned>" << std::endl;
                ss << "    </Moment>" << std::endl;
            } else {
                ss << "    <Moment id=\"&lt;E&gt;\">" << std::endl;
                ss << "       <value id=\"&lt;E&gt;\">" << std::fixed << std::setw(6) << std::setprecision(3) << acentric_first() << "</value>" << std::endl;
                ss << "       <untwinned id=\"&lt;E&gt;\">" << std::fixed << std::setw(6) << std::setprecision(3) << theo_untwinned_acentric_first() << " </untwinned>" << std::endl;
                ss << "       <twinned id=\"&lt;E&gt;\">" << std::fixed << std::setw(6) << std::setprecision(3) << theo_perfect_acentric_first() << " </twinned>" << std::endl;
                ss << "    </Moment>" << std::endl;
                ss << "    <Moment id=\"&lt;|E^2-1|&gt;\">" << std::fixed << std::setw(6) << std::setprecision(3) << std::endl;
                ss << "       <value id=\"&lt;|E^2-1|&gt;\">" << std::fixed << std::setw(6) << std::setprecision(3) << acentric_variance() << "</value>" << std::endl;
                ss << "       <untwinned id=\"&lt;|E^2-1|&gt;\">" << std::fixed << std::setw(6) << std::setprecision(3) << theo_untwinned_acentric_variance() << " </untwinned>" << std::endl;
                ss << "       <twinned id=\"&lt;|E^2-1|&gt;\">" << std::fixed << std::setw(6) << std::setprecision(3) << theo_perfect_acentric_variance() << " </twinned>" << std::endl;
                ss << "    </Moment>" << std::endl;
                
                ss << "    <Moment id=\"&lt;F^3&gt;/&lt;F&gt;^3\">" << std::endl;
                ss << "       <value id=\"&lt;F^3&gt;/&lt;F&gt;^3\">" << std::fixed << std::setw(6) << std::setprecision(3) << acentric_third() << "</value>" << std::endl;
                ss << "       <untwinned id=\"&lt;F^3&gt;/&lt;F&gt;^3\">" << std::fixed << std::setw(6) << std::setprecision(3) << theo_untwinned_acentric_third() << " </untwinned>" << std::endl;
                ss << "       <twinned id=\"&lt;F^3&gt;/&lt;F&gt;^3\">" << std::fixed << std::setw(6) << std::setprecision(3) << theo_perfect_acentric_third() << " </twinned>" << std::endl;
                ss << "    </Moment>" << std::endl;
                ss << "    <Moment id=\"&lt;Z^2&gt;\">" << std::endl;
                ss << "       <value id=\"&lt;Z^2&gt;\">" << std::fixed << std::setw(6) << std::setprecision(3) << acentric_fourth() << "</value>" << std::endl;
                ss << "       <untwinned id=\"&lt;Z^2&gt;\">" << std::fixed << std::setw(6) << std::setprecision(3) << theo_untwinned_acentric_fourth() << " </untwinned>" << std::endl;
                ss << "       <twinned id=\"&lt;Z^2&gt;\">" << std::fixed << std::setw(6) << std::setprecision(3) << theo_perfect_acentric_fourth() << " </twinned>" << std::endl;
                ss << "    </Moment>" << std::endl;
                
            }
            ss << "  </Moments>" << std::endl;
        }
        if (_c_reso.size() != 1) {
            ss << "  <Moments id=\"centric\" type=\""<< ((is_intensity() ) ? "intensity" : "amplitudes") << "\">" << std::endl;
            if (is_intensity() ) {
                ss << "    <Moment id=\"&lt;I^2&gt;/&lt;I&gt;^2\">" << std::endl;
                ss << "       <value id=\"&lt;I^2&gt;/&lt;I&gt;^2\">" << std::fixed << std::setw(6) << std::setprecision(3) << centric_second() << "</value>" << std::endl;
                ss << "       <untwinned id=\"&lt;I^2&gt;/&lt;I&gt;^2\">" << std::fixed << std::setw(6) << std::setprecision(3) << theo_untwinned_centric_second() << " </untwinned>" << std::endl;
                ss << "       <twinned id=\"&lt;I^2&gt;/&lt;I&gt;^2\">" << std::fixed << std::setw(6) << std::setprecision(3) << theo_perfect_centric_second() << " </twinned>" << std::endl;
                ss << "    </Moment>" << std::endl;
                ss << "    <Moment id=\"&lt;I^3&gt;/&lt;I&gt;^3\">" << std::endl;
                ss << "       <value id=\"&lt;I^3&gt;/&lt;I&gt;^3\">" << std::fixed << std::setw(6) << std::setprecision(3) << centric_third() << "</value>" << std::endl;
                ss << "       <untwinned id=\"&lt;I^3&gt;/&lt;I&gt;^3\">" << std::fixed << std::setw(6) << std::setprecision(3) << theo_untwinned_centric_third() << " </untwinned>" << std::endl;
                ss << "       <twinned id=\"&lt;I^3&gt;/&lt;I&gt;^3\">" << std::fixed << std::setw(6) << std::setprecision(3) << theo_perfect_centric_third() << " </twinned>" << std::endl;
                ss << "    </Moment>" << std::endl;
                ss << "    <Moment id=\"&lt;I^4&gt;/&lt;I&gt;^4\">" << std::endl;
                ss << "       <value id=\"&lt;I^4&gt;/&lt;I&gt;^4\">" << std::fixed << std::setw(6) << std::setprecision(3) << centric_fourth() << "</value>" << std::endl;
                ss << "       <untwinned id=\"&lt;I^4&gt;/&lt;I&gt;^4\">" << std::fixed << std::setw(6) << std::setprecision(3) << theo_untwinned_centric_fourth() << " </untwinned>" << std::endl;
                ss << "       <twinned id=\"&lt;I^4&gt;/&lt;I&gt;^4\">" << std::fixed << std::setw(6) << std::setprecision(3) << theo_perfect_centric_fourth() << " </twinned>" << std::endl;
                ss << "    </Moment>" << std::endl;
            } else {
                ss << "    <Moment id=\"&lt;E&gt;\">" << std::endl;
                ss << "       <value id=\"&lt;E&gt;\">" << std::fixed << std::setw(6) << std::setprecision(3) << centric_first() << "</value>" << std::endl;
                ss << "       <untwinned id=\"&lt;E&gt;\">" << std::fixed << std::setw(6) << std::setprecision(3) << theo_untwinned_centric_first() << " </untwinned>" << std::endl;
                ss << "       <twinned id=\"&lt;E&gt;\">" << std::fixed << std::setw(6) << std::setprecision(3) << theo_perfect_centric_first() << " </twinned>" << std::endl;
                ss << "    </Moment>" << std::endl;
                ss << "    <Moment id=\"&lt;|E^2-1|&gt;\">" << std::fixed << std::setw(6) << std::setprecision(3) << std::endl;
                ss << "       <value id=\"&lt;|E^2-1|&gt;\">" << std::fixed << std::setw(6) << std::setprecision(3) << centric_variance() << "</value>" << std::endl;
                ss << "       <untwinned id=\"&lt;|E^2-1|&gt;\">" << std::fixed << std::setw(6) << std::setprecision(3) << theo_untwinned_centric_variance() << " </untwinned>" << std::endl;
                ss << "       <twinned id=\"&lt;|E^2-1|&gt;\">" << std::fixed << std::setw(6) << std::setprecision(3) << theo_perfect_centric_variance() << " </twinned>" << std::endl;
                ss << "    </Moment>" << std::endl;
                
                ss << "    <Moment id=\"&lt;F^3&gt;/&lt;F&gt;^3\">" << std::endl;
                ss << "       <value id=\"&lt;F^3&gt;/&lt;F&gt;^3\">" << std::fixed << std::setw(6) << std::setprecision(3) << centric_third() << "</value>" << std::endl;
                ss << "       <untwinned id=\"&lt;F^3&gt;/&lt;F&gt;^3\">" << std::fixed << std::setw(6) << std::setprecision(3) << theo_untwinned_centric_third() << " </untwinned>" << std::endl;
                ss << "       <twinned id=\"&lt;F^3&gt;/&lt;F&gt;^3\">" << std::fixed << std::setw(6) << std::setprecision(3) << theo_perfect_centric_third() << " </twinned>" << std::endl;
                ss << "    </Moment>" << std::endl;
                ss << "    <Moment id=\"&lt;Z^2&gt;\">" << std::endl;
                ss << "       <value id=\"&lt;Z^2&gt;\">" << std::fixed << std::setw(6) << std::setprecision(3) << centric_fourth() << "</value>" << std::endl;
                ss << "       <untwinned id=\"&lt;Z^2&gt;\">" << std::fixed << std::setw(6) << std::setprecision(3) << theo_untwinned_centric_fourth() << " </untwinned>" << std::endl;
                ss << "       <twinned id=\"&lt;Z^2&gt;\">" << std::fixed << std::setw(6) << std::setprecision(3) << theo_perfect_centric_fourth() << " </twinned>" << std::endl;
                ss << "    </Moment>" << std::endl;
            }
            ss << "  </Moments>" << std::endl;
        }
        
		float invresolsq(1.0);
		if (_a_reso.size() != 1 && _c_reso.size() != 1) invresolsq = std::max(_a_reso[_a_reso.size()-1],_c_reso[_c_reso.size()-1]) + 0.01;
		else if (_a_reso.size() != 1) invresolsq = _a_reso[_a_reso.size()-1] + 0.01;
		else if (_c_reso.size() != 1) invresolsq = _c_reso[_c_reso.size()-1] + 0.01;
		
        if (_a_reso.size() != 1) {
            ss << "<CCP4Table groupID=\"graphMoments\" id=\"acentricMoments\" title=\""<< ((is_intensity() ) ? "intensity" : "amplitudes") << " moments vs resolution\">" << std::endl;
            if (is_intensity() ) {
                ss << "<plot>" << std::endl;
                ss << "<title>Intensity second moment &lt;|I^2|&gt;: " << std::fixed << std::setprecision(2) << acentric_second() << "</title>" << std::endl;
                ss << "<xscale>oneoversqrt</xscale>" << std::endl;
                ss << "<yrange min=\"0\" max=\"None\"/>" << std::endl;
                ss << "<plotline xcol=\"1\" ycol=\"2\" >" <<  std::endl;
                ss << "<symbolsize>  0</symbolsize>" << std::endl;
                ss << "<linestyle>-</linestyle>" << std::endl;
                ss << "<colour>red</colour>" << std::endl;
                ss << "</plotline>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << theo_untwinned_acentric_second() << "\" y2=\"" << std::setw(9) << theo_untwinned_acentric_second() <<"\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << theo_perfect_acentric_second() << "\" y2=\"" << std::setw(9) << theo_perfect_acentric_second() <<"\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << acentric_second() << "\" y2=\"" << std::setw(9) << acentric_second() <<"\" linestyle=\"-\" linecolour=\"red\"/>" << std::endl;
                ss << "</plot>" << std::endl;
                ss << "<plot>" << std::endl;
                ss << "<title>Intensity third moment &lt;|I^3|&gt;: " << std::fixed << std::setprecision(2) << acentric_third() << "</title>" << std::endl;
                ss << "<xscale>oneoversqrt</xscale>" << std::endl;
                ss << "<yrange min=\"0\" max=\"None\"/>" << std::endl;
                ss << "<plotline xcol=\"1\" ycol=\"3\" >" <<  std::endl;
                ss << "<symbolsize>  0</symbolsize>" << std::endl;
                ss << "<linestyle>-</linestyle>" << std::endl;
                ss << "<colour>red</colour>" << std::endl;
                ss << "</plotline>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << theo_untwinned_acentric_third() << "\" y2=\"" << std::setw(9) << theo_untwinned_acentric_third() <<"\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << theo_perfect_acentric_third() << "\" y2=\"" << std::setw(9) << theo_perfect_acentric_third() <<"\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << acentric_third() << "\" y2=\"" << std::setw(9) << acentric_third() <<"\" linestyle=\"-\" linecolour=\"red\"/>" << std::endl;
                ss << "</plot>" << std::endl;
                ss << "<plot>" << std::endl;
                ss << "<title>Intensity fourth moment &lt;|I^4|&gt;: " << std::fixed << std::setprecision(2) << acentric_fourth() << "</title>" << std::endl;
                ss << "<xscale>oneoversqrt</xscale>" << std::endl;
                ss << "<yrange min=\"0\" max=\"None\"/>" << std::endl;
                ss << "<plotline xcol=\"1\" ycol=\"4\" >" <<  std::endl;
                ss << "<symbolsize>  0</symbolsize>" << std::endl;
                ss << "<linestyle>-</linestyle>" << std::endl;
                ss << "<colour>red</colour>" << std::endl;
                ss << "</plotline>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << theo_untwinned_acentric_fourth() << "\" y2=\"" << std::setw(9) << theo_untwinned_acentric_fourth() <<"\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << theo_perfect_acentric_fourth() << "\" y2=\"" << std::setw(9) << theo_perfect_acentric_fourth() <<"\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << acentric_fourth() << "\" y2=\"" << std::setw(9) << acentric_fourth() <<"\" linestyle=\"-\" linecolour=\"red\"/>" << std::endl;
                ss << "</plot>" << std::endl;
                ss << "<headers separator=\" \">\n 1/resol^2 &lt;I^2&gt; &lt;I^3&gt; &lt;I^4&gt;\n </headers>" << std::endl;
            } else {
                ss << "<plot>" << std::endl;
                ss << "<title>Amplitude first moment &lt;|E|&gt;: " << std::fixed << std::setprecision(2) << acentric_first() << "</title>" << std::endl;
                ss << "<xscale>oneoversqrt</xscale>" << std::endl;
                ss << "<yrange min=\"0\" max=\"None\"/>" << std::endl;
                ss << "<plotline xcol=\"1\" ycol=\"2\" >" <<  std::endl;
                ss << "<symbolsize>  0</symbolsize>" << std::endl;
                ss << "<linestyle>-</linestyle>" << std::endl;
                ss << "<colour>red</colour>" << std::endl;
                ss << "</plotline>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << theo_untwinned_acentric_first() << "\" y2=\"" << std::setw(9) << theo_untwinned_acentric_first() <<"\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << theo_perfect_acentric_first() << "\" y2=\"" << std::setw(9) << theo_perfect_acentric_first() <<"\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << acentric_first() << "\" y2=\"" << std::setw(9) << acentric_first() <<"\" linestyle=\"-\" linecolour=\"red\"/>" << std::endl;
                ss << "</plot>" << std::endl;
                ss << "<plot>" << std::endl;
                ss << "<title>Amplitude variance &lt;|E^2-1|&gt;: " << std::fixed << std::setprecision(2) << acentric_variance() << "</title>" << std::endl;
                ss << "<xscale>oneoversqrt</xscale>" << std::endl;
                ss << "<yrange min=\"0\" max=\"None\"/>" << std::endl;
                ss << "<plotline xcol=\"1\" ycol=\"2\" >" <<  std::endl;
                ss << "<symbolsize>  0</symbolsize>" << std::endl;
                ss << "<linestyle>-</linestyle>" << std::endl;
                ss << "<colour>red</colour>" << std::endl;
                ss << "</plotline>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << theo_untwinned_acentric_variance() << "\" y2=\"" << std::setw(9) << theo_untwinned_acentric_variance() <<"\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << theo_perfect_acentric_first() << "\" y2=\"" << std::setw(9) << theo_perfect_acentric_first() <<"\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << acentric_variance() << "\" y2=\"" << std::setw(9) << acentric_variance() <<"\" linestyle=\"-\" linecolour=\"red\"/>" << std::endl;
                ss << "</plot>" << std::endl;
                ss << "<plot>" << std::endl;
                ss << "<title>Amplitude third moment &lt;|E^3|&gt;: " << std::fixed << std::setprecision(2) << acentric_first() << "</title>" << std::endl;
                ss << "<xscale>oneoversqrt</xscale>" << std::endl;
                ss << "<yrange min=\"0\" max=\"None\"/>" << std::endl;
                ss << "<plotline xcol=\"1\" ycol=\"2\" >" <<  std::endl;
                ss << "<symbolsize>  0</symbolsize>" << std::endl;
                ss << "<linestyle>-</linestyle>" << std::endl;
                ss << "<colour>red</colour>" << std::endl;
                ss << "</plotline>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << theo_untwinned_acentric_third() << "\" y2=\"" << std::setw(9) << theo_untwinned_acentric_third() <<"\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << theo_perfect_acentric_third() << "\" y2=\"" << std::setw(9) << theo_perfect_acentric_third() <<"\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << acentric_third() << "\" y2=\"" << std::setw(9) << acentric_first() <<"\" linestyle=\"-\" linecolour=\"red\"/>" << std::endl;
                ss << "</plot>" << std::endl;
                ss << "<title>Amplitude fourth moment &lt;|I^3|&gt;: " << std::fixed << std::setprecision(2) << acentric_first() << "</title>" << std::endl;
                ss << "<xscale>oneoversqrt</xscale>" << std::endl;
                ss << "<yrange min=\"0\" max=\"None\"/>" << std::endl;
                ss << "<plotline xcol=\"1\" ycol=\"2\" >" <<  std::endl;
                ss << "<symbolsize>  0</symbolsize>" << std::endl;
                ss << "<linestyle>-</linestyle>" << std::endl;
                ss << "<colour>red</colour>" << std::endl;
                ss << "</plotline>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << theo_untwinned_acentric_fourth() << "\" y2=\"" << std::setw(9) << theo_untwinned_acentric_fourth() <<"\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << theo_perfect_acentric_fourth() << "\" y2=\"" << std::setw(9) << theo_perfect_acentric_fourth() <<"\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << acentric_fourth() << "\" y2=\"" << std::setw(9) << acentric_fourth() <<"\" linestyle=\"-\" linecolour=\"red\"/>" << std::endl;
                ss << "</plot>" << std::endl;
                ss << "<headers separator=\" \">\n  1/resol^2 &lt;E&gt; &lt;|E^2-1|&gt; &lt;E^3&gt; &lt;E^4&gt;\n  </headers>" << std::endl;
            }
            
            ss << "<data>" << std::endl;;
            for(int i=0;i!=_a_reso.size();++i){
                if (is_intensity() )  {
                    ss << std::fixed << std::setw(8) << std::setprecision(4) << _a_reso[i]<< " " << I2a[i] << " " << " " <<  I3a[i] << " "  << I4a[i] << std::endl;
                } else {
                    ss << std::fixed << std::setw(8) << std::setprecision(4) << _a_reso[i]<< " " << I1a[i] << " " <<  I5a[i] << " " << " " <<  I3a[i] << " "  << I4a[i] << std::endl;
                }
            }
            ss << "</data>" << std::endl;
            ss << "</CCP4Table>" << std::endl;
        }
        if (_c_reso.size() != 1) {
            ss << "<CCP4Table groupID=\"graphMoments\" id=\"centricMoments\" title=\""<< ((is_intensity() ) ? "intensity" : "amplitudes") << " moments vs resolution\">" << std::endl;
            if (is_intensity() ) {
                ss << "<plot>" << std::endl;
                ss << "<title>Intensity second moment &lt;|I^2|&gt;: " << std::fixed << std::setprecision(2) << centric_second() << "</title>" << std::endl;
                ss << "<xscale>oneoversqrt</xscale>" << std::endl;
                ss << "<yrange min=\"0\" max=\"None\"/>" << std::endl;
                ss << "<plotline xcol=\"1\" ycol=\"2\" >" <<  std::endl;
                ss << "<symbolsize>  0</symbolsize>" << std::endl;
                ss << "<linestyle>-</linestyle>" << std::endl;
                ss << "<colour>red</colour>" << std::endl;
                ss << "</plotline>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << theo_untwinned_centric_second() << "\" y2=\"" << std::setw(9) << theo_untwinned_centric_second() <<"\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << theo_perfect_centric_second() << "\" y2=\"" << std::setw(9) << theo_perfect_centric_second() <<"\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << centric_second() << "\" y2=\"" << std::setw(9) << centric_second() <<"\" linestyle=\"-\" linecolour=\"red\"/>" << std::endl;
                ss << "</plot>" << std::endl;
                ss << "<plot>" << std::endl;
                ss << "<title>Intensity third moment &lt;|I^3|&gt;: " << std::fixed << std::setprecision(2) << centric_third() << "</title>" << std::endl;
                ss << "<xscale>oneoversqrt</xscale>" << std::endl;
                ss << "<yrange min=\"0\" max=\"None\"/>" << std::endl;
                ss << "<plotline xcol=\"1\" ycol=\"3\" >" <<  std::endl;
                ss << "<symbolsize>  0</symbolsize>" << std::endl;
                ss << "<linestyle>-</linestyle>" << std::endl;
                ss << "<colour>red</colour>" << std::endl;
                ss << "</plotline>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << theo_untwinned_centric_third() << "\" y2=\"" << std::setw(9) << theo_untwinned_centric_third() <<"\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << theo_perfect_centric_third() << "\" y2=\"" << std::setw(9) << theo_perfect_centric_third() <<"\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << centric_third() << "\" y2=\"" << std::setw(9) << centric_third() <<"\" linestyle=\"-\" linecolour=\"red\"/>" << std::endl;
                ss << "</plot>" << std::endl;
                ss << "<plot>" << std::endl;
                ss << "<title>Intensity fourth moment &lt;|I^4|&gt;: " << std::fixed << std::setprecision(2) << centric_fourth() << "</title>" << std::endl;
                ss << "<xscale>oneoversqrt</xscale>" << std::endl;
                ss << "<yrange min=\"0\" max=\"None\"/>" << std::endl;
                ss << "<plotline xcol=\"1\" ycol=\"4\" >" <<  std::endl;
                ss << "<symbolsize>  0</symbolsize>" << std::endl;
                ss << "<linestyle>-</linestyle>" << std::endl;
                ss << "<colour>red</colour>" << std::endl;
                ss << "</plotline>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << theo_untwinned_centric_fourth() << "\" y2=\"" << std::setw(9) << theo_untwinned_centric_fourth() <<"\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << theo_perfect_centric_fourth() << "\" y2=\"" << std::setw(9) << theo_perfect_centric_fourth() <<"\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << centric_fourth() << "\" y2=\"" << std::setw(9) << centric_fourth() <<"\" linestyle=\"-\" linecolour=\"red\"/>" << std::endl;
                ss << "</plot>" << std::endl;
                ss << "<headers separator=\" \">\n 1/resol^2 &lt;I^2&gt; &lt;I^3&gt; &lt;I^4&gt;\n </headers>" << std::endl;
            } else {
                ss << "<plot>" << std::endl;
                ss << "<title>Amplitude first moment &lt;|E|&gt;: " << std::fixed << std::setprecision(2) << centric_first() << "</title>" << std::endl;
                ss << "<xscale>oneoversqrt</xscale>" << std::endl;
                ss << "<yrange min=\"0\" max=\"None\"/>" << std::endl;
                ss << "<plotline xcol=\"1\" ycol=\"2\" >" <<  std::endl;
                ss << "<symbolsize>  0</symbolsize>" << std::endl;
                ss << "<linestyle>-</linestyle>" << std::endl;
                ss << "<colour>red</colour>" << std::endl;
                ss << "</plotline>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << theo_untwinned_centric_first() << "\" y2=\"" << std::setw(9) << theo_untwinned_centric_first() <<"\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << theo_perfect_centric_first() << "\" y2=\"" << std::setw(9) << theo_perfect_centric_first() <<"\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << centric_first() << "\" y2=\"" << std::setw(9) << centric_first() <<"\" linestyle=\"-\" linecolour=\"red\"/>" << std::endl;
                ss << "</plot>" << std::endl;
                ss << "<plot>" << std::endl;
                ss << "<title>Amplitude variance &lt;|E^2-1|&gt;: " << std::fixed << std::setprecision(2) << centric_variance() << "</title>" << std::endl;
                ss << "<xscale>oneoversqrt</xscale>" << std::endl;
                ss << "<yrange min=\"0\" max=\"None\"/>" << std::endl;
                ss << "<plotline xcol=\"1\" ycol=\"2\" >" <<  std::endl;
                ss << "<symbolsize>  0</symbolsize>" << std::endl;
                ss << "<linestyle>-</linestyle>" << std::endl;
                ss << "<colour>red</colour>" << std::endl;
                ss << "</plotline>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << theo_untwinned_centric_variance() << "\" y2=\"" << std::setw(9) << theo_untwinned_centric_variance() <<"\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << theo_perfect_centric_first() << "\" y2=\"" << std::setw(9) << theo_perfect_centric_first() <<"\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << centric_variance() << "\" y2=\"" << std::setw(9) << centric_variance() <<"\" linestyle=\"-\" linecolour=\"red\"/>" << std::endl;
                ss << "</plot>" << std::endl;
                ss << "<plot>" << std::endl;
                ss << "<title>Amplitude third moment &lt;|E^3|&gt;: " << std::fixed << std::setprecision(2) << centric_first() << "</title>" << std::endl;
                ss << "<xscale>oneoversqrt</xscale>" << std::endl;
                ss << "<yrange min=\"0\" max=\"None\"/>" << std::endl;
                ss << "<plotline xcol=\"1\" ycol=\"2\" >" <<  std::endl;
                ss << "<symbolsize>  0</symbolsize>" << std::endl;
                ss << "<linestyle>-</linestyle>" << std::endl;
                ss << "<colour>red</colour>" << std::endl;
                ss << "</plotline>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << theo_untwinned_centric_third() << "\" y2=\"" << std::setw(9) << theo_untwinned_centric_third() <<"\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << theo_perfect_centric_third() << "\" y2=\"" << std::setw(9) << theo_perfect_centric_third() <<"\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << centric_third() << "\" y2=\"" << std::setw(9) << centric_first() <<"\" linestyle=\"-\" linecolour=\"red\"/>" << std::endl;
                ss << "</plot>" << std::endl;
                ss << "<title>Amplitude fourth moment &lt;|I^3|&gt;: " << std::fixed << std::setprecision(2) << centric_first() << "</title>" << std::endl;
                ss << "<xscale>oneoversqrt</xscale>" << std::endl;
                ss << "<yrange min=\"0\" max=\"None\"/>" << std::endl;
                ss << "<plotline xcol=\"1\" ycol=\"2\" >" <<  std::endl;
                ss << "<symbolsize>  0</symbolsize>" << std::endl;
                ss << "<linestyle>-</linestyle>" << std::endl;
                ss << "<colour>red</colour>" << std::endl;
                ss << "</plotline>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << theo_untwinned_centric_fourth() << "\" y2=\"" << std::setw(9) << theo_untwinned_centric_fourth() <<"\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << theo_perfect_centric_fourth() << "\" y2=\"" << std::setw(9) << theo_perfect_centric_fourth() <<"\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"" << std::fixed << std::setw(8) << std::setprecision(4) << invresolsq << "\" y1=\"" << std::fixed << std::setw(9) << std::setprecision(4) << centric_fourth() << "\" y2=\"" << std::setw(9) << centric_fourth() <<"\" linestyle=\"-\" linecolour=\"red\"/>" << std::endl;
                ss << "</plot>" << std::endl;
                ss << "<headers separator=\" \">\n  1/resol^2 &lt;E&gt; &lt;|E^2-1|&gt; &lt;E^3&gt; &lt;E^4&gt;\n  </headers>" << std::endl;
            }
            
            ss << "<data>" << std::endl;;
            for(int i=0;i!=_c_reso.size();++i){
                if (is_intensity() )  {
                    ss << std::fixed << std::setw(8) << std::setprecision(4) << _c_reso[i]<< " " << I2c[i] << " " << " " <<  I3c[i] << " "  << I4c[i] << std::endl;
                } else {
                    ss << std::fixed << std::setw(8) << std::setprecision(4) << _c_reso[i]<< " " << I1c[i] << " " <<  I5c[i] << " " << " " <<  I3c[i] << " "  << I4c[i] << std::endl;
                }
            }
            ss << "</data>" << std::endl;
            ss << "</CCP4Table>" << std::endl;
        }
        ss << "<Moments>" << std::endl;
        ss << "  <Comment id='MomentsReso'>" << std::endl;
        ss << "The acentric moments plot can also give an indication of the quality of the high resolution data.  Large deviations often indicate bias in the data." << std::endl;
        ss << "  </Comment>" << std::endl;
        ss << "</Moments>" << std::endl;
        return ss;
    }
    
    //-------instantiate templates----------------------------------------------
    
    //template void Moments::operator()<float, clipper::datatypes::I_sigI>(const clipper::HKL_data<clipper::datatypes::I_sigI<float> >&,const clipper::Range<clipper::ftype>&);
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
