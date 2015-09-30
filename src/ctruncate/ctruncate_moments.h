//
//     CTRUNCATE
//     Copyright (C) 2006-2011 Norman Stein, Charles Ballard
//
//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Application.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//

#ifndef __CTRUNCATE_MOMENTS_H
#define __CTRUNCATE_MOMENTS_H

#include "clipper/clipper.h"
#include "clipper/clipper-ccp4.h"
#include "intensity_target.h"

namespace ctruncate {
			
	class Moments {
	public:
		Moments() {}
		template <class T, template <class> class D> Moments( clipper::HKL_data<D<T> >&);
		template <class T, template <class> class D> Moments( clipper::HKL_data<D<T> >&, clipper::Range<clipper::ftype>&);
		~Moments();
		
		template <class T, template <class> class D> void operator()( clipper::HKL_data<D<T> >&,clipper::Range<clipper::ftype>&);
		
		const clipper::String type() { return _base->type(); }
		
		inline bool is_intensity() { return (type().compare("F_sigF") || type().compare("F_sigF_ano") ); }
		
		void resolution(const clipper::Range<clipper::ftype>& reso) { calc_range = reso; for (int i=0; i !=5; ++i ) acache[i]=ccache[i]=0.0;}
		
		clipper::ftype acentric_first();
		clipper::ftype acentric_second();
		clipper::ftype acentric_third();
		clipper::ftype acentric_fourth();
		clipper::ftype acentric_variance();
		int num_acentric() { return na; }
		
		clipper::ftype centric_first();
		clipper::ftype centric_second();
		clipper::ftype centric_third();
		clipper::ftype centric_fourth();
		clipper::ftype centric_variance();
		int num_centric() { return nc; }
		
		clipper::ftype fraction();
		
		void summary();
		
		void loggraph();
		
        void output();
        
		std::stringstream& xml_output(std::stringstream&);
		
		clipper::ftype theo_untwinned_acentric_first();
		clipper::ftype theo_perfect_acentric_first();
		clipper::ftype theo_untwinned_acentric_second();
		clipper::ftype theo_perfect_acentric_second();
		clipper::ftype theo_untwinned_acentric_third();
		clipper::ftype theo_perfect_acentric_third();
		clipper::ftype theo_untwinned_acentric_fourth();
		clipper::ftype theo_perfect_acentric_fourth();
		clipper::ftype theo_untwinned_acentric_variance();
        clipper::ftype theo_perfect_acentric_variance();
		
		clipper::ftype theo_untwinned_centric_first();
		clipper::ftype theo_perfect_centric_first();
		clipper::ftype theo_untwinned_centric_second();
		clipper::ftype theo_perfect_centric_second();
		clipper::ftype theo_untwinned_centric_third();
		clipper::ftype theo_perfect_centric_third();
		clipper::ftype theo_untwinned_centric_fourth();
		clipper::ftype theo_perfect_centric_fourth();
		clipper::ftype theo_untwinned_centric_variance();
        clipper::ftype theo_perfect_centric_variance();
		
	private:
		static const clipper::ftype E_THEO_UNTWINNED_CENTRIC_FIRST=0.798;
		static const clipper::ftype E_THEO_PERFECT_CENTRIC_FIRST=0.886;
		static const clipper::ftype E_THEO_PERFECT_ACENTRIC_FIRST=0.94;
		static const clipper::ftype E_THEO_UNTWINNED_CENTRIC_SECOND=0.0;
		static const clipper::ftype E_THEO_PERFECT_CENTRIC_SECOND=0.0;
		static const clipper::ftype E_THEO_PERFECT_ACENTRIC_SECOND=0.0;
		static const clipper::ftype E_THEO_UNTWINNED_CENTRIC_THIRD=1.596;
		static const clipper::ftype E_THEO_PERFECT_CENTRIC_THIRD=1.329;
		static const clipper::ftype E_THEO_PERFECT_ACENTRIC_THIRD=1.175;
		static const clipper::ftype E_THEO_UNTWINNED_CENTRIC_FOURTH=3.0;
		static const clipper::ftype E_THEO_PERFECT_CENTRIC_FOURTH=2.0;
		static const clipper::ftype E_THEO_PERFECT_ACENTRIC_FOURTH=1.5;
		static const clipper::ftype E_THEO_UNTWINNED_CENTRIC_VAR=0.968;
		static const clipper::ftype E_THEO_UNTWINNED_ACENTRIC_VAR=0.736;
		static const clipper::ftype Z_THEO_UNTWINNED_ACENTRIC_VAR=0.0;
		static const clipper::ftype Z_THEO_UNTWINNED_CENTRIC_FIRST=0.0;
		static const clipper::ftype Z_THEO_PERFECT_CENTRIC_FIRST=0.0;
		static const clipper::ftype Z_THEO_PERFECT_ACENTRIC_FIRST=0.0;
		static const clipper::ftype Z_THEO_UNTWINNED_CENTRIC_THIRD=15.0;
		static const clipper::ftype Z_THEO_PERFECT_CENTRIC_THIRD=6.0;
		static const clipper::ftype Z_THEO_PERFECT_ACENTRIC_THIRD=3.0;
		static const clipper::ftype Z_THEO_UNTWINNED_CENTRIC_FOURTH=105.0;
		static const clipper::ftype Z_THEO_PERFECT_CENTRIC_FOURTH=24.0;
		static const clipper::ftype Z_THEO_PERFECT_ACENTRIC_FOURTH=7.5;
		static const clipper::ftype E_THEO_PERFECT_ACENTRIC_VAR=0.541;
		
		const clipper::HKL_data_base *_base; //!< data
		clipper::HKL_data<clipper::data64::I_sigI> _normalised; //!< store norm
		//clipper::Generic_ordinal _ord_c;  //!< store centric bin resolution 
		//clipper::Generic_ordinal _ord_a;  //!< store acentric bin resolution
		clipper::Range<clipper::ftype> calc_range;
		//clipper::BasisFn_base* acentric_norm_base;
		//std::vector<clipper::ftype> acentric_norm;
		//clipper::BasisFn_base* centric_norm_base;
		//std::vector<clipper::ftype> centric_norm;
		int na; //!< number acentric in sets
		int nc; //!< number centric in sets
		//int abins;
		//int cbins;
		std::vector<clipper::ftype> _a_reso;  //!< acentric resolution bins
		std::vector<clipper::ftype> _c_reso;  //!< centric resolution bins
		std::vector<clipper::ftype> I1a, I2a, I3a, I4a, I5a; //!< store acentric moment data
		std::vector<clipper::ftype> I1c, I2c, I3c, I4c, I5c; //!< store centric moment data
		
		//bool mode;  //!< dealing with F or I
		//clipper::ftype pn;  //!< power associated with intensity or amplitudes
		
		clipper::ftype acache[5];
		clipper::ftype ccache[5];
		
		//template<class T, template<class> class D > void norm(const clipper::HKL_data<D<T> >& );
		template <class T, template<class> class D> void init(const clipper::HKL_data<D<T> >&, clipper::Range<clipper::ftype>& );
		clipper::ftype calc(clipper::ftype pow, bool centric);
		clipper::ftype calc_var(bool centric);
		inline clipper::ftype pn() { return ( is_intensity() ) ? 1.0 : 2.0; }

		template <class T > const T&    obs( const clipper::datatypes::F_sigF<T>& f ) const { return f.f(); }
		template <class T > const T&    obs( const clipper::datatypes::I_sigI<T>& f ) const { return f.I(); }
		template <class T > const T&    obs( const clipper::datatypes::F_sigF_ano<T>& f ) const { return f.f(); }
		template <class T > const T&    obs( const clipper::datatypes::I_sigI_ano<T>& f ) const { return f.I(); }
		template <class T > const T& sigobs( const clipper::datatypes::F_sigF<T>& f ) { return f.sigf(); }
		template <class T > const T& sigobs( const clipper::datatypes::I_sigI<T>& f ) { return f.sigI(); }		
		template <class T > const T& sigobs( const clipper::datatypes::F_sigF_ano<T>& f ) { return f.sigf(); }
		template <class T > const T& sigobs( const clipper::datatypes::I_sigI_ano<T>& f ) { return f.sigI(); }		
		
	};
	
	template <class T, template<class> class D> Moments::Moments(clipper::HKL_data<D<T> >& data) {
		clipper::Range<clipper::ftype> range;
		init(data,range);
	}
	
	template <class T, template<class> class D> Moments::Moments(clipper::HKL_data<D<T> >& data, clipper::Range<clipper::ftype>& range) {	
		init(data,range);
	 }
	
	template <class T, template<class> class D> void Moments::operator()(clipper::HKL_data<D<T> >& data, clipper::Range<clipper::ftype>& range) {
		calc_range = range;
		init(data,range);
	}	
	
	
//#include "intensity_target.h"
	
	template <class T, template<class> class D> void Moments::init(const clipper::HKL_data<D<T> >& data,clipper::Range<clipper::ftype>& range) {
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		_base = &data;
		calc_range = range;
		for (int i=0; i !=5; ++i ) acache[i]=ccache[i]=0.0;
		na=nc=0;
		
		clipper::HKL_data<D<T> > fa(data.hkl_info() ), fc(data.hkl_info() );
		for ( HRI ih = data.first(); !ih.last(); ih.next() ) {
			if ( !data[ih].missing() ) {
				if ( !ih.hkl_class().centric() ) {
					++na;
					fa[ih] = D<T>(obs(data[ih]),sigobs(data[ih]) );
				} else {
					++nc;
					fc[ih] = D<T>(obs(data[ih]),sigobs(data[ih]) );
				}
			}
		}
		
		clipper::Resolution_ordinal _ord_a;
		_ord_a.init(fa,1.0);
		clipper::Resolution_ordinal _ord_c;
		_ord_c.init(fc,1.0);
		
		int nreflns(500);
		int abins = std::max(1,na/nreflns);
		int cbins = std::max(1,nc/nreflns);
        if ( abins < 10 ) abins = std::max(1,na/200);
        if ( cbins < 10 ) cbins = std::max(1,nc/200);
        if ( abins < 10 ) abins = std::max(1,na/100);
        if ( cbins < 10 ) cbins = std::max(1,nc/100);
        
		_a_reso.resize(abins,0.0);
		_c_reso.resize(cbins,0.0);
		
		std::vector<clipper::ftype> sumova(abins,0.0), sumovc(cbins,0.0);
		for ( HRI ih = data.first(); !ih.last(); ih.next() ) {
			if ( !data[ih].missing() ) {
				//clipper::ftype eps = ( mode ) ? 1.0/ih.hkl_class().epsilonc() : 1.0/std::sqrt(ih.hkl_class().epsilonc());
				clipper::ftype s = ih.invresolsq();
				int bin(0);
				if ( !ih.hkl_class().centric() ) {
					bin = clipper::Util::bound( 0,clipper::Util::intf( clipper::ftype(abins) * _ord_a.ordinal( s ) ), abins-1 );
					_a_reso[bin] += s;
					sumova[bin] += 1.0;
				} else {
					bin = clipper::Util::bound( 0,clipper::Util::intf( clipper::ftype(cbins) * _ord_c.ordinal( s ) ), cbins-1 );
					_c_reso[bin] += s;
					sumovc[bin] += 1.0;
				}
			}
		}
		
		for (int i=0 ; i != abins; ++i ) _a_reso[i]/=sumova[i];
		for (int i=0 ; i != cbins; ++i ) _c_reso[i]/=sumovc[i];
		
		std::vector<clipper::ftype> params_init( abins, 1.0 );
		clipper::BasisFn_binner acentric_norm(fa, abins, 1.0);
		clipper::ftype pn=this->pn();
		
		TargetFn_meanInth<D<T> >  atarget_fn( fa, pn );
		clipper::ResolutionFn af( data.hkl_info() , acentric_norm, atarget_fn, params_init ); 
		
		clipper::BasisFn_binner centric_norm(fc, cbins, 1.0);
		
		TargetFn_meanInth<D<T> > ctarget_fn( fc, pn );
		clipper::ResolutionFn cf( data.hkl_info() , centric_norm, ctarget_fn, params_init );
		
		
		_normalised.init(data );
		
		pn=1.0/this->pn();
		for ( HRI ih = data.first(); !ih.last(); ih.next() ) {
			if ( !data[ih].missing() ) {
				if ( !ih.hkl_class().centric() ) {
					//na += ih.hkl_class().epsilonc();
					_normalised[ih] = clipper::data64::I_sigI((double)obs(data[ih])/std::pow(acentric_norm.f_s(ih.invresolsq(), af.params() ),pn),(double)sigobs(data[ih]) );
				} else {
					//nc += ih.hkl_class().epsilonc();
					_normalised[ih] = clipper::data64::I_sigI((double)obs(data[ih])/std::pow(centric_norm.f_s(ih.invresolsq(), cf.params() ),pn),(double)sigobs(data[ih]) );
				}
			}
		} 
		
		I1a.resize(abins,0.0);
		I1c.resize(cbins,0.0);
		I2a.resize(abins,0.0);
		I2c.resize(cbins,0.0);
		I3a.resize(abins,0.0); 
		I3c.resize(cbins,0.0);
		I4a.resize(abins,0.0); 
		I4c.resize(cbins,0.0);
		I5a.resize(abins,0.0);
		I5c.resize(cbins,0.0);
				
		bool is(is_intensity() );
		for ( HRI ih = _normalised.first(); !ih.last(); ih.next() ) {
			if ( !_normalised[ih].missing() ) {
				clipper::ftype eps = ( is ) ? 1.0/ih.hkl_class().epsilon() : 1.0/std::sqrt(ih.hkl_class().epsilon());
				double val = eps*obs(_normalised[ih]);
				clipper::ftype s = ih.invresolsq();
				int bin = clipper::Util::bound( 0,clipper::Util::intf( clipper::ftype(abins) * _ord_a.ordinal( s ) ), abins-1 );
				int cbin = clipper::Util::bound( 0,clipper::Util::intf( clipper::ftype(cbins) * _ord_c.ordinal( s ) ), cbins-1 );
				if (!ih.hkl_class().centric()) {
					if (val > 0.0) {
						I1a[bin] += val;
						I2a[bin] += val*val;
						I3a[bin] += val*val*val;
						I4a[bin] += val*val*val*val;
					}
				}
				else if (cbins != 0) {
					if (val > 0.0) {
						I1c[cbin] += val;
						I2c[cbin] += val*val;
						I3c[cbin] += val*val*val;
						I4c[cbin] += val*val*val*val;
					}
				}
			}
		}
		
		for ( int bin=0 ; bin != abins ; bin++ ) {
			if ( sumova[bin] != 0 ) {
				I1a[bin] /= sumova[bin];
				I2a[bin] /= sumova[bin];
				I3a[bin] /= sumova[bin];
				I4a[bin] /= sumova[bin];
				I5a[bin] = I2a[bin] - I1a[bin]*I1a[bin];	
			}
		}
		
		for ( int bin=0 ; bin != cbins ; bin++ ) {
			if ( sumovc[bin] != 0 ) {
				I1c[bin] /= sumovc[bin];
				I2c[bin] /= sumovc[bin];
				I3c[bin] /= sumovc[bin];
				I4c[bin] /= sumovc[bin];
				I5c[bin] = I2c[bin] - I1c[bin]*I1c[bin];				
			} 
		}

		calc(1.0, false);
		calc(2.0, false);
		calc(3.0, false);
		calc(4.0, false);
		calc_var(false);
		calc(1.0, true);
		calc(2.0, true);
		calc(3.0, true);
		calc(4.0, true);
		calc_var(true);
	}
	
}

#endif

