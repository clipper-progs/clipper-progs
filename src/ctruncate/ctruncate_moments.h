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

namespace ctruncate {
		
	template<class D> class Moments {
	public:
		Moments( const clipper::HKL_data<D>&);
		Moments( const clipper::HKL_data<D>&, const clipper::Range<clipper::ftype>&);
		~Moments();
		
		const clipper::String type() { return hkl_data->type(); }
		
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
		
		clipper::ftype theo_untwinned_acentric_first();
		clipper::ftype theo_perfect_acentric_first();
		clipper::ftype theo_untwinned_acentric_second();
		clipper::ftype theo_perfect_acentric_second();
		clipper::ftype theo_untwinned_acentric_third();
		clipper::ftype theo_perfect_acentric_third();
		clipper::ftype theo_untwinned_acentric_fourth();
		clipper::ftype theo_perfect_acentric_fourth();
		clipper::ftype theo_untwinned_acentric_variance();
		
		clipper::ftype theo_untwinned_centric_first();
		clipper::ftype theo_perfect_centric_first();
		clipper::ftype theo_untwinned_centric_second();
		clipper::ftype theo_perfect_centric_second();
		clipper::ftype theo_untwinned_centric_third();
		clipper::ftype theo_perfect_centric_third();
		clipper::ftype theo_untwinned_centric_fourth();
		clipper::ftype theo_perfect_centric_fourth();
		clipper::ftype theo_untwinned_centric_variance();
		
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

		const clipper::HKL_data<D>* hkl_data;
		clipper::HKL_data<D> normalised;
		clipper::Range<clipper::ftype> calc_range;
		//clipper::BasisFn_base* acentric_norm_base;
		//std::vector<clipper::ftype> acentric_norm;
		//clipper::BasisFn_base* centric_norm_base;
		//std::vector<clipper::ftype> centric_norm;
		int na;
		int nc;
		int abins;
		int cbins;
		
		bool mode;
		clipper::ftype pn;
		
		clipper::ftype acache[5];
		clipper::ftype ccache[5];
		
		void norm(const clipper::HKL_data<D>& );
		clipper::ftype calc(clipper::ftype pow, bool centric);
		clipper::ftype calc_var(bool centric);
		
		template <class T > const T&    obs( const clipper::datatypes::F_sigF<T>& f ) const { return f.f(); }
		template <class T > const T&    obs( const clipper::datatypes::I_sigI<T>& f ) const { return f.I(); }
		template <class T > const T& sigobs( const clipper::datatypes::F_sigF<T>& f ) { return f.sigf(); }
		template <class T > const T& sigobs( const clipper::datatypes::I_sigI<T>& f ) { return f.sigI(); }		
		
	};
}

#endif

