//
//     CTRUNCATE
//     Copyright (C) 2006-2011 Norman Stein, Charles Ballard
//
//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Application.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//

#ifndef __CTRUNCATE_TWIN
#define __CTRUNCATE_TWIN

#include "clipper/clipper.h"
#include "clipper/clipper-ccp4.h"
#include <string>

namespace ctruncate {
	
    
    void twin_summary(clipper::ftype, clipper::ftype);
	
    //! compute the twinning operators
    /*! Calculate either by first principes or
        via tables
     */
    class TwinSymops {
    public:
        enum MODE {FP,TABLE}; //!< use first principles or tabulated twinning operators
        // constructor
        explicit TwinSymops(const clipper::Cell& cell, const clipper::Spacegroup& spg, MODE mode=FP );
        // return ith twinop as ISymop
        const clipper::Isymop& operator()(int i) const { return _twinops[i]; }
        const clipper::Isymop& operator[](int i) const { return _twinops[i]; }
        // return number of twinning operators
        int size() const { return _twinops.size(); }
        // return summary
        void summary() const;
        
        
    private:
        // set up tabulated symops
        void table(const clipper::Cell&, const clipper::Spacegroup&);
        // set up first principle symops
        void principles(const clipper::Cell&, const clipper::Spacegroup&);

        MODE _mode; //record mode
        std::vector<clipper::Isymop> _twinops; //!< twinning operators
    };
    
    //! Compute L-test scores on intensity
    /*! Compute the L-test, store results for summary and loggraph.
     Padilla and Yeates Acta Cryst. D59 1124 (2003) 
     */
    class L_test {
    public:
        // constructor, setup up L-test
        L_test( int nbins=20) : _cdf(nbins,0.0) {}
        // constructor, calculates L statistics if true
        template<class T> explicit L_test( clipper::HKL_data< clipper::datatypes::I_sigI<T> >& isig, clipper::Resolution reso = clipper::Resolution(3.0), int nbins=20);
        // perform calculation, or recalculate value
        template<class T> clipper::ftype operator()(clipper::HKL_data< clipper::datatypes::I_sigI<T> >& isig, clipper::Resolution reso = clipper::Resolution(3.0)) ;
        // return L-statistic
        clipper::ftype statistic() const { return _Lav; }
		// return L2-statistic
		clipper::ftype L2statistic() const { return _L2av; }
		// estimate alpha from L
		clipper::ftype estimateAlphafromL() const;
		// estimate alpha from L2
		clipper::ftype estimateAlphafromL2() const;
		// estimate alpha from N(L)
		clipper::ftype estimateAlphafromNL() const;
		// estimate twin fraction
		        const clipper::ftype fraction() { 
					if (_alpha == -1.0 ) _alpha = std::max(std::max(this->estimateAlphafromL(),this->estimateAlphafromL2() ),this->estimateAlphafromNL());
					return _alpha; 
				}
        // print summary
        void summary() const ;
        // print loggraph
        void loggraph() const ;
        
    private:
        clipper::Resolution _reso; //!< resolution used in calculation
        clipper::ftype _Lav; //!< L-statistic
        clipper::ftype _L2av; //!< L-statistic squared
		clipper::ftype _alpha; //!< alpha values from H-test
        clipper::ftype _NLT; //!< Number of data points
        std::vector<clipper::ftype> _cdf; //!< cumulative density function
        
        // perform calculation
        template<class T> clipper::ftype calc(clipper::HKL_data< clipper::datatypes::I_sigI<T> >& isig, clipper::Resolution& reso );
		//interpolate
		inline clipper::ftype interpolate(const clipper::ftype y, const clipper::ftype y1, const clipper::ftype y2, const clipper::ftype x1, const clipper::ftype x2) const {
			return	x1+(y-y1)*(x2-x1)/(y2-y1);
		}
		// mean L for alpha
		// P&Y eqn. 31
		inline clipper::ftype meanL(const clipper::ftype a) const {
			clipper::ftype a12 = 1.0-2.0*a;
			clipper::ftype a1 = 1.0-a;
			clipper::ftype a2 = a*a;
			return (a12*a12*(1.0-6.0*(a-a2))-8.0*a1*a1*a2*log(4.0*a*a1))/(2.0*a12*a12*a12*a12);
		}
		// mean L2 for alpha
		// P&Y eqn. 32
		inline clipper::ftype meanL2(const clipper::ftype a) const {
			clipper::ftype a12 = 1.0-2.0*a;
			clipper::ftype a1 = 1.0-a;
			clipper::ftype a2 = a*a;
			return (a12*(1.0-12.0*a-4.0*a2+32.0*a*a2-16*a2*a2)+24.0*a1*a1*a2*log(a1/a))/(3.0*std::pow(a12,5.0));
		}
		// expected N(L) for alpha
		// P&Y eqn. 30 integrated
		inline clipper::ftype meanNL(const clipper::ftype a, const clipper::ftype L) const {
			clipper::ftype a12 = 1.0-2.0*a;
			clipper::ftype a1 = 1.0-a;
			clipper::ftype a2 = a*a;
			return (2.0/(a12*a12))*(0.5*L*(a2+a1*a1) + a*a1*a1*( -1.0/(1.0+a12*L) + 1.0/(a12*(1.0+a12*L)) ) - a*a*a1*( 1.0/(1.0-a12*L) + 1.0/(a12*(1.0-a12*L)) ) );
		}
		
    };
        
    //! Compute H-test on intensities
    /*! Compute the H-test using tabulated or first principles twinning
     operators.
     Yeates Acta Cryst. A44 142 (1980)
     */
    class H_test {
    public:
        // constructor, setup class
        H_test(int nbins=20) { _nbins=nbins; return; };
        // constructor, calculates H tests if true
        template<class T> explicit H_test( const clipper::HKL_data< clipper::datatypes::I_sigI<T> >& isig, const clipper::Isymop& twinop, clipper::Resolution reso = clipper::Resolution(3.0), int nbins=20);
        // perform calculation, or recalculate value
        template<class T> const clipper::ftype operator()(const clipper::HKL_data< clipper::datatypes::I_sigI<T> >& isig, const clipper::Isymop& twinop, clipper::Resolution reso = clipper::Resolution(3.0)) ;
        // return Hav from H-test
        const clipper::ftype statistics() const { return _Hav; }   
        // return twin fraction
        const clipper::ftype fraction() const { return _alpha; }
        // return tested symop
        const clipper::Isymop& twinops() const { return _twinop; };
        // return operator string
        const std::string description() const;
        // print summary
        void summary() const ;
        // print loggraph
        void loggraph() const ;
        
    private:
        clipper::Resolution _reso; //!< resolution used in calculation
        int _nbins; //!< number of bins in cdfs and pdfs
        clipper::Isymop _twinop; //!< twinning operator
        std::vector<int> _cdf; //!< cumulative density functions
        clipper::ftype _alpha; //!< alpha values from H-test 
        clipper::ftype _Hav;   //!< average value from H-test
        clipper::ftype _H2av;  //!< H statistic squared
        
        template<class T> void calc(const clipper::HKL_data< clipper::datatypes::I_sigI<T> >& isig, const clipper::Isymop& twinop, const clipper::Resolution& reso);
    };
    
    //! Compute Britton-test on intensities
    /*! Compute the Britton-test using tabulated or first principles twinning
     operators.
     Britton Acta Cryst. A28 296 (1972)
     */
    class Britton_test {
    public:
        // constructor, setup class
        Britton_test(int nbins=20) { _nbins=nbins; return; };
        // constructor, calculates Britton tests using fp if true
        template<class T> explicit Britton_test( const clipper::HKL_data< clipper::datatypes::I_sigI<T> >& isig, const clipper::Isymop& twinop, clipper::Resolution reso = clipper::Resolution(3.0), int nbins=20);
        // perform calculation, or recalculate value
        template<class T> const clipper::ftype operator()(const clipper::HKL_data< clipper::datatypes::I_sigI<T> >& isig, const clipper::Isymop& twinop,clipper::Resolution reso = clipper::Resolution(3.0)) ;
        // return score
        const clipper::ftype statistic() const { return _alpha; }
        // return twin fraction
        const clipper::ftype fraction() const { return _alpha; }
        // return tested symop
        const clipper::Isymop& twinop() const { return _twinop; }
        // return text description of twinop
        const std::string description() const;
        // print summary
        void summary() const ;
        // print loggraph
        void loggraph() const ;
        
    private:
        clipper::Isymop _twinop; //!< twinning operator
        int _nbins; //!< number of bins for distributions
        clipper::Resolution _reso; //!< resolution used in calculation
        std::vector <clipper::ftype> _pdf; //!< Britton plot
        std::vector <clipper::ftype> _zpdf; //!< zero probability density functions
        clipper::ftype _NTL; //!< number of data points
        clipper::ftype _alpha; //!< alpha values from Britton-tests
        clipper::ftype _beta; //!< offsets for Murray-Rust plot
        
        template<class T> void calc(const clipper::HKL_data< clipper::datatypes::I_sigI<T> >& isig, const clipper::Isymop& twinop, const clipper::Resolution& reso);
    };


    //! Compute Britton-test on intensities, including covariance
    /*! Compute the Britton-test using tabulated or first principles twinning
     operators.
     Britton Acta Cryst. A28 296 (1972)
     */
    class MLBritton_test  {
    public:
        // constructor, setup class
        MLBritton_test(int nbins=20) {  _nbins=nbins; return; };
        // constructor, calculates Britton tests using fp if true
        template<class T> explicit MLBritton_test( const clipper::HKL_data< clipper::datatypes::I_sigI<T> >& isig, const clipper::Isymop& twinop,  clipper::Coord_frac ncs=clipper::Coord_frac(0.0,0.0,0.0), clipper::Resolution reso = clipper::Resolution(3.0), int nbins=20);
        // perform calculation, or recalculate value
        template<class T> const clipper::ftype operator()( const clipper::HKL_data< clipper::datatypes::I_sigI<T> >& isig, const clipper::Isymop& twinop,  clipper::Coord_frac ncs=clipper::Coord_frac(0.0,0.0,0.0), clipper::Resolution reso = clipper::Resolution(3.0)) ;
        // return score
        const clipper::ftype statistics() const { return _llk; }  
        // return twin fraction
        const clipper::ftype fraction() const { return _alpha; } 
        // return error in coordinates
        const clipper::ftype deltaR() const { return _dr; }
        // return tested symops
        const clipper::Isymop& twinop() const { return _twinop; }
        // return operator string
        const std::string description() const;
        // print summary
        void summary() const ;
        // print loggraph
        void loggraph() const ;
        
    private:
        clipper::Isymop _twinop; //!< twinning operator
        clipper::Coord_frac _ncs; //!< ncs vector
        clipper::Resolution _reso; //!< resolution used in calculation
        int _nbins; //!< number of bins for distributions
        std::vector<clipper::ftype> _pdf; //!< zero probability density functions
        clipper::ftype _alpha; //!< alpha values from Britton-tests
        clipper::ftype _llk; //!< score
        clipper::ftype _dr; //!< delta R
        
        
        template<class T> void calc(const clipper::HKL_data< clipper::datatypes::I_sigI<T> >& isig, const clipper::Isymop& twinop, const clipper::Resolution& reso, const clipper::ResolutionFn& norm, const clipper::Coord_frac& ncs);
        
        template <class T> const clipper::ftype point(const clipper::HKL_data< clipper::datatypes::I_sigI<T> >& isig, const clipper::ResolutionFn& norm, const clipper::Isymop& twinop, const clipper::ftype& alpha, const clipper::ftype& Dr);
    };

}

#endif

