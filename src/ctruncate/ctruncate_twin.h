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
#include "ctruncate_moments.h"
#include "ctruncate_analyse.h"
#include "intensity_target.h"
#include <string>

namespace ctruncate {
    
    
    void twin_summary(clipper::ftype, clipper::ftype);
    
    //! compute the twinning operators
    /*! Calculate either by first principes or
     via tables
     */
    /*class TwinSymops {
     public:
     enum MODE {FP,TABLE}; //!< use first principles or tabulated twinning operators
     // constructor
     TwinSymops(MODE mode=FP ) : _mode(mode) {};
     explicit TwinSymops(const clipper::Cell& cell, const clipper::Spacegroup& spg, MODE mode=FP );
     // return ith twinop as ISymop
     int operator()(const clipper::Cell& cell, const clipper::Spacegroup& spg );
     const clipper::Isymop& operator[](int i) const { return _twinops[i]; }
     // return number of twinning operators
     int size() const { return _twinops.size(); }
     // return description of i'th symop
     const std::string description(const int&) const;
     // return summary
     void summary() const;
     
     
     private:
     // set up tabulated symops
     void table(const clipper::Cell&, const clipper::Spacegroup&);
     // set up first principle symops
     void principles(const clipper::Cell&, const clipper::Spacegroup&);
     
     MODE _mode; //record mode
     std::vector<clipper::Isymop> _twinops; //!< twinning operators
     }; */
    
    class CumulativeAnalysis {
    public:
        // default constructore
        inline CumulativeAnalysis() : _ntw(0) {}
        // constructure with calc
        template <class T, template<class> class D> CumulativeAnalysis(const clipper::HKL_data<D<T> >&)  ;
        // operator()
        template <class T, template<class> class D> bool operator()(const clipper::HKL_data<D<T> >&);
        // output to std::out
        void output() const;
        // output in xml
        std::stringstream& xml_output(std::stringstream&) const;
        // number of possible twinned
        int numbertwinned() { return _ntw; }
        
    private:
        template<class T> const T&    obs( const clipper::datatypes::F_sigF<T>& f ) { return f.f()*f.f(); }
        template<class T> const T&    obs( const clipper::datatypes::I_sigI<T>& f ) { return f.I(); }
        template<class T> const T&    obs( const clipper::datatypes::F_sigF_ano<T>& f ) { return f.f()*f.f(); }
        template<class T> const T&    obs( const clipper::datatypes::I_sigI_ano<T>& f ) { return f.I(); }
        
        static double _acentricideal[51];  //!< ideal data for acentric reflections, also twinned centric
        static double _centricideal[51];   //!< ideal data for centric reflections
        static double _acentrictwin[51];   //!< ideal data for twinned acentric reflections
        
        int _ntw;
        
        clipper::Generic_ordinal _intensity_ord_c;  //!< store centric cumulative
        clipper::Generic_ordinal _intensity_ord_a;  //!< store acentric cumulative
    };
    
    //TwinSymops-----------------------------------------------------
    
    //! compute the twinning operators
    /*! Calculate either by first principes
     */
    class TwinSymops {
    public:
        class Twinop : public clipper::Isymop {
        public:
            inline Twinop() {}
            inline explicit Twinop( const clipper::RTop<int>& rt, const clipper::ftype sc ) : clipper::Isymop(rt), _sc(sc) {}
            Twinop( const clipper::Symop& symop, const clipper::Grid& grid, const clipper::ftype sc ) : clipper::Isymop(symop,grid), _sc(sc) {}
            std::string description() const;
            clipper::ftype score() const {
                return _sc;
            }
            std::string type(const clipper::Spacegroup& spgr) const;
        private:
            clipper::ftype _sc; //!< score
            //std::string _type;  //!< type of operator
        };
        
        // constructor
        inline TwinSymops() {}
        TwinSymops(TwinSymops& parent){
            _twinops = parent._twinops;
            _cell = parent._cell;
            _spg = parent._spg;
        }
        template<class T, template<class> class D> TwinSymops(clipper::HKL_data<D<T> >&);
        TwinSymops(const clipper::Cell&, const clipper::Spacegroup&);
        // operator() return number of sumops
        template<class T, template<class> class D> int operator()(clipper::HKL_data<D<T> >&);
        int operator()(const clipper::Cell&, const clipper::Spacegroup&);
        // destructor
        ~TwinSymops() {}
        // return number of twinning operators
        int size() const { return _twinops.size(); }
        // return i'th twinning operator as symop
        const Twinop& operator[](const int& i) const { return _twinops[i]; }
        // return description of i'th symop
        std::string description(const int&) const;
        // return summary
        void output() const;
        //return xml
        std::stringstream& xml_output(std::stringstream&);
        
    private:
        // set up first principle symops
        void principles(const clipper::Cell&, const clipper::Spacegroup&);
        
        clipper::Cell _cell;  //!< cell
        clipper::Spacegroup _spg; //!< spacegroup
        std::vector<Twinop> _twinops; //!< twinning operators
    };
    
    //L_test----------------------------------------------------------
    //! Compute L-test scores on intensity
    /*! Compute the L-test, store results for summary and loggraph.
     Padilla and Yeates Acta Cryst. D59 1124 (2003)
     */
    class L_test {
    public:
        // constructor, setup up L-test
        L_test( int nbins=20) : _cdf(nbins,0.0) {}
        // constructor, calculates L statistics if true
        template<class T, template<class> class D> explicit L_test( clipper::HKL_data< D<T> >&, clipper::Range<clipper::ftype> reso = clipper::Range<clipper::ftype>(), int nbins=20);
        template<class T, template<class> class D> explicit L_test( clipper::HKL_data< D<T> >&, std::vector< clipper::Symop>& tncs, clipper::Range<clipper::ftype> reso = clipper::Range<clipper::ftype>(), int nbins=20);
        // perform calculation, or recalculate value
        template<class T, template<class> class D> clipper::ftype operator()( clipper::HKL_data< D<T> >&, clipper::Range<clipper::ftype> reso = clipper::Range<clipper::ftype>() ) ;
        template<class T, template<class> class D> clipper::ftype operator()( clipper::HKL_data< D<T> >&, std::vector< clipper::Symop>& tncs, clipper::Range<clipper::ftype> reso = clipper::Range<clipper::ftype>() ) ;
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
        clipper::ftype fraction() const { return _alpha; }
        // print summary
        void summary() const ;
        // print loggraph
        void loggraph() const ;
        std::stringstream& xml_output(std::stringstream&);
        std::stringstream& xml_graph(std::stringstream&) const;
        
    private:
        clipper::Range<clipper::ftype> _reso; //!< resolution used in calculation
        clipper::ftype _Lav; //!< L-statistic
        clipper::ftype _L2av; //!< L-statistic squared
        clipper::ftype _alpha; //!< alpha values from H-test
        clipper::ftype _NLT; //!< Number of data points
        std::vector<clipper::ftype> _cdf; //!< cumulative density function
        
        // perform calculation
        template<class T, template<class> class D> clipper::ftype calc(clipper::HKL_data< D<T> >& isig, std::vector<clipper::Symop>& tncs );
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
    
    //H_test---------------------------------------------------------
    //! Compute H-test on intensities
    /*! Compute the H-test using tabulated or first principles twinning
     operators.
     Yeates Acta Cryst. A44 142 (1980)
     */
    class H_test {
    public:
        // default constructor, setup class
        H_test(int nbins=20) : _cdf(nbins,0.0) { }
        // constructor, calculates H tests if true
        template<class T, template<class> class D> explicit H_test( const clipper::HKL_data< D<T> >& isig, const clipper::Isymop& twinop, clipper::Range<clipper::ftype> reso = clipper::Range<clipper::ftype>(), int nbins=20);
        // perform calculation, or recalculate value
        template<class T, template<class> class D> clipper::ftype operator()(const clipper::HKL_data< D<T> >& isig, const clipper::Isymop& twinop, clipper::Range<clipper::ftype> reso = clipper::Range<clipper::ftype>()) ;
        /*H_test& operator=(const H_test& htest) {
            if ( this != &htest) {
                _reso = htest._reso;
                _twinop = htest._twinop;
                _cdf.assign(htest._cdf.begin(),htest._cdf.end() );
                _alpha = htest._alpha;
                _Hav = htest._Hav;
                _H2av = htest._H2av;
            }
            return *(this);
        }*/
        // return Hav from H-test
        clipper::ftype statistics() const { return _Hav; }
        // return twin fraction
        clipper::ftype fraction() const { return _alpha; }
        // return tested symop
        const clipper::Isymop& twinops() const { return _twinop; }
        // return operator string
        const std::string description() const;
        // print summary
        void summary() const ;
        // print loggraph
        void loggraph() const ;
        // return array of cdf
        const std::vector<clipper::ftype> array() { return _cdf; }
        
    private:
        clipper::Range<clipper::ftype> _reso; //!< resolution used in calculation
        clipper::Isymop _twinop; //!< twinning operator
        std::vector<clipper::ftype> _cdf; //!< cumulative density functions
        clipper::ftype _alpha; //!< alpha values from H-test
        clipper::ftype _Hav;   //!< average value from H-test
        clipper::ftype _H2av;  //!< H statistic squared
        
        template<class T, template<class> class D> void calc(const clipper::HKL_data< D<T> >& isig, const clipper::Isymop& twinop);
    };
    
    //Britton-------------------------------------------------------
    
    //! Compute Britton-test on intensities
    /*! Compute the Britton-test using tabulated or first principles twinning
     operators.
     Britton Acta Cryst. A28 296 (1972)
     */
    class Britton_test {
    public:
        // constructor, setup class
        Britton_test(int nbins=20) : _pdf(nbins,0.0), _zpdf(nbins,0.0) { }
        // constructor, calculates Britton tests using fp if true
        template<class T, template<class> class D> explicit Britton_test( const clipper::HKL_data< D<T> >& isig, const clipper::Isymop& twinop, clipper::Range<clipper::ftype> reso = clipper::Range<clipper::ftype>(), int nbins=20);
        // perform calculation, or recalculate value
        template<class T, template<class> class D> clipper::ftype operator()(const clipper::HKL_data< D<T> >& isig, const clipper::Isymop& twinop,clipper::Range<clipper::ftype> reso = clipper::Range<clipper::ftype>()) ;
        // return score
        clipper::ftype statistic() const { return _alpha; }
        // return twin fraction
        clipper::ftype fraction() const { return _alpha; }
        // return tested symop
        const clipper::Isymop& twinop() const { return _twinop; }
        // return text description of twinop
        const std::string description() const;
        // print summary
        void summary() const ;
        // print loggraph
        void loggraph() const ;
        // return array of zpdf
        const std::vector<clipper::ftype> array() { return _zpdf; }
        
    private:
        clipper::Isymop _twinop; //!< twinning operato
        clipper::Range<clipper::ftype> _reso; //!< resolution used in calculation
        std::vector <clipper::ftype> _pdf; //!< Britton plot
        std::vector <clipper::ftype> _zpdf; //!< zero probability density functions
        clipper::ftype _NTL; //!< number of data points
        clipper::ftype _alpha; //!< alpha values from Britton-tests
        clipper::ftype _beta; //!< offsets for Murray-Rust plot
        
        template<class T, template<class> class D> void calc(const clipper::HKL_data< D<T> >& isig, const clipper::Isymop& twinop);
    };
    
    //ML-Britton-----------------------------------------------------
    //! Compute Britton-test on intensities, including covariance
    /*! Compute the Britton-test using tabulated or first principles twinning
     operators.
     Britton Acta Cryst. A28 296 (1972)
     */
    class MLBritton_test  {
    public:
        // constructor, setup class
        MLBritton_test(int nbins=20) : _pdf(nbins,0.0) { }
        // constructor, calculates Britton tests using fp if true
        template<class T, template<class> class D> explicit MLBritton_test( const clipper::HKL_data< D<T> >& isig, const clipper::Isymop& twinop,  clipper::Range<clipper::ftype> reso = clipper::Range<clipper::ftype>(),  clipper::Symop ncs=clipper::Symop(clipper::Symop::identity() ), int nbins=20);
        //
        //template<class T, template<class> class D> explicit MLBritton_test( const clipper::HKL_data< D<T> >& isig, clipper::Symop ncs=clipper::Symop(clipper::Symop::identity() ), clipper::Range<clipper::ftype> reso = clipper::Range<clipper::ftype>(), int nbins=20) : _pdf(nbins,0.0);
        // perform calculation, or recalculate valu
        template<class T, template<class> class D> clipper::ftype operator()( const clipper::HKL_data< D<T> >& isig, const clipper::Isymop& twinop,
                                                                             clipper::Range<clipper::ftype> reso = clipper::Range<clipper::ftype>(),
                                                                             clipper::Symop ncs=clipper::Symop(clipper::Symop::identity() ) );
        // return score
        clipper::ftype statistics() const { return _llk; }
        // return twin fraction
        clipper::ftype fraction() const { return _alpha; }
        // return error in coordinates
        clipper::ftype deltaR() const { return _dr; }
        // return tested symops
        const clipper::Isymop& twinop() const { return _twinop; }
        // return operator string
        const std::string description() const;
        // print summary
        void summary() const ;
        // print loggraph
        void loggraph() const ;
        // return array of zpdf
        const std::vector<clipper::ftype> array() { return _pdf; }
        
    private:
        clipper::Isymop _twinop; //!< twinning operator
        clipper::Symop _ncs; //!< ncs vector
        clipper::Range<clipper::ftype> _reso; //!< resolution used in calculation
        std::vector<clipper::ftype> _pdf; //!< zero probability density functions
        clipper::ftype _alpha; //!< alpha values from Britton-tests
        clipper::ftype _llk; //!< score
        clipper::ftype _dr; //!< delta R
        
        
        template<class T, template<class> class D> void calc(const clipper::HKL_data< D<T> >& isig, const clipper::Isymop& twinop, const clipper::ResolutionFn& norm, const clipper::Symop& ncs);
        
        template<class T, template<class> class D> clipper::ftype point(const clipper::HKL_data< D<T> >& isig, const clipper::ResolutionFn& norm, const clipper::Isymop& twinop, const clipper::ftype& alpha, const clipper::ftype& Dr);
    };
    
    //R_test---------------------------------------------------------
    //! Compute twinning R on intensities
    /*! Compute the R values using tabulated or first principles twinning
     operators.
     */
    class R_test : public ResolStats_base {
    public:
        // default constructor, setup class
        R_test() { };
        // constructor, calculates R if true
        template<class T, template<class> class D> explicit R_test( const clipper::HKL_data< D<T> >& isig, const clipper::Isymop& twinop, clipper::Range<clipper::ftype> reso = clipper::Range<clipper::ftype>() );
        // perform calculation, or recalculate value
        template<class T, template<class> class D> clipper::ftype operator()(const clipper::HKL_data< D<T> >& isig, const clipper::Isymop& twinop, clipper::Range<clipper::ftype> reso = clipper::Range<clipper::ftype>()) ;
        // return value in shell i
        clipper::ftype operator[](int index) const { return _r[index]; }
        // return Rav from R-test
        clipper::ftype statistics() const { return _Rav; }
        // return tested symop
        const clipper::Isymop& twinops() const { return _twinop; };
        // return operator string
        const std::string description() const;
        // print summary
        void summary() const ;
        // print loggraph
        void loggraph() const ;
        // return array of cdf
        const std::vector<clipper::ftype> array() { return _r; }
        
    private:
        clipper::Range<clipper::ftype> _reso; //!< resolution used in calculation
        clipper::Isymop _twinop; //!< twinning operator
        std::vector<clipper::ftype> _r; //!< cumulative density functions
        clipper::ftype _Rav;   //!< average value from R-test
    };
    
    //TwinAnalysis-----------------------------------------------------
    //! Wrap twin tests
    /*! Wrap the twin tests and output
     */
    class TwinAnalysis {
    public:
        //construct with data
        template<class T, template<class> class D >  TwinAnalysis(clipper::HKL_data< D<T> >&);
        template<class T, template<class> class D >  TwinAnalysis(clipper::HKL_data< D<T> >&, std::vector<clipper::Symop>&);
        template<class T, template<class> class D, class TT >  TwinAnalysis(clipper::HKL_data< D<T> >&, std::vector<clipper::Symop>&, clipper::Range<TT>&);
        //template<class D, class TT >  TwinAnalysis(clipper::HKL_data< D >&, std::vector<clipper::Symop>&, clipper::Range<TT>&);
        //print to std::out
        void output();
        //return xml
        std::stringstream& xml_output(std::stringstream& );
    private:
        template<class T, template<class> class D, class TT> void init(clipper::HKL_data< D<T> >&, std::vector<clipper::Symop>&, clipper::Range<TT>&);
        
        clipper::HKL_data_base *_base;              //!< data
        clipper::Range<clipper::ftype> _range;     //!< active resolution range
        std::vector<clipper::Symop> _ncsops;                     //!< effective symops
        TwinSymops _ts;                                          //!< possible twin operators
        CumulativeAnalysis _ca;                     //!< cumulative analysis
        L_test _ltest;                              //!< L-test
        Moments _moments;                          //!< moments calc object
        std::vector<MLBritton_test> _mlbtests;     //!< ML Britton object
        std::vector<Britton_test> _btests;         //!< Britton object
        std::vector<H_test> _htests;               //!< H-test object
        //std::vector<std::vector<R_test> > _rtests;               //!< |R| objects
    };
    
    //definitions-----------------------------------------------------
    
    template<class T, template<class> class D> TwinAnalysis::TwinAnalysis(clipper::HKL_data<D<T> >& hkldata)
    {
        std::vector<clipper::Symop> ncsops;
        clipper::Range<clipper::ftype> range;
        init(hkldata,ncsops,range);
    }
    
    template<class T, template<class> class D> TwinAnalysis::TwinAnalysis(clipper::HKL_data<D<T> >& hkldata, std::vector<clipper::Symop>& ncsops)
    {
        clipper::Range<clipper::ftype> range;
        init(hkldata,ncsops,range);
    }
    
    template<class T, template<class> class D, class TT> TwinAnalysis::TwinAnalysis(clipper::HKL_data<D<T> >& hkldata, std::vector<clipper::Symop>& ncsops, clipper::Range<TT>& range)
    {
        init(hkldata,ncsops,range);
    }
    
    template<class T, template<class> class D, class TT> void TwinAnalysis::init(clipper::HKL_data<D<T> >& hkldata, std::vector<clipper::Symop>& ncsops, clipper::Range<TT>& range)
    {
        _base = &hkldata;
        _range = ((range.min() > range.max() ) ? hkldata.invresolsq_range() : range );
        _ncsops = ncsops;
        
        _ca(hkldata);
        
        _ts(hkldata );
        
        _ltest = (L_test(hkldata,_ncsops,range) );
        
        _moments(hkldata,range);
        
        _htests.reserve(_ts.size() );
        //_htests.resize(_ts.size(),H_test() );
        for (int i = 0; i != _ts.size() ; ++i ) {
            _htests.push_back(H_test(hkldata,_ts[i],range) );
            //(_htests[i])(hkldata,_ts[i],range);
        }

        _btests.reserve(_ts.size() );
        for (int i = 0; i != _ts.size() ; ++i ) {
            _btests.push_back(Britton_test(hkldata,_ts[i],range) );
        }
        
        _mlbtests.reserve(_ts.size() );
        for (int i = 0; i != _ts.size() ; ++i ) {
            clipper::ftype product(0.0);
            int jp;
            if ( _ncsops.size() ) {
                for (int j=0; j != _ncsops.size() ; ++j) {
                    clipper::Vec3<clipper::ftype> vect = _ncsops[j].rtop_orth(hkldata.cell() ).trn();
                    clipper::Mat33<int> tmp = _ts[i].rot();
                    clipper::Rotation rot(clipper::Mat33<clipper::ftype>(
                                                                         clipper::ftype(tmp(0,0) )/12.0,
                                                                         clipper::ftype(tmp(0,1) )/12.0,
                                                                         clipper::ftype(tmp(0,2) )/12.0,
                                                                         clipper::ftype(tmp(1,0) )/12.0,
                                                                         clipper::ftype(tmp(1,1) )/12.0,
                                                                         clipper::ftype(tmp(1,2) )/12.0,
                                                                         clipper::ftype(tmp(2,0) )/12.0,
                                                                         clipper::ftype(tmp(2,1) )/12.0,
                                                                         clipper::ftype(tmp(2,2) )/12.0));
                    clipper::Vec3<clipper::ftype> vecr(rot.x(),rot.y(),rot.z() );
                    clipper::ftype p = vecr.unit()*vect.unit();
                    if (p > product) {
                        product = p;
                        jp = j;
                    }
                }
            }
            if (product > 0.95 ) _mlbtests.push_back(MLBritton_test(hkldata,_ts[i],range,_ncsops[jp]) );
            else {
                _mlbtests.push_back(MLBritton_test(hkldata,_ts[i],range) );
            }
        }
    }
    
    
    /*! constructor for list of twinning operators.
     This uses either first principles or tabulated values.
     \param HDL_data type
     \return type TwinSymops */
    template<class T, template<class> class D> TwinSymops::TwinSymops(clipper::HKL_data<D<T> >& hkldata)
    {
        _cell = hkldata.hkl_info().cell();
        _spg = hkldata.hkl_info().spacegroup();
        principles(_cell,_spg );
    }
    
    /*! operator() to perform calculation
     \param HKL_data object
     \return in number of symops */
    template<class T, template<class> class D> int TwinSymops::operator()(clipper::HKL_data<D<T> >& hkldata)
    {
        _cell = hkldata.hkl_info().cell();
        _spg = hkldata.hkl_info().spacegroup();
        principles(_cell,_spg );
        return size();
    }
    
    /*! constructor for the L-test.
     \param isig Experimental intensities
     \param reso Resolution limit for calculation (default 3.0A)
     \param nbins Number of bins for the CDF (default 20)
     \return type L_test */
    template<class T, template<class> class D> L_test::L_test( clipper::HKL_data< D<T> >& isig, clipper::Range<clipper::ftype> reso, int nbins) {
        _reso = reso;
        if (_reso.max() < 0 ) _reso = clipper::Range<clipper::ftype>(isig.hkl_info().invresolsq_range() );
        _alpha = -1.0;
        //_cdf.resize(nbins);
        std::vector<clipper::Symop> tncs;
        calc(isig,tncs);
        return;
    }
    
    /*! constructor for the L-test.
     \param isig Experimental intensities
     \param tncs Information on existing ncs operators
     \param reso Resolution limit for calculation (default 3.0A)
     \param nbins Number of bins for the CDF (default 20)
     \return type L_test */
    template<class T, template<class> class D> L_test::L_test( clipper::HKL_data< D<T> >& isig, std::vector<clipper::Symop>& tncs, clipper::Range<clipper::ftype> reso, int nbins) {
        _reso = reso;
        if (_reso.max() < 0 ) _reso = clipper::Range<clipper::ftype>(isig.hkl_info().invresolsq_range() );
        _cdf.resize(nbins,0.0);
        _alpha = -1.0;
        
        calc(isig,tncs);
        return;
    }
    
    /*! recalculate L-test.
     \param reso Resolution limit for calculation (default 3.0A)
     \return L-statistic */
    template<class T, template<class> class D> clipper::ftype L_test::operator() (clipper::HKL_data< D<T> >& isig, clipper::Range<clipper::ftype> reso ) {
        _reso = reso;
        if (_reso.max() < 0 ) _reso = clipper::Range<clipper::ftype>(isig.hkl_info().invresolsq_range() );
        _alpha = -1.0;
        std::vector<clipper::Symop> tncs;
        
        return calc(isig,tncs);
    }
    
    template<class T, template<class> class D> clipper::ftype L_test::operator() (clipper::HKL_data< D<T> >& isig, std::vector<clipper::Symop>& tncs, clipper::Range<clipper::ftype> reso ) {
        _reso = reso;
        if (_reso.max() < 0 ) _reso = clipper::Range<clipper::ftype>(isig.hkl_info().invresolsq_range() );
        _alpha = -1.0;
        
        return calc(isig,tncs);
    }
    
    
    /*! constructor for the H-test.
     \param isig Experimental intensities
     \param twinop Twinning operator to investigate
     \param reso Resolution limit for calculation (default 3.0A)
     \param nbins Number of bins for the CDF (default 20)
     \return type H_test */
    template<class T, template<class> class D> H_test::H_test( const clipper::HKL_data< D<T> >& isig, const clipper::Isymop& twin, clipper::Range<clipper::ftype> reso , int nbins) {
        _cdf.resize(nbins,0.0);
        _reso = reso;
        if (_reso.max() < 0 ) _reso = clipper::Range<clipper::ftype>(isig.hkl_info().invresolsq_range() );
        _twinop = twin;
        calc(isig,twin);
    }
    
    /*! recalculate H-test.
     \param reso Resolution limit for calculation (default 3.0A)
     \return alphas */
    template<class T, template<class> class D> clipper::ftype H_test::operator() (const clipper::HKL_data< D<T> >& isig, const clipper::Isymop& twin, clipper::Range<clipper::ftype> reso ) {
        _reso=reso;
        if (_reso.max() < 0 ) _reso = clipper::Range<clipper::ftype>(isig.hkl_info().invresolsq_range() );
        _twinop=twin;
        calc(isig,twin);
        return _alpha;
    }
    
    /*! constructor for the Britton-test.
     \param isig Experimental intensities
     \param twinop Twinning operator
     \param reso Resolution limit for calculation (default 3.0A)
     \param nbins Number of bins for the CDF (default 50)
     \return type Britton_test */
    template<class T, template<class> class D> Britton_test::Britton_test( const clipper::HKL_data< D<T> >& isig, const clipper::Isymop& twinop, clipper::Range<clipper::ftype> reso , int nbins) {
        _reso=reso;
        if (_reso.max() < 0 ) _reso = clipper::Range<clipper::ftype>(isig.hkl_info().invresolsq_range() );
        _twinop=twinop;
        _pdf.resize(nbins,0.0);
        _zpdf.resize(nbins,0.0);
        calc(isig,twinop);
    }
    
    /*! recalculate Britton-test.
     \param reso Resolution limit for calculation (default 3.0A)
     \return alphas */
    template<class T, template<class> class D> clipper::ftype Britton_test::operator() (const clipper::HKL_data< D<T> >& isig, const clipper::Isymop& twinop, clipper::Range<clipper::ftype> reso ) {
        _reso=reso;
        if (_reso.max() < 0 ) _reso = clipper::Range<clipper::ftype>(isig.hkl_info().invresolsq_range() );
        _twinop=twinop;
        calc(isig,twinop);
        return _alpha;
    }
    
    /*! constructor for the ML Britton test.
     \param hkldata Experimental intensities
     \param twinop Twinning operator to be tested
     \param ncs tNCS vector as a Coord_frac (default NULL)
     \param reso Resolution limit for calculation (default 3.0A)
     \param nbins Number of bins for the CDF (default 20)
     \return type MLBritton_test */
    template<class T, template<class> class D> MLBritton_test::MLBritton_test( const clipper::HKL_data< D<T> >& hkldata, const clipper::Isymop& twinop, clipper::Range<clipper::ftype> reso, clipper::Symop ncs, int nbins ) {
        _reso=reso;
        if (_reso.max() < 0 ) _reso = clipper::Range<clipper::ftype>(hkldata.hkl_info().invresolsq_range() );
        _twinop = twinop;
        _ncs = ncs;
        _pdf.resize(nbins,0.0);
        //norm using lots of bins
        int nprm(0);
        int Nreflections(hkldata.num_obs() );
        const int nreflns(500);
        {
            //for ( clipper::HKL_data_base::HKL_reference_index ih = hkldata.first(); !ih.last(); ih.next() ) {
            //		 if (!hkldata.missing(ih.index() ) )++Nreflections;
            //}
            if ( nprm == 0 && nreflns != 0 ) {
                nprm = std::max( Nreflections/nreflns , 1);
                //} else if ( nreflns == 0 && nprm2 != 0 ) {
                //nprm = nprm;
            } else {
                //nprm2 = std::max( Nreflections/nreflns , nprm2);
                double np1(nprm+0.499);
                double np2(Nreflections/nreflns);
                double np(std::sqrt(np1*np1*np2*np2/(np1*np1+np2*np2) ) );
                nprm = std::max( int(np), 1 );
            }
        }
        std::vector<double> params( nprm, 1.0 );
        clipper::BasisFn_spline basis( hkldata, nprm, 2.0 );
        //TargetFn_meanInth<clipper::data32::I_sigI> target( hkldata, 1 );
        TargetFn_meanInth<D <T> > target( hkldata, 1 );
        clipper::ResolutionFn norm( hkldata.hkl_info(), basis, target, params );
        
        calc(hkldata,twinop,norm,ncs);
    }
    
    /*! calculate Britton-test.
     \param hkldata Intensity data
     \param twinop Twinning operator as Isymop
     \param ncs tNCS vector as Coord_frac (default NULL)
     \param reso Resolution limit for calculation (default 3.0A)
     \return alpha */
    template<class T, template<class> class D> clipper::ftype MLBritton_test::operator()( const clipper::HKL_data< D<T> >& hkldata, const clipper::Isymop& twinop,  clipper::Range<clipper::ftype> reso, clipper::Symop ncs)
    {
        _reso=reso;
        if (_reso.max() < 0 ) _reso = clipper::Range<clipper::ftype>(hkldata.hkl_info().invresolsq_range() );
        _twinop = twinop;
        _ncs = ncs;
        
        int nprm(0);
        int Nreflections(hkldata.num_obs() );
        const int nreflns(500);
        {
            //for ( clipper::HKL_data_base::HKL_reference_index ih = hkldata.first(); !ih.last(); ih.next() ) {
            //		 if (!hkldata.missing(ih.index() ) )++Nreflections;
            //}
            if ( nprm == 0 && nreflns != 0 ) {
                nprm = std::max( Nreflections/nreflns , 1);
                //} else if ( nreflns == 0 && nprm2 != 0 ) {
                //nprm = nprm;
            } else {
                //nprm2 = std::max( Nreflections/nreflns , nprm2);
                double np1(nprm+0.499);
                double np2(Nreflections/nreflns);
                double np(std::sqrt(np1*np1*np2*np2/(np1*np1+np2*np2) ) );
                nprm = std::max( int(np), 1 );
            }
        }
        
        std::vector<double> params( nprm, 1.0 );
        clipper::BasisFn_spline basis( hkldata, nprm, 2.0 );
        //TargetFn_meanInth<clipper::datatypes::I_sigI<T> > target( hkldata, 1 );
        TargetFn_meanInth<D <T> > target( hkldata, 1 );
        clipper::ResolutionFn norm( hkldata.hkl_info(), basis, target, params );
        
        calc(hkldata,twinop,norm,ncs);
        return _alpha;
    }
    
    /*! constructor for the Cumulative plots.
     \param hkldata Experimental intensities
     \return type Cumulative */
    template<class T, template<class> class D> CumulativeAnalysis::CumulativeAnalysis( const clipper::HKL_data< D<T> >& hkldata ) : _ntw(0) {
        this->operator()(hkldata);
    }
    
    /*! calculate Cumlative distribution.
     \param hkldata Intensity data
     \return bool Indicator of twinning */
    template<class T, template<class> class D> bool CumulativeAnalysis::operator() (const clipper::HKL_data< D<T> >& hkldata) {
        typedef clipper::HKL_data_base::HKL_reference_index HRI;
        
        int nprm(0);
        int Nreflections(hkldata.num_obs() );
        const int nreflns(500);
        {
            //for ( clipper::HKL_data_base::HKL_reference_index ih = hkldata.first(); !ih.last(); ih.next() ) {
            //		 if (!hkldata.missing(ih.index() ) )++Nreflections;
            //}
            if ( nprm == 0 && nreflns != 0 ) {
                nprm = std::max( Nreflections/nreflns , 1);
                //} else if ( nreflns == 0 && nprm2 != 0 ) {
                //nprm = nprm;
            } else {
                //nprm2 = std::max( Nreflections/nreflns , nprm2);
                double np1(nprm+0.499);
                double np2(Nreflections/nreflns);
                double np(std::sqrt(np1*np1*np2*np2/(np1*np1+np2*np2) ) );
                nprm = std::max( int(np), 1 );
            }
        }
        
        std::vector<double> params( nprm, 1.0 );
        clipper::BasisFn_spline basis( hkldata, nprm, 2.0 );
        //TargetFn_meanInth<clipper::datatypes::I_sigI<T> > target( hkldata, 1 );
        TargetFn_meanInth<D<T> > target( hkldata, 1 );
        clipper::ResolutionFn norm( hkldata.hkl_info(), basis, target, params );
        
        // construct cumulative distribution function for intensity (using Z rather than E)
        clipper::Range<double> intensity_range_centric;
        clipper::Range<double> intensity_range_acentric;
        for ( HRI ih = hkldata.first(); !ih.last(); ih.next() ) {
            if ( !hkldata[ih].missing() ) {
                clipper::ftype eps=ih.hkl_class().epsilon();
                if ( ih.hkl_class().centric() ) intensity_range_centric.include( (obs(hkldata[ih])/eps )/norm.f(ih) );
                else intensity_range_acentric.include( (obs(hkldata[ih])/eps ) /norm.f(ih) );
            }
        }
        
        int ncentric = 0;
        _intensity_ord_c.init( intensity_range_centric );
        _intensity_ord_a.init( intensity_range_acentric );
        for ( HRI ih = hkldata.first(); !ih.last(); ih.next() ) {
            if ( !hkldata[ih].missing() ) {
                clipper::ftype eps=ih.hkl_class().epsilon();
                if ( ih.hkl_class().centric() ) {
                    _intensity_ord_c.accumulate( ( obs(hkldata[ih])/eps )/norm.f(ih) );
                    ncentric++;
                }
                else _intensity_ord_a.accumulate( ( obs(hkldata[ih])/eps )/norm.f(ih) );
            }
        }
        _intensity_ord_c.prep_ordinal();
        _intensity_ord_a.prep_ordinal();
        
        clipper::ftype x = 0.08;
        clipper::ftype deltax = 0.12;
        for (int i=0;i!=4; ++i) {
            if ( (_acentricideal[3*i+2] - _intensity_ord_a.ordinal(x))/_acentricideal[3*i+2] > 0.4 ) _ntw ++;
            x += deltax;
        }
        
        return (_ntw > 2);
        
    }
    
    template<class T, template<class> class D> R_test::R_test(const clipper::HKL_data< D<T> >& hkldata, const clipper::Isymop& twin, clipper::Range<clipper::ftype> reso)
    {
        clipper::Range<clipper::ftype> range(hkldata.invresolsq_range() );
        this->operator()(hkldata, twin,range);
    }
    
    /*! calculate R-test.
     \param reso Resolution limit for calculation (default 3.0A)
     \return R average */
    template<class T, template<class> class D> clipper::ftype R_test::operator() (const clipper::HKL_data< D<T> >& isig, const clipper::Isymop& twin, clipper::Range<clipper::ftype> reso ) {
        typedef clipper::HKL_data_base::HKL_reference_index HRI;
        init(isig,false);
        _reso=reso;
        if (_reso.max() < 0 ) _reso = clipper::Range<clipper::ftype>(isig.hkl_info().invresolsq_range() );
        _twinop=twin;
        
        clipper::Spacegroup spgr = isig.hkl_info().spacegroup();
        
        double HT(0.0), NT(0.0);
        double scalefac = 12.0;
        _r.resize(this->size(),0.0);
        std::vector<double> nt(this->size(),0.0);
        clipper::Vec3<int> jhkl, jhkl2;
        for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
            if ( !isig[ih].missing() && !ih.hkl_class().centric() ) {
                if (reso.contains(ih.invresolsq() ) ) {
                    clipper::HKL hkl = ih.hkl();
                    jhkl[0] = hkl.h();
                    jhkl[1] = hkl.k();
                    jhkl[2] = hkl.l();
                    clipper::HKL twin;
                    jhkl2 = _twinop*jhkl;
                    twin.h() = jhkl2[0]/scalefac;
                    twin.k() = jhkl2[1]/scalefac;
                    twin.l() = jhkl2[2]/scalefac;
                    if (!isig[twin].missing()) {
                        double I1 = isig[ih].I()/ih.hkl_class().epsilon();
                        double I2 = isig[twin].I()/spgr.hkl_class(twin).epsilon();
                        clipper::ftype eps = (this->is_intensity() ) ? 1.0/ih.hkl_class().epsilon() : 1.0/std::sqrt(ih.hkl_class().epsilon());
                        clipper::ftype weight = 2.0/ih.hkl_class().epsilonc();
                        int bin = clipper::Util::bound( 0,clipper::Util::intf( clipper::ftype(this->size() ) * ResolStats_base::operator()(ih.invresolsq() ) ), this->size() );
                        double H = 0.0;
                        if ( I1 != 0.0 && I2 != 0.0) H = (I2-I1)/(I2+I1);
                        if (fabs(H) < 1){
                            HT += weight*fabs(H)*eps;
                            NT += weight;
                            _r[bin] += weight*fabs(H)*eps;
                            nt[bin] += weight;
                        }
                    }
                }
            }
        }
        for (int i = 0; i != this->size(); ++i ) _r[i] /= nt[i];
        _Rav = HT/NT;
        return _Rav;
    }
}

#endif

