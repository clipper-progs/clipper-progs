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
#include "intensity_scale.h"
#include "intensity_target.h"
#include "ctruncate_wilson.h"
#include "cpsf_utils.h"

namespace ctruncate {
	
	int cumulative_plot(clipper::HKL_data<clipper::data32::I_sigI>& isig, clipper::ResolutionFn& Sigma);
	
	int cumulative_plot(clipper::HKL_data<clipper::data32::I_sigI>& isig, clipper::HKL_data<clipper::data32::I_sigI>& Sigma);
			
	//--------------------------------------------------------------
	
	template<class D> class Completeness
	{
	public:
		enum TYPE { F, I };
		Completeness( int nbins=60 ) : 
		_nbins(nbins), 
		_compsig(nbins,0.0), _compsig1(nbins,0.0), _compsig2(nbins,0.0), _compsig3(nbins,0.0), _standard(nbins,0.0)
		{}
		void operator() (clipper::HKL_data<D>& fo, clipper::Resolution reso=clipper::Resolution() );
		clipper::ftype completeness(const clipper::ftype invresolsq) { return _compsig[int( double(_nbins) * invresolsq / _reso.invresolsq_limit() - 0.001)]; }
        clipper::ftype completeness2(const clipper::ftype invresolsq) { return _compsig2[int( double(_nbins) * invresolsq / _reso.invresolsq_limit() - 0.001)]; }
		clipper::ftype completeness3(const clipper::ftype invresolsq) { return _compsig3[int( double(_nbins) * invresolsq / _reso.invresolsq_limit() - 0.001)]; }
		clipper::ftype standard(const clipper::ftype invresolsq) { return _standard[int( double(_nbins) * invresolsq / _reso.invresolsq_limit() - 0.001)]; }
		clipper::ftype bin2invresolsq(const int i) { return _reso.invresolsq_limit()*(double(i)+0.5)/double(_nbins);; }
		void plot();
		
	private:
		template<class T> const T&    obs( const clipper::datatypes::F_sigF<T>& f ) { return f.f(); }
		template<class T> const T&    obs( const clipper::datatypes::I_sigI<T>& f ) { return f.I(); }
		template<class T> const T& sigobs( const clipper::datatypes::F_sigF<T>& f ) { return f.sigf(); }
		template<class T> const T& sigobs( const clipper::datatypes::I_sigI<T>& f ) { return f.sigI(); }
		template<class T> TYPE type( const clipper::datatypes::F_sigF<T>& f ) { return F; }
		template<class T> TYPE type( const clipper::datatypes::I_sigI<T>& f ) { return I; }
		
		TYPE _t;
		clipper::Resolution _reso;       //max resolution in calculation
		int _nbins;                    //number of bins
		
		std::vector<clipper::ftype> _compsig;
		std::vector<clipper::ftype> _compsig1;
		std::vector<clipper::ftype> _compsig2;
		std::vector<clipper::ftype> _compsig3;
		std::vector<clipper::ftype> _standard;
	};
	
		
	//--------------------------------------------------------------
	//! Check for tNCS using patterson map
	/*! Use patterson map to check for tNCA
	 \ingroup g_funcobj */
	
	class tNCS
	{
	public:
		//! operator
		template <class T, template<class> class D> const std::vector<clipper::Symop>& operator()(clipper::HKL_data<D<T> >& hkldata, clipper::Resolution reso=clipper::Resolution(4.0));
		//! has tNCS test
		bool hasNCS() const { return _peaks.size() != 0; }
		//! number of operators
		int numOps() const { return _peaks.size(); } 
		//! return peak height
		clipper::ftype height(int i) const { return _peak_height[i]; }
		//! return peak prob
		clipper::ftype prob(int i) const { return _peak_prob[i]; }
		//! print out text summary
		void output()const;
		//! output stringstream with xml
		std::stringstream& xml_output(std::stringstream&) const ;
        //! return tNCS operator as fractional coordinate
        const clipper::Symop& operator[](const int i) const { return _peaks[i]; }
		
	private:		
		template<class T> const T&    I( const clipper::datatypes::F_sigF<T>& f ) { return f.f()*f.f(); }
		template<class T> const T&    I( const clipper::datatypes::I_sigI<T>& f ) { return f.I(); }
		template<class T> const T&    I( const clipper::datatypes::F_sigF_ano<T>& f ) { return f.f()*f.f(); }
		template<class T> const T&    I( const clipper::datatypes::I_sigI_ano<T>& f ) { return f.I(); }
		
		clipper::HKL_data_base *_base;; //!< dataset
		clipper::Resolution _reso; //!< resolution limits
		std::vector<clipper::Symop> _peaks; //!< located peaks
		std::vector<clipper::ftype> _peak_height; //!< height of located peaks as fraction of origin
		std::vector<clipper::ftype> _peak_prob; //!< Zwartz estimation of probability
		
	};
	
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
	
	//ResolStats_base------------------------------------------------
	//! Base class for anomalous diffraction statistics
    /*! Gives resolution bins
     */
	class ResolStats_base {
	public:
		ResolStats_base(const clipper::HKL_data_base& hkldata, const bool missing = true, const int nrefln = 500) {
            clipper::Range<clipper::ftype> range;
            if (missing) {
                //_s_ord.init( hkldata, 1.0 );
                for (clipper::HKL_data_base::HKL_reference_index  ih = hkldata.first_data(); !ih.last(); hkldata.next_data(ih) )
                    range.include( ih.invresolsq() );
            } else {
                //_s_ord.init( hkldata.hkl_info(), 1.0);
                for (clipper::HKL_data_base::HKL_reference_index  ih = hkldata.hkl_info().first(); !ih.last(); ih.next() )
                    range.include( ih.invresolsq() );
            }
			init(hkldata,range,missing,nrefln);
		}
        template <class T> ResolStats_base(const clipper::HKL_data_base& hkldata, const clipper::Range<T>& range, const bool missing = true, const int nrefln = 500) {
            clipper::Range<clipper::ftype> tmp = range;
            tmp.include(range.max() ); tmp.include(range.min() );
            init(hkldata,tmp,missing,nrefln);
        }
		ResolStats_base(const ResolStats_base& base) {
			init(base);
		}
		ResolStats_base() {}
		~ResolStats_base() {}
		int operator()(const clipper::ftype) const;
		virtual clipper::ftype operator[](int index) const;
		ResolStats_base& operator=(const ResolStats_base&); 
		int size() const { return _b_reso.size(); }
		const clipper::HKL_data_base* parent() { return _base; }
        clipper::Range<clipper::ftype> binRange(const int i) {
            return _b_range[i];
        }
        clipper::ftype binContains(const int i) {
            return _b_contains[i];
        }
	protected:
		void init(const clipper::HKL_data_base& hkldata, const clipper::Range<clipper::ftype>&, const bool, const int);
		void init(const ResolStats_base&);
		bool is_intensity() { return _base->type() == "I_sigI_ano" || _base->type() == "J_sigJ_ano" || _base->type() == "I_sigI"; }
		bool is_anomalous() { return _base->type() == "I_sigI_ano" || _base->type() == "J_sigJ_ano" || _base-> type() == "F_sigF_ano"; }
	private:
		clipper::Generic_ordinal _s_ord;    //<! ordinal
        clipper::Range<clipper::ftype> _range; //<! range of analysis
		std::vector<clipper::ftype> _b_reso;   //<! mean resolution in bin
        std::vector<clipper::ftype> _b_contains; //<! total number in bin
        std::vector<clipper::Range<clipper::ftype> > _b_range; //< ! reso range of bin
		clipper::HKL_data_base *_base;
	};
	
	class AnomStats_measurability : public ResolStats_base {
	public:
		AnomStats_measurability(const clipper::HKL_data_base& hkldata);
		AnomStats_measurability(const ResolStats_base&);
		AnomStats_measurability() {}
		~AnomStats_measurability() {}
		clipper::ftype operator[](int index) const;
		AnomStats_measurability& operator=(const AnomStats_measurability&); 
	private:
		std::vector<clipper::ftype> _meas;
		
		void calc(const clipper::HKL_data_base& hkldata);
	};
	
	class AnomStats_bijveot : public ResolStats_base {
	public:
		AnomStats_bijveot(const clipper::HKL_data_base& hkldata);
		AnomStats_bijveot(const ResolStats_base&);
		AnomStats_bijveot() {}
		~AnomStats_bijveot() {}
		clipper::ftype operator[](int index) const;
		AnomStats_bijveot& operator=(const AnomStats_bijveot&); 
	private:
		std::vector<clipper::ftype> _meandI;
		
		void calc(const clipper::HKL_data_base& hkldata);
	};
	
	class AnomStats_signoise : public ResolStats_base {
	public:
		AnomStats_signoise(const clipper::HKL_data_base&);
		AnomStats_signoise(const ResolStats_base&);
		AnomStats_signoise() {}
		~AnomStats_signoise() {}
		clipper::ftype operator[](int index) const;
		AnomStats_signoise& operator=(const AnomStats_signoise&); 
	private:
		std::vector<clipper::ftype> _meandI;
		
		void calc(const clipper::HKL_data_base& hkldata);
	};
	
//! Perform analysis of Anomalous data
	/*! Wrap anomalous output
	 */
class AnomStats
{
public:
	//construct with data
    AnomStats(const clipper::HKL_data_base& hkldata);
	//return bijveot range
	const clipper::Range<clipper::ftype>& bijveot_range() { return _bij_range; }
	//return anomalous signal to noise
	const clipper::Range<clipper::ftype>& anom_signal_range() { return _signoise_range; }
	//return measureability
	const clipper::Range<clipper::ftype>& measurability_range() { return _meas_range; }
	//return output
	void output();
	//return xml
	std::stringstream& xml_output(std::stringstream&);
	//return data type (intensity or amplitudes)
	bool is_intensity() { return _base->type() == "I_sigI_ano" || _base->type() == "J_sigJ_ano"; }
private:
	clipper::HKL_data_base *_base;  //!< parent data object
	clipper::Range<clipper::ftype> _bij_range, _signoise_range, _meas_range;  //!< computed ranges
	ResolStats_base _binner; //!< resolution bins
	AnomStats_measurability _meas;  //!< measureability object
	AnomStats_bijveot _bij;         //!< biojvet object
	AnomStats_signoise _signoise;   //!< signal to noise object
};

	//--------------------------------------------------------------
	
	template<class D1, class D2> class ResoCorrel
	{
	public:
		enum TYPE { F, I };
		ResoCorrel( int nbins=60 ) : 
		_nbins(nbins), _comp(nbins,0.0)
		{}
		void operator() (clipper::HKL_data<D1>& fo, clipper::HKL_data<D2>&fc, clipper::Range<clipper::ftype> reso=clipper::Range<clipper::ftype>() );
		clipper::ftype CC_r(const clipper::ftype invresolsq) { return _comp[int( clipper::ftype(_nbins) * (invresolsq - _reso.min() ) / _reso.range() - 0.001)]; }
		clipper::ftype bin2invresolsq(const int i) { return _reso.min()+_reso.range()*(double(i)+0.5)/double(_nbins); }
		clipper::ftype CC() { return _cc; }
		int nbins() { return _nbins; }
		void plot();
		
	private:
		template<class T> const T&    obs( const clipper::datatypes::F_sigF<T>& f ) { return f.f(); }
		template<class T> const T&    obs( const clipper::datatypes::I_sigI<T>& f ) { return f.I(); }
		template<class T> const T& sigobs( const clipper::datatypes::F_sigF<T>& f ) { return f.sigf(); }
		template<class T> const T& sigobs( const clipper::datatypes::I_sigI<T>& f ) { return f.sigI(); }
		template<class T> TYPE type( const clipper::datatypes::F_sigF<T>& f ) { return F; }
		template<class T> TYPE type( const clipper::datatypes::I_sigI<T>& f ) { return I; }
		
		TYPE _t; //!< type intensity or structure factor
		clipper::Range<clipper::ftype> _reso;       //!< resolution limits in calculation, inverse resolution squared
		int _nbins;                    //!<number of bins
		
		std::vector<clipper::ftype> _comp; //!< correlation coeff array
		clipper::ftype _cc; //!< average CC
	};
	
	//NonAnomAnalysis-------------------------------------------------
	//! analyse data for completeness
	class HKLStats_completeness : public ResolStats_base {
	public: 
		HKLStats_completeness(const clipper::HKL_data_base& hkldata, clipper::ftype val=-99.0);
		HKLStats_completeness(const ResolStats_base&, clipper::ftype val=-99.0);
		HKLStats_completeness(clipper::ftype val=-99.0) { _val=val; }
		~HKLStats_completeness() {}
		clipper::ftype operator[](int index) const;
		clipper::ftype IoversigI() { return _val; }
		HKLStats_completeness& operator=(const HKLStats_completeness&); 
	private:
		std::vector<clipper::ftype> _completeness;  //!< completeness stats
		clipper::ftype _val;                //!< reference value
		
		void calc(const clipper::HKL_data_base& hkldata);
	};

	//NonAnomAnalysis-------------------------------------------------
	//! analyse data for completeness
	class HKLStats_Rstandard : public ResolStats_base {
	public: 
		HKLStats_Rstandard(const clipper::HKL_data_base& hkldata, clipper::ftype val=-99.0);
		HKLStats_Rstandard(const ResolStats_base&, clipper::ftype val=-99.0);
		HKLStats_Rstandard() {  }
		~HKLStats_Rstandard() {}
		clipper::ftype operator[](int index) const;
		HKLStats_Rstandard& operator=(const HKLStats_Rstandard&); 
	private:
		std::vector<clipper::ftype> _Rstandard;  //!< completeness stats
		
		void calc(const clipper::HKL_data_base& hkldata);
	};
	
    //----bins analysis----------------------------------------------
	//! Analyse for bins
	class Rings_analyse
	{
	public:
		//! contructor
		Rings_analyse() { }
        ~Rings_analyse() { }
		//! check for presence of  rings
		template <class T, template <class> class D> bool operator()(const clipper::HKL_data<D<T> >&, Rings&);
        //! check for presence of  rings vs WilsonB
        template <class T, template <class> class D> bool operator()(const clipper::HKL_data<D<T> >&, Rings&, ctruncate::WilsonB&);
		
	private:
		template<class T> const T&    obs( const clipper::datatypes::F_sigF<T>& f ) { return f.f(); }
		template<class T> const T&    obs( const clipper::datatypes::I_sigI<T>& f ) { return f.I(); }
		template<class T> const T&    obs( const clipper::datatypes::F_sigF_ano<T>& f ) { return f.f(); }
		template<class T> const T&    obs( const clipper::datatypes::I_sigI_ano<T>& f ) { return f.I(); }
		template<class T> const T& sigobs( const clipper::datatypes::F_sigF<T>& f ) { return f.sigf(); }
		template<class T> const T& sigobs( const clipper::datatypes::I_sigI<T>& f ) { return f.sigI(); }
		template<class T> const T& sigobs( const clipper::datatypes::F_sigF_ano<T>& f ) { return f.sigf(); }
		template<class T> const T& sigobs( const clipper::datatypes::I_sigI_ano<T>& f ) { return f.sigI(); }
		
    protected:
		ctruncate::Rings _ideal_rings;              //!< expected values in rings
		std::vector<clipper::ftype > _comp;        //!< completeness
		
		const clipper::HKL_data_base *_data;              //!< pointer to data
		ctruncate::Rings* _rings;                  //!< pointer to rings data
	};

	//----Rings analysis----------------------------------------------
	//! Analyse for ice rings
	class IceRings_analyse : public Rings_analyse
	{
	public:
		//! contruc
        IceRings_analyse(clipper::ftype tol=4.0, clipper::ftype ratioI=2.0, clipper::ftype ratioC=1.25 ) : _zTolerance(tol), _ratioI(ratioI), _ratioC(ratioC), _wB(NULL)  {
            _ice.DefaultIceRings();
            _ice.ClearSums();
        }
        ~IceRings_analyse() { }
		//! check for presence of  rings
		template <class T, template <class> class D> bool operator()(const clipper::HKL_data<D<T> >& data);
        //! check for presence of  rings against wilson B
        template <class T, template <class> class D> bool operator()(const clipper::HKL_data<D<T> >& data, ctruncate::WilsonB&);
		//!output summary
		std::string output();
		//!output xml
		std::stringstream& xml_output(std::stringstream&);
		//!ice rings?
		bool present();
        //!return rings data
        ctruncate::Rings& rings() { return _ice; }
        //!return or set tolerance for z score
        clipper::ftype tolerance() const { return _zTolerance; }
        void tolerance(clipper::ftype tol) { _zTolerance=tol; }
        //!return or set test ratio of intensities
        clipper::ftype ratioIntensity() const { return _ratioI; }
        void ratioIntensity(clipper::ftype ratio) { _ratioI=ratio; }
        //!return or set test ratio of completenesses
        clipper::ftype ratioComp() const { return _ratioC; }
        void ratioComp(clipper::ftype ratio) { _ratioC=ratio; }
		
	private:
		ctruncate::Rings _ice;                  //!< pointer to rings data
        ctruncate::WilsonB* _wB;                 //!< pointer to optional WilsonB
        clipper::ftype _zTolerance;              //!< tolerance for Z-score
        clipper::ftype _ratioI;                  //!< ratio of intensities
        clipper::ftype _ratioC;                  //!< ratio of completeness
	};
    
    //----Rings analysis----------------------------------------------
	//! Analyse for outlier rings
	class OutlierRings_analyse : public Rings_analyse
	{
	public:
		//! contructor
		OutlierRings_analyse(clipper::ftype tol=6.0, clipper::ftype ratioI=2.0, clipper::ftype ratioC=99.0 ) : _zTolerance(tol), _ratioI(ratioI), _ratioC(ratioC), _wB(NULL) {}
        ~OutlierRings_analyse() { }
		//! check for presence of  rings
        template <class T, template <class> class D> bool operator()(const clipper::HKL_data<D<T> >& data);
        //! check for presence of rings using wilson data
        template <class T, template <class> class D> bool operator()(const clipper::HKL_data<D<T> >& data, ctruncate::WilsonB&);
		//!output summary
		std::string output();
		//!output xml
		std::stringstream& xml_output(std::stringstream&);
		//!ice rings?
		bool present();
        //!return rings data
        ctruncate::Rings& rings() { return _outliers; }
        //!return or set tolerance for z score
        clipper::ftype tolerance() const { return _zTolerance; }
        void tolerance(clipper::ftype tol) { _zTolerance=tol; }
        //!return or set test ratio of intensities
        clipper::ftype ratioIntensity() const { return _ratioI; }
        void ratioIntensity(clipper::ftype ratio) { _ratioI=ratio; }
        //!return or set test ratio of completenesses
        clipper::ftype ratioComp() const { return _ratioC; }
        void ratioComp(clipper::ftype ratio) { _ratioC=ratio; }

		
	private:
		ctruncate::Rings _outliers;                  //!< pointer to rings data
        ctruncate::WilsonB* _wB;                 //!< pointer to optional WilsonB
        clipper::ftype _zTolerance;              //!< tolerance for Z-score
        clipper::ftype _ratioI;                  //!< ratio of intensities
        clipper::ftype _ratioC;                  //!< ratio of completeness
	};

	//NonAnomAnalysis-------------------------------------------------
	//! analyse data for completeness, wilson temp and ice rings
	class HKLAnalysis {
	public:
		//constructor
        HKLAnalysis() {}
		template<class T, template<class> class D> HKLAnalysis(const clipper::HKL_data< D<T> >&);
        ~HKLAnalysis() {}
		//return an estimate of the active range
		const clipper::Range<clipper::ftype>& active_range();
		//output to std::out
		void output();
		//return xml
		std::stringstream& xml_output(std::stringstream&);
		//return data type (intensity or amplitudes)
		bool is_intensity() { return _data->type() == "I_sigI_ano" || _data->type() == "J_sigJ_ano" || _data->type() == "I_sigI"; }
		bool is_anomalous() { return _data->type() == "I_sigI_ano" || _data->type() == "J_sigJ_ano" || _data-> type() == "F_sigF_ano"; }
		//return has ice rings
		bool is_iced() { return _ira.present();}
		//return ice rings
		Rings& ice_rings() { return _ira.rings(); }
        //return wilson B factor
        clipper::ftype wilson_Bfactor() { return _wilsonB.B(); }
        //return wilson intercept
        clipper::ftype wilson_intercept() { return _wilsonB.intercept(); }
	private:
		const clipper::HKL_data_base *_data;              //!< pointer to data
		
		ResolStats_base _binner; //!< resolution bins
		std::vector<HKLStats_completeness> _completeness;  //!< completeness object
		std::vector<clipper::Range<clipper::ftype> > _activerange; //!< resolution range over which we have 85% of data
		HKLStats_Rstandard _Rstandard; //!< rstandard object
        ResolStats_base _lb; //!< low resolution bins
        HKLStats_completeness _lcompleteness; //!< low resolution completeness object
		
		IceRings_analyse _ira; //!< analysis of ice rings
        
        OutlierRings_analyse _ora; //!< analysis for outliers
		
		static clipper::ftype ACCEPTABLE;       //!< acceptable stat for active range
		clipper::Range<clipper::ftype> _active; //!< store active range
		
		WilsonB _wilsonB; //!<< wilson B plot
	};
	
	
	   
	//template members----------------------------------------------
	/*! constructor for the HKL analysis.
     \param hkldata Experimental intensities
     \return type HKLAnalysis */
    template<class T, template<class> class D> HKLAnalysis::HKLAnalysis( const clipper::HKL_data< D<T> >& hkldata ) {
		_data = &hkldata;
		_completeness.resize(7);
		_activerange.resize(7);
		
		bool missing=false;
		_binner  = ResolStats_base(hkldata,missing);
		_completeness[0]    = HKLStats_completeness(_binner);
		_completeness[1]    = HKLStats_completeness(_binner,1.0);
		_completeness[2]    = HKLStats_completeness(_binner,2.0);
		_completeness[3]    = HKLStats_completeness(_binner,3.0);
		_completeness[4]    = HKLStats_completeness(_binner,5.0);
		_completeness[5]    = HKLStats_completeness(_binner,10.0);
		_completeness[6]    = HKLStats_completeness(_binner,15.0);
		_Rstandard = HKLStats_Rstandard(_binner);
				
        int nbins = _binner.size();
		
		int NBINS = _binner.size();
		for (int ii=0; ii != _completeness.size() ; ++ii) {
			if ((_completeness[ii])[0] >= ACCEPTABLE) _activerange[ii].include(hkldata.invresolsq_range().min());
			if ((_completeness[ii])[NBINS-1] >= ACCEPTABLE) _activerange[ii].include(hkldata.invresolsq_range().max());
			for ( int i=1 ; i != NBINS ; ++i) {
				int i1=i-1;
				if ( (_completeness[ii])[i] >= ACCEPTABLE && (_completeness[ii])[i1] < ACCEPTABLE ) {
                    //_activerange[ii].include(0.5*(_binner[i]+_binner[i1]) );
					_activerange[ii].include(
                                             _binner[i1]+(ACCEPTABLE-(_completeness[ii])[i1])/((_completeness[ii])[i]-(_completeness[ii])[i1])*(_binner[i]-_binner[i1])
                    );
				} else if ( (_completeness[ii])[i] < ACCEPTABLE && (_completeness[ii])[i1] >= ACCEPTABLE ) {
					//_activerange[ii].include(0.5*(_binner[i]+_binner[i1]) );
                    _activerange[ii].include(
                                             _binner[i1]+(ACCEPTABLE-(_completeness[ii])[i1])/((_completeness[ii])[i]-(_completeness[ii])[i1])*(_binner[i]-_binner[i1]) );
				}
			}
		}
		
        {
            clipper::Range<clipper::ftype> lrange;
            lrange.include((_binner.binRange(0) ).min() );
            if (_binner.size() > 2 ) lrange.include( (_binner.binRange(1) ).max() );
            else lrange.include( (_binner.binRange(0) ).max() );
            //low resolution completeness
            _lb = ResolStats_base(hkldata,lrange,missing,100);
            _lcompleteness = HKLStats_completeness(_lb);
        }
        
		//ice rings
		_ira(hkldata);
		
        Rings ice=_ira.rings();
		//wilson plot (how to handle contents)
		_wilsonB(hkldata,&(this->active_range()),&ice );
        
        //IceRings_analyse test;
        //test(hkldata,_wilsonB);
        
        _ora(hkldata,_wilsonB);
		
		return;
	}
	
    //----Ice Rings analysis----------------------------------------------
	/*! operator to do ice rings analysis
	 \param data Reflection data
	 \param rings rings to be analysed
	 \return rings rejected?
	 */
	template<class T, template<class> class D> bool IceRings_analyse::operator()(const clipper::HKL_data< D<T> >& data)
	{
        for (int i = 0; i != _ice.Nrings(); ++i) _ice.SetReject(i, true);
        
        this->Rings_analyse::operator()(data,_ice);
        
        if (_zTolerance < 1.0) _zTolerance = 1.0/_zTolerance;
        if (_ratioI < 1.0) _ratioI = 1.0/_ratioI;
        if (_ratioC < 1.0) _ratioI = 1.0/_ratioC;
        
        for (int i = 0; i != _ice.Nrings(); ++i) {
            bool reject = false;
            float reso = _ice.MeanSSqr(i);
            if ( reso > 0.0 ) {
                clipper::ftype r1 = _ice.MeanI(i)/_ideal_rings.MeanI(i);
                clipper::ftype r2 = _ice.Comp(i)/_comp[i];
                if ((std::abs(_ice.MeanI(i)-_ideal_rings.MeanI(i))/_ice.MeanSigI(i) > _zTolerance) ||
                    ( r1 >= _ratioI || 1.0/r1 > _ratioI) ||
                    ( r2 >= _ratioC || 1.0/r2 > _ratioC) ) reject = true;
            }
            _ice.SetReject(i, reject);
        }
        
        for ( int i = 0; i != _ice.Nrings(); ++i)
            if (_rings->Reject(i) ) return true;
        
        return false;
    }
    
    //----Ice Rings analysis----------------------------------------------
    /*! operator to do ice rings analysis
     \param data Reflection data
     \param rings rings to be analysed
     \param wilsonB object
     \return rings rejected?
     */
    template<class T, template<class> class D> bool IceRings_analyse::operator()(const clipper::HKL_data< D<T> >& data, ctruncate::WilsonB& wilson)
    {
        _wB = &wilson;
        
        for (int i = 0; i != _ice.Nrings(); ++i) _ice.SetReject(i, true);
        
        this->Rings_analyse::operator()(data,_ice,wilson);
        
        if (_zTolerance < 1.0) _zTolerance = 1.0/_zTolerance;
        if (_ratioI < 1.0) _ratioI = 1.0/_ratioI;
        if (_ratioC < 1.0) _ratioI = 1.0/_ratioC;
        
        for (int i = 0; i != _ice.Nrings(); ++i) {
            bool reject = false;
            float reso = _ice.MeanSSqr(i);
            if ( reso > 0.0 ) {
                clipper::ftype r1 = _ice.MeanI(i)/_ideal_rings.MeanI(i);
                clipper::ftype r2 = _ice.Comp(i)/_comp[i];
                if ((std::abs(_ice.MeanI(i)-_ideal_rings.MeanI(i))/_ice.MeanSigI(i) > _zTolerance) ||
                    ( r1 >= _ratioI || 1.0/r1 > _ratioI) ||
                    ( r2 >= _ratioC || 1.0/r2 > _ratioC) ) reject = true;
            }
            _ice.SetReject(i, reject);
        }
        
        for ( int i = 0; i != _ice.Nrings(); ++i)
            if (_rings->Reject(i) ) return true;
        
        return false;
    }

    
	//----Rings analysis----------------------------------------------
	/*! operator to do rings analysis
	 \param data Reflection data
	 \param rings rings to be analysed
	 \return rings rejected?
	 */
	template<class T, template<class> class D> bool Rings_analyse::operator()(const clipper::HKL_data< D<T> >& data, Rings& rings)
	{
        _rings = &rings;
        
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		
		int nr = _rings->Nrings();
		_ideal_rings.Copy(*_rings);
		_data = &data;
		_comp.resize(_rings->Nrings());
		
		//_rings->ClearSums();
		_ideal_rings.ClearSums();
		clipper::HKL_data<D<T> > xsig(data.hkl_info() );
		// dataset with all ice rings removed
		for ( HRI ih = data.first(); !ih.last(); ih.next() ) {
			double reso = ih.invresolsq();
			xsig[ih] = D<T>( obs(data[ih]), sigobs(data[ih]) );
            int ring = rings.InRing(reso);
            if ( ring != -1 ) {
				if ( rings.Reject( ring ) ) xsig[ih].I() = xsig[ih].sigI() = clipper::Util::nan();
                clipper::ftype mult=data.hkl_info().spacegroup().num_symops()/ih.hkl_class().epsilon();
                if (ih.hkl_class().centric() ) {
                   _rings->AddObs(ring,D<T>(0.5*data[ih].I(),0.5*data[ih].sigI() ),reso,mult);
                } else {
                    _rings->AddObs(ring,data[ih],reso,mult);
                }
            }
		}
        
        clipper::Range<clipper::ftype> range=data.invresolsq_range();
        clipper::Generic_ordinal s_ord;
        s_ord.init( range, 1000 );
        for (clipper::HKL_data_base::HKL_reference_index ih = xsig.first_data(); !ih.last(); xsig.next_data(ih) ) {
            s_ord.accumulate( ih.invresolsq() );
        }
        s_ord.prep_ordinal();
        
        int nbins(0);
        int Nreflections(0);
        int nreflns(500);
        
        for (clipper::HKL_data_base::HKL_reference_index ih = xsig.first_data(); !ih.last(); xsig.next_data(ih) )
            if (!xsig[ih].missing()  ) ++Nreflections;
        
        nbins = std::max( (Nreflections/nreflns) , 1);
        
		std::vector<float> summeas(nbins,0.0), sumov(nbins,0.0), mean1(nbins,0.0), means(nbins,0.0), var(nbins,0.0);
        // completeness
        {
            for ( HRI ih = data.hkl_info().first(); !ih.last(); ih.next() ) {
                clipper::ftype reso = ih.invresolsq();
                if ( _rings->InRing(reso) == -1 ) {
                    clipper::ftype mult=data.hkl_info().spacegroup().num_symops()/ih.hkl_class().epsilon();
                    int bin = clipper::Util::bound( 0,clipper::Util::intf( clipper::ftype(nbins) * s_ord.ordinal( reso ) ), nbins-1 );
                    clipper::ftype s = ih.invresolsq();
                    sumov[bin] += mult;
                    if ( !data[ih].missing() ) {
                        if (ih.hkl_class().centric() ) {
                            var[bin] += 0.5*mult*data[ih].I();
                            means[bin] += 0.5*mult*data[ih].sigI()*data[ih].sigI();
                        } else {
                            var[bin] += mult*data[ih].I();
                            means[bin] += mult*data[ih].sigI()*data[ih].sigI();
                        }
                        summeas[bin] += mult;
                    }
                }
            }
            for (int i=0 ; i != _rings->Nrings() ; ++ i) {
                float reso = _rings->MeanSSqr(i);
                int bin = clipper::Util::bound( 0,clipper::Util::intf( clipper::ftype(nbins) * s_ord.ordinal( reso ) ), nbins-1 );
                _comp[i] = summeas[bin]/sumov[bin];
            }
        }
        
		std::vector<double> params_ice( nbins, 1.0 );
		clipper::BasisFn_spline basis_ice( xsig, nbins, 2.0 );
		TargetFn_meanInth<D<T> > target_ice( xsig, 1 );
		clipper::ResolutionFn mean( xsig.hkl_info(), basis_ice, target_ice, params_ice );
		
		
		// repeat with ibest
		for ( HRI ih = xsig.hkl_info().first(); !ih.last(); ih.next() ) {
			clipper::ftype reso = ih.invresolsq();
            clipper::ftype eps = ih.hkl_class().epsilon();
			clipper::ftype mult=xsig.hkl_info().spacegroup().num_symops()/eps;
			int ring=_ideal_rings.InRing(reso);
			if ( ring != -1 ) {
                if (ih.hkl_class().centric() ) {
                    _ideal_rings.AddObs(ring,D<T>(0.5*eps*mean.f(ih),0.0 ),reso,mult);
                } else {
                    _ideal_rings.AddObs(ring,D<T>(eps*mean.f(ih),0.0 ),reso,mult);
                }
            }
		}
        
        for ( HRI ih = xsig.hkl_info().first(); !ih.last(); ih.next() ) {
            clipper::ftype reso = ih.invresolsq();
            if ( _rings->InRing(reso) == -1 ) {
                if ( !xsig[ih].missing() ) {
                    int bin = clipper::Util::bound( 0,clipper::Util::intf( clipper::ftype(nbins) * s_ord.ordinal( reso ) ), nbins-1 );
                    clipper::ftype mult=data.hkl_info().spacegroup().num_symops();
                    if (ih.hkl_class().centric() ) {
                        mean1[bin] += mult*0.5*mean.f(ih);
                    } else {
                        mean1[bin] += mult*mean.f(ih);
                    }
                }
            }
        }
        clipper::ftype v1(0.0),z1(0.0);
        for (int i=0 ; i != nbins; ++i ) {
            //std::cout << mean1[i]/summeas[i] << " - " << var[i]/summeas[i] << " - " << std::sqrt(means[i]/summeas[i]) << " -- " << ((mean1[i]-var[i])/summeas[i])/std::sqrt(means[i]/summeas[i])  << std::endl;
            clipper::ftype t =((mean1[i]-var[i])/summeas[i])/std::sqrt(means[i]/summeas[i]);
            z1 += t;
            v1 += t*t;
        }
        //std::cout << "global stats: " << z1/nbins << " --- " << std::sqrt(nbins*v1-z1*z1)/nbins << std::endl;
        return true;
	}
	
    //----Rings analysis----------------------------------------------
    /*! operator to do rings analysis
     \param data Reflection data
     \param rings rings to be analysed
     \param wilson WilsonB object (must be after calculation)
     \return rings rejected?
     */
    template<class T, template<class> class D> bool Rings_analyse::operator()(const clipper::HKL_data< D<T> >& data, Rings& rings, ctruncate::WilsonB& wilson)
    {
        _rings = &rings;
        
        typedef clipper::HKL_data_base::HKL_reference_index HRI;
        
        int nr = _rings->Nrings();
        _ideal_rings.Copy(*_rings);
        _data = &data;
        _comp.resize(_rings->Nrings());
        
        //_rings->ClearSums();
        _ideal_rings.ClearSums();
        
        int nbins(0);
        int Nreflections(0);
        int nreflns(500);
        clipper::Range<clipper::ftype> range=data.invresolsq_range();
        clipper::Generic_ordinal s_ord;
        s_ord.init( range, 1000 );
        for (clipper::HKL_data_base::HKL_reference_index ih = data.first_data(); !ih.last(); data.next_data(ih) ) {
            clipper::ftype reso = ih.invresolsq();
            int ring = _rings->InRing(reso);
            if ( ring == -1 ) {
                s_ord.accumulate( reso );
            } else {
                if ( !rings.Reject( ring ) ) s_ord.accumulate( reso );
            }
        }
        s_ord.prep_ordinal();
        
        for (clipper::HKL_data_base::HKL_reference_index ih = data.first_data(); !ih.last(); data.next_data(ih) ) {
            clipper::ftype reso = ih.invresolsq();
            int ring = _rings->InRing(reso);
            if ( ring == -1 ) {
                if (!data[ih].missing()  ) ++Nreflections;
            } else {
                if ( !rings.Reject( ring ) && !data[ih].missing()  ) ++Nreflections;
            }
        }
        nbins = std::max( (Nreflections/nreflns) , 1);
        
        std::vector<float> summeas(nbins,0.0), sumov(nbins,0.0), mean1(nbins,0.0), means(nbins,0.0), var(nbins,0.0);
        // repeat with ibest
        for ( HRI ih = data.hkl_info().first(); !ih.last(); ih.next() ) {
            clipper::ftype reso = ih.invresolsq();
            clipper::ftype eps = ih.hkl_class().epsilon();
            clipper::ftype mult=data.hkl_info().spacegroup().num_symops()/eps;
            int ring=_ideal_rings.InRing(reso);
            if ( ring != -1 ) {
                if (ih.hkl_class().centric() ) {
                    _ideal_rings.AddObs(ring,D<T>(0.5*eps*wilson.f(ih),0.0 ),reso,mult);
                    _rings->AddObs(ring,D<T>(0.5*data[ih].I(),0.5*data[ih].sigI() ),reso,mult);
                    
                } else {
                    _ideal_rings.AddObs(ring,D<T>(eps*wilson.f(ih),0.0 ),reso,mult);
                    _rings->AddObs(ring,data[ih],reso,mult);
                }
                if ( !rings.Reject( ring )  && !data[ih].missing() ) {
                    int bin = clipper::Util::bound( 0,clipper::Util::intf( clipper::ftype(nbins) * s_ord.ordinal( reso ) ), nbins-1 );
                    clipper::ftype mult=data.hkl_info().spacegroup().num_symops();
                    if (ih.hkl_class().centric() ) {
                        mean1[bin] += mult*0.5*wilson.f(ih);
                    } else {
                        mean1[bin] += mult*wilson.f(ih);
                    }
                }
            } else {
                if ( !data[ih].missing() ) {
                    int bin = clipper::Util::bound( 0,clipper::Util::intf( clipper::ftype(nbins) * s_ord.ordinal( reso ) ), nbins-1 );
                    clipper::ftype mult=data.hkl_info().spacegroup().num_symops();
                    if (ih.hkl_class().centric() ) {
                        mean1[bin] += mult*0.5*wilson.f(ih);
                    } else {
                        mean1[bin] += mult*wilson.f(ih);
                    }
                }
            }
        }

        // completeness
        {
            for ( HRI ih = data.hkl_info().first(); !ih.last(); ih.next() ) {
                clipper::ftype reso = ih.invresolsq();
                int ring = _rings->InRing(reso);
                if ( ring == -1 ) {
                    clipper::ftype mult=data.hkl_info().spacegroup().num_symops()/ih.hkl_class().epsilon();
                    int bin = clipper::Util::bound( 0,clipper::Util::intf( clipper::ftype(nbins) * s_ord.ordinal( reso ) ), nbins-1 );
                    clipper::ftype s = ih.invresolsq();
                    sumov[bin] += mult;
                    if ( !data[ih].missing() ) {
                        if (ih.hkl_class().centric() ) {
                            var[bin] += 0.5*mult*data[ih].I();
                            means[bin] += 0.5*mult*data[ih].sigI()*data[ih].sigI();
                        } else {
                            var[bin] += mult*data[ih].I();
                            means[bin] += mult*data[ih].sigI()*data[ih].sigI();
                        }
                        summeas[bin] += mult;
                    }
                } else {
                    if ( !rings.Reject( ring ) ) {
                        clipper::ftype mult=data.hkl_info().spacegroup().num_symops()/ih.hkl_class().epsilon();
                        int bin = clipper::Util::bound( 0,clipper::Util::intf( clipper::ftype(nbins) * s_ord.ordinal( reso ) ), nbins-1 );
                        clipper::ftype s = ih.invresolsq();
                        sumov[bin] += mult;
                        if ( !data[ih].missing() ) {
                            if (ih.hkl_class().centric() ) {
                                var[bin] += 0.5*mult*data[ih].I();
                                means[bin] += 0.5*mult*data[ih].sigI()*data[ih].sigI();
                            } else {
                                var[bin] += mult*data[ih].I();
                                means[bin] += mult*data[ih].sigI()*data[ih].sigI();
                            }
                            summeas[bin] += mult;
                        }
                    }
                }
            }
            for (int i=0 ; i != _rings->Nrings() ; ++ i) {
                float reso = _rings->MeanSSqr(i);
                int bin = clipper::Util::bound( 0,clipper::Util::intf( clipper::ftype(nbins) * s_ord.ordinal( reso ) ), nbins-1 );
                if (bin == 0 ) _comp[0] = 0.5*(summeas[1]/sumov[1]+summeas[2]/sumov[2]);
                else if (bin == (nbins-1)) _comp[nbins-1] = 0.5*(summeas[nbins-2]/sumov[nbins-2]+summeas[nbins-3]/sumov[nbins-3]);
                else _comp[bin] = 0.5*(summeas[bin-1]/sumov[bin-1]+summeas[bin+1]/sumov[bin+1]);
            }
        }
        clipper::ftype v1(0.0),z1(0.0);
        for (int i=0 ; i != nbins; ++i ) {
            //std::cout << mean1[i]/summeas[i] << " - " << var[i]/summeas[i] << " - " << std::sqrt(means[i]/summeas[i]) << " -- " << ((mean1[i]-var[i])/summeas[i])/std::sqrt(means[i]/summeas[i])  << std::endl;
            clipper::ftype t =((mean1[i]-var[i])/summeas[i])/std::sqrt(means[i]/summeas[i]);
            z1 += t;
            v1 += t*t;
        }
        //std::cout << "global stats: " << z1/nbins << " --wilson-- " << std::sqrt(nbins*v1-z1*z1)/nbins << std::endl;
        return true;
    }

	//-------tNCS peak search-------------------------------------------
	
	template <class T, template<class> class D> const std::vector<clipper::Symop>& tNCS::operator() (clipper::HKL_data<D<T> >& hkldata, clipper::Resolution r)
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
		hklp.init( pspgr, cell, _reso, true );
		
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
		clipper::Grid_sampling grid( pspgr, cell, _reso );
		
		// make xmap
		clipper::Xmap<float> patterson( pspgr, cell, grid );
		patterson.fft_from( fphi );
		
		//peak search in patterson
		PeakSearch pksch;                      // peak search object
		
		int npeak = 10;
		
		const std::vector<int>& ppks = pksch( patterson );
		
		if (ppks.size() == 0 || ppks.size() == 1 ) return _peaks;
		
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
	}
	
    //----Rings analysis----------------------------------------------
	/*! operator to do ice rings analysis
	 \param data Reflection data
	 \param rings rings to be analysed
	 \return rings rejected?
	 */
    template<class T, template<class> class D> bool OutlierRings_analyse::operator()(const clipper::HKL_data< D<T> >& data)
	{
        clipper::Range<clipper::ftype> range=data.invresolsq_range();
        clipper::Generic_ordinal s_ord;
        s_ord.init( range, 1000 );
        for (clipper::HKL_data_base::HKL_reference_index ih = data.hkl_info().first(); !ih.last(); ih.next() ) {
                s_ord.accumulate( ih.invresolsq() );
        }
        s_ord.prep_ordinal();
        
        int nbins(0);
        int Nreflections(0);
        int nreflns(200);
        
        for (clipper::HKL_data_base::HKL_reference_index ih = data.hkl_info().first(); !ih.last(); ih.next() )
            ++Nreflections;
        
        {
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
        
        std::vector<clipper::ftype> bmax(nbins,0.0), bmin(nbins,99999.0);
        for (clipper::HKL_data_base::HKL_reference_index ih = data.hkl_info().first(); !ih.last(); ih.next() ) {
            clipper::ftype s = ih.invresolsq();
            int bin = clipper::Util::bound( 0,clipper::Util::intf( clipper::ftype(nbins) * s_ord.ordinal( s ) ), nbins-1 );
            if ( s > bmax[bin] ) bmax[bin] = s;
            if ( s < bmin[bin] ) bmin[bin] = s;
        }
        
        for (int i=0 ; i != nbins; ++i ) _outliers.AddRing( (0.5*(1.0/std::sqrt(bmax[i])+1.0/std::sqrt(bmin[i]) ) ), std::fabs(0.5*(bmax[i]-bmin[i])) );
        
        _outliers.ClearSums();
        this->Rings_analyse::operator()(data,_outliers);
    
        if (_zTolerance < 1.0) _zTolerance = 1.0/_zTolerance;
        if (_ratioI < 1.0) _ratioI = 1.0/_ratioI;
        if (_ratioC < 1.0) _ratioI = 1.0/_ratioC;
        
        for (int i = 0; i != _outliers.Nrings(); ++i) {
            bool reject = false;
            float reso = _outliers.MeanSSqr(i);
            if ( reso > 0.0 ) {
                clipper::ftype r1 = _outliers.MeanI(i)/_ideal_rings.MeanI(i);
                clipper::ftype r2 = _outliers.Comp(i)/_comp[i];
                if ((std::abs(_outliers.MeanI(i)-_ideal_rings.MeanI(i))/_outliers.MeanSigI(i) > _zTolerance) &&
                    ( r1 >= _ratioI || 1.0/r1 > _ratioI) &&
                    ( r2 >= _ratioC || 1.0/r2 > _ratioC) ) reject = true;
            }
            _outliers.SetReject(i, reject);
        }

        
        for ( int i = 0; i != _outliers.Nrings(); ++i)
            if (_rings->Reject(i) ) return true;
        
        return false;

    }

    //----Rings analysis----------------------------------------------
    /*! operator to do ice rings analysis
     \param data Reflection data
     \param rings rings to be analysed
     \return rings rejected?
     */
    template<class T, template<class> class D> bool OutlierRings_analyse::operator()(const clipper::HKL_data< D<T> >& data, ctruncate::WilsonB& wilson)
    {
        _wB = &wilson;
        
        clipper::Range<clipper::ftype> range=data.invresolsq_range();
        clipper::Generic_ordinal s_ord;
        s_ord.init( range, 1000 );
        for (clipper::HKL_data_base::HKL_reference_index ih = data.hkl_info().first(); !ih.last(); ih.next() ) {
            s_ord.accumulate( ih.invresolsq() );
        }
        s_ord.prep_ordinal();
        
        int nbins(0);
        int Nreflections(0);
        int nreflns(200);
        
        for (clipper::HKL_data_base::HKL_reference_index ih = data.hkl_info().first(); !ih.last(); ih.next() )
            ++Nreflections;
        
        {
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
        
        std::vector<clipper::ftype> bmax(nbins,0.0), bmin(nbins,99999.0);
        for (clipper::HKL_data_base::HKL_reference_index ih = data.hkl_info().first(); !ih.last(); ih.next() ) {
            clipper::ftype s = ih.invresolsq();
            int bin = clipper::Util::bound( 0,clipper::Util::intf( clipper::ftype(nbins) * s_ord.ordinal( s ) ), nbins-1 );
            if ( s > bmax[bin] ) bmax[bin] = s;
            if ( s < bmin[bin] ) bmin[bin] = s;
        }
        
        for (int i=0 ; i != nbins; ++i ) _outliers.AddRing( (0.5*(1.0/std::sqrt(bmax[i])+1.0/std::sqrt(bmin[i]) ) ), std::fabs(0.5*(bmax[i]-bmin[i])) );
        
        _outliers.ClearSums();
        this->Rings_analyse::operator()(data,_outliers,wilson);
        
        if (_zTolerance < 1.0) _zTolerance = 1.0/_zTolerance;
        if (_ratioI < 1.0) _ratioI = 1.0/_ratioI;
        if (_ratioC < 1.0) _ratioI = 1.0/_ratioC;
        
        for (int i = 0; i != _outliers.Nrings(); ++i) {
            bool reject = false;
            float reso = _outliers.MeanSSqr(i);
            if ( reso > 0.0 ) {
                clipper::ftype r1 = _outliers.MeanI(i)/_ideal_rings.MeanI(i);
                clipper::ftype r2 = _outliers.Comp(i)/_comp[i];
                if ((std::abs(_outliers.MeanI(i)-_ideal_rings.MeanI(i))/_outliers.MeanSigI(i) > _zTolerance) &&
                    ( r1 >= _ratioI || 1.0/r1 > _ratioI) &&
                    ( r2 >= _ratioC || 1.0/r2 > _ratioC) ) reject = true;
            }
            _outliers.SetReject(i, reject);
        }
        
        for ( int i = 0; i != _outliers.Nrings(); ++i)
            if (_rings->Reject(i) ) return true;
        
        return false;
        
    }
    
}

#endif
