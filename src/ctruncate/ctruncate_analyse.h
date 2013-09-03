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

namespace ctruncate {
	
	int cumulative_plot(clipper::HKL_data<clipper::data32::I_sigI>& isig, clipper::ResolutionFn& Sigma);
	
	int cumulative_plot(clipper::HKL_data<clipper::data32::I_sigI>& isig, clipper::HKL_data<clipper::data32::I_sigI>& Sigma);
	
	void yorgo_modis_plot(clipper::HKL_data<clipper::data32::F_sigF>& fsig, float maxres, int nbins, CCP4Program& prog, clipper::U_aniso_orth uao = clipper::U_aniso_orth(1.0));
	
	//void yorgo_modis_plot(clipper::HKL_data<clipper::data32::I_sigI>& isig, float maxres, int nbins, CCP4Program& prog);
	
	void yorgo_modis_plot(clipper::HKL_data<clipper::data32::I_sigI>& isig, float maxres, int nbins, CCP4Program& prog, clipper::U_aniso_orth uao = clipper::U_aniso_orth(1.0) );
	
	//--------------------------------------------------------------
	
	template<class D> class YorgoModis
	{
	public:
		enum TYPE { F, I };
		YorgoModis( clipper::ftype maxres=3.0, int nbins=60, clipper::U_aniso_orth uao = clipper::U_aniso_orth(1.0) ) : _nbins(nbins), _uao(uao),
		_somov(nbins,0.0), _somsdov(nbins,0.0), _numov(nbins,0), _enumov(nbins,0.0),
		_somdir(3*nbins,0.0), _somsddir(3*nbins,0.0), _numdir(3*nbins,0.0), _enumdir(3*nbins,0.0)
		{}
		void operator() (clipper::HKL_data<D>& fo, clipper::Resolution reso=clipper::Resolution());
		void plot();
		
	private:
		template<class T> const T&    obs( const clipper::datatypes::F_sigF<T>& f ) { return f.f(); }
		template<class T> const T&    obs( const clipper::datatypes::I_sigI<T>& f ) { return f.I(); }
		template<class T> const T& sigobs( const clipper::datatypes::F_sigF<T>& f ) { return f.sigf(); }
		template<class T> const T& sigobs( const clipper::datatypes::I_sigI<T>& f ) { return f.sigI(); }
		template<class T> TYPE type( const clipper::datatypes::F_sigF<T>& f ) { return F; }
		template<class T> TYPE type( const clipper::datatypes::I_sigI<T>& f ) { return I; }
		
		TYPE _t;
		clipper::Resolution _reso;       //resolution in calculation
		int _nbins;                    //number of bins
		clipper::U_aniso_orth _uao;   //anisotropy matrix for directions
		
		clipper::Mat33<clipper::ftype> _e123; //directions
		
		std::vector<clipper::ftype> _somov;
		std::vector<clipper::ftype> _somsdov;
		std::vector<int> _numov;
		std::vector<clipper::ftype> _enumov;
		std::vector<clipper::ftype> _somdir;    
		std::vector<clipper::ftype> _somsddir;
		std::vector<int>  _numdir;
		std::vector<clipper::ftype>  _enumdir;
			
		int _nzerosigma;
	};
		
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
	
	template<class T> class AnisoPlot
	{
	public:
		AnisoPlot( clipper::U_aniso_orth& uao);
		AnisoPlot( clipper::ftype scale, clipper::U_aniso_orth& uao);
		void loggraph();
		
	private:
		std::vector<T> _eigen; //eigenvalues
		clipper::Mat33<T> _e123; //directions
		
		std::vector<std::vector<clipper::Coord_orth> > _isoline1;
		std::vector<std::vector<clipper::Coord_orth> > _isoline2;
		std::vector<std::vector<clipper::Coord_orth> > _isoline3;
		
		std::vector<clipper::Coord_orth> isoline(clipper::ftype offsetx, clipper::ftype offsety, clipper::ftype sigmau, 
					 clipper::ftype sigmav, clipper::ftype angleuv, clipper::ftype frac=1.0/1.17741, int steps=60);
		
		clipper::Coord_orth point(clipper::ftype offsetx, clipper::ftype offsety, clipper::ftype sigmau, 
								  clipper::ftype sigmav, clipper::ftype angleuv, clipper::ftype theta, clipper::ftype frac=1.0/1.17741);
		
		void ellipse(clipper::ftype sig1, clipper::ftype sig2, clipper::ftype cov, 
					 clipper::ftype& angle, clipper::ftype& sigmau, clipper::ftype& sigmav);
	};
		
	//--------------------------------------------------------------
	//! Check for tNCS using patterson map
	/*! Use patterson map to check for tNCA
	 \ingroup g_funcobj */
	
	template<class T> class tNCS
	{
	public:
		//! operator
		const std::vector<clipper::Coord_frac>& operator()(clipper::HKL_data<clipper::datatypes::I_sigI<T> >& I, clipper::Resolution reso=clipper::Resolution(4.0));
		//! has tNCS test
		bool hasNCS() { return peaks.size() != 0; }
		//! number of operators
		int numOps() { return peaks.size(); } 
		//! return peak height
		clipper::ftype height(int i) { return peak_height[i]; }
		//! return peak prob
		clipper::ftype prob(int i) { return peak_prob[i]; }
		//! print out text summary
		void summary();
        //! return tNCS operator as fractional coordinate
        const clipper::Coord_frac& operator[](const int i) { return peaks[i]; }
		
	private:
		clipper::HKL_data<clipper::datatypes::I_sigI<T> >* intensity; //!< dataset
		clipper::Resolution reso; //!< resolution limits
		std::vector<clipper::Coord_frac> peaks; //!< located peaks
		std::vector<T> peak_height; //!< height of located peaks as fraction of origin
		std::vector<T> peak_prob; //!< Zwartz estimation of probability
		
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
	
    //---------Calculate anisotropy-------------------------------------
    
    //! Compute anisotropy correction
    /*! Perform calculation using ML with spherical restraint or hybrid wilson plus correction method
     */
    template <class SCALER, class DATA, class T> class AnisoCorr {
    public:
        //enum TYPE { F, I };
        //! empty constructor
        AnisoCorr() {}
        //! constructor using observation, I_sigI or F_sigF
        AnisoCorr(clipper::HKL_data<DATA>& observed, bool protein=true, bool rna=false) : _observed(&observed), _is_protein(protein), _is_nucl(rna) { calc(observed, protein, rna); }
        //! destructor
        ~AnisoCorr() { }
        //! perform calculation returing the U_aniso_orth
        const clipper::U_aniso_orth& operator()(clipper::HKL_data<DATA>& observed, bool protein=true, bool rna=false);
        //! return the U_aniso_orth
        const clipper::U_aniso_orth& u_aniso_orth( Scaling::TYPE t ) const;
        //! return the scale factor
		const T kscale() const { return _iscale.kscale(); }

    protected:
        //! routines to allow use of F_sigF and I_sigI
        const T&    obs( const clipper::datatypes::F_sigF<T>& f ) { return f.f(); }
        const T&    obs( const clipper::datatypes::I_sigI<T>& f ) { return f.I(); }
        const T& sigobs( const clipper::datatypes::F_sigF<T>& f ) { return f.sigf(); }
        const T& sigobs( const clipper::datatypes::I_sigI<T>& f ) { return f.sigI(); }
        
        //! calculate anisotropy correction, eigenvalues and directional vector
        void calc(clipper::HKL_data<DATA>& observed, bool protein, bool rna);
                  
    private:   
        clipper::HKL_data<DATA>* _observed; //!< pointer for observed data
        bool _is_protein;                //!< cell contains protein
        bool _is_nucl;                    //!< cell contains rna/dna
        SCALER _iscale;                        //!< scaling object
        //clipper::U_aniso_orth _U_f;      //!< computed correction structure factors
        //clipper::U_aniso_orth _U_i;      //!< computed correction intensity
        //clipper::ftype _scale;           //!< computed scale factor
    };
    
    //! Compute eigenvalues and eigenvectors associated with an U_aniso_orth
    /*! Store eigenvalues and eigenvectors for U_ansio_orth decomposition.
     Store in a*, b*, c* order, but also reference max eigenvalue
     */
    template <class T> class AnisoDirection  {
    public:
        //! construct from UAO
        explicit AnisoDirection(clipper::U_aniso_orth& uao);
        //! return eigenvalues, sorted closest to a*,b*,c*
        const std::vector<T>& eigenValues() const 
        { return _eigenvalues; }
        //! return eigenvectors, sorted closest to a*,b*,c*
        const std::vector<clipper::Vec3<T> >& eigenVectors() const 
        { return _eigenvectors; }
        //! return max eigenvalue
         T max() { return _max; }
        
    private:
        clipper::U_aniso_orth* _uao;       //!< reference setup UAO
        std::vector<T> _eigenvalues; //!< sorted eigenvalues
        std::vector<clipper::Vec3<T> > _eigenvectors; //!< sorted eigenvalues
        T _max;               //!< max eigenvalue
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
    const T& sigobs_mi( const clipper::datatypes::J_sigJ_ano<T>& f ) { return f.sigI_mi(); }
    const T&    obs_pl( const clipper::datatypes::G_sigG_ano<T>& f ) { return f.f_pl(); }
    const T&    obs_mi( const clipper::datatypes::G_sigG_ano<T>& f ) { return f.f_mi(); }
    const T& sigobs_pl( const clipper::datatypes::G_sigG_ano<T>& f ) { return f.sigf_pl(); }
    const T& sigobs_mi( const clipper::datatypes::G_sigG_ano<T>& f ) { return f.sigf_mi(); }

private:
    int _nbins;
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
	
}
#endif
