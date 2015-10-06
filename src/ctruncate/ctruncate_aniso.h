//
//     CTRUNCATE
//     Copyright (C) 2006-2011 Norman Stein, Charles Ballard
//
//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Application.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//

#ifndef __CTRUNCATE_ANISO_H
#define __CTRUNCATE_ANISO_H

#include "clipper/clipper.h"
#include "clipper/clipper-ccp4.h"
#include "alt_hkl_datatypes.h"
#include "intensity_scale.h"
#include "ctruncate_wilson.h"
#include "ctruncate_analyse.h"

namespace ctruncate {

	void yorgo_modis_plot(clipper::HKL_data<clipper::data32::F_sigF>& fsig, float maxres, int nbins, CCP4Program& prog, clipper::U_aniso_orth uao = clipper::U_aniso_orth(1.0));
	
	//void yorgo_modis_plot(clipper::HKL_data<clipper::data32::I_sigI>& isig, float maxres, int nbins, CCP4Program& prog);
	
	void yorgo_modis_plot(clipper::HKL_data<clipper::data32::I_sigI>& isig, float maxres, int nbins, CCP4Program& prog, clipper::U_aniso_orth uao = clipper::U_aniso_orth(1.0) );
	
	//--------------------------------------------------------------
	
	class YorgoModis
	{
	public:
		YorgoModis() { }
		template<class T, template<class> class D> YorgoModis(const clipper::HKL_data<D<T> >&, clipper::U_aniso_orth uao = clipper::U_aniso_orth(1.0));
		template<class T, template<class> class D, class TT> YorgoModis(const clipper::HKL_data<D<T> >&, const clipper::Range<TT>&, clipper::U_aniso_orth uao = clipper::U_aniso_orth(1.0));
		~YorgoModis() { }
		template<class T, template<class> class D> void operator() (const clipper::HKL_data<D<T> >&, clipper::U_aniso_orth uao = clipper::U_aniso_orth(1.0));
		template<class T, template<class> class D, class TT> void operator() (const clipper::HKL_data<D<T> >&, const clipper::Range<TT>&, clipper::U_aniso_orth uao = clipper::U_aniso_orth(1.0));
		void output() const;
		std::stringstream& xml_output(std::stringstream&) const;
		bool is_intensity() const { return _base->type() == "I_sigI_ano" || _base->type() == "J_sigJ_ano" || _base->type() == "I_sigI"; }
		bool is_anomalous() const { return _base->type() == "I_sigI_ano" || _base->type() == "J_sigJ_ano" || _base-> type() == "F_sigF_ano"; }
		int size() const { return _b_reso.size(); }
		
	private:
		template<class T> const T&    obs( const clipper::datatypes::F_sigF<T>& f ) { return f.f(); }
		template<class T> const T&    obs( const clipper::datatypes::I_sigI<T>& f ) { return f.I(); }
		template<class T> const T&    obs( const clipper::datatypes::F_sigF_ano<T>& f ) { return f.f(); }
		template<class T> const T&    obs( const clipper::datatypes::I_sigI_ano<T>& f ) { return f.I(); }
		template<class T> const T& sigobs( const clipper::datatypes::F_sigF<T>& f ) { return f.sigf(); }
		template<class T> const T& sigobs( const clipper::datatypes::I_sigI<T>& f ) { return f.sigI(); }
		template<class T> const T& sigobs( const clipper::datatypes::F_sigF_ano<T>& f ) { return f.sigf(); }
		template<class T> const T& sigobs( const clipper::datatypes::I_sigI_ano<T>& f ) { return f.sigI(); }
		
		const clipper::HKL_data_base *_base;
		clipper::Range<clipper::ftype> _reso;       //resolution range in calculation
		clipper::U_aniso_orth _uao;   //!<anisotropy matrix for directions
		
		clipper::Mat33<clipper::ftype> _e123; //!<directions of eigevectors
		
        ctruncate::ResolutionRange_ordinal _s_ord;    //<! ordinal
		std::vector<clipper::ftype> _b_reso;  //!< resolution bins
		std::vector<clipper::ftype> _somov;
		std::vector<clipper::ftype> _somsdov;
		std::vector<clipper::ftype> _numov;
		std::vector<clipper::ftype> _somdir;    
		std::vector<clipper::ftype> _somsddir;
		std::vector<clipper::ftype> _numdir;
        std::vector<clipper::Range<clipper::ftype> > _b_range; //< ! reso range of bin
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

    //---------Calculate anisotropy-------------------------------------
    
    //! Compute anisotropy correction
    /*! Perform calculation using ML with spherical restraint or hybrid wilson plus correction method
     */
    template <class SCALER, class DATA, class T> class AnisoCorr {
    public:
        //enum TYPE { F, I };
        //! empty constructor
        AnisoCorr() : _is_protein(false), _is_nucl(false) {}
		//! set reso range
		explicit AnisoCorr(clipper::Range<clipper::ftype>& reso) : _is_protein(false), _is_nucl(false), _range(reso) {}
		//! set reso range
		explicit AnisoCorr(bool protein, bool rna) : _is_protein(protein), _is_nucl(rna) {}
        //! constructor using observation, I_sigI or F_sigF
        AnisoCorr(clipper::HKL_data<DATA>& observed, bool protein=false, bool rna=false, clipper::Range<clipper::ftype> reso = clipper::Range<clipper::ftype>() ) 
		: _observed(&observed), _is_protein(protein), _is_nucl(rna), _range(reso) { calc(observed); }
        //! destructor
        ~AnisoCorr() { }
        //! perform calculation returing the U_aniso_orth
        const clipper::U_aniso_orth& operator()(clipper::HKL_data<DATA>& observed, bool protein=true, bool rna=false, 
												clipper::Range<clipper::ftype> reso = clipper::Range<clipper::ftype>());
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
        void calc(clipper::HKL_data<DATA>& observed);
		
    private:   
        clipper::HKL_data<DATA>* _observed; //!< pointer for observed data
        bool _is_protein;                //!< cell contains protein
        bool _is_nucl;                    //!< cell contains rna/dna
		clipper::Range<clipper::ftype> _range; //!< active reso range
        SCALER _iscale;                        //!< scaling object
        //clipper::U_aniso_orth _U_f;      //!< computed correction structure factors
        //clipper::U_aniso_orth _U_i;      //!< computed correction intensity
        //clipper::ftype _scale;           //!< computed scale factor
    };
    
    //! Compute eigenvalues and eigenvectors associated with an U_aniso_orth
    /*! Store eigenvalues and eigenvectors for U_ansio_orth decomposition.
     Store in a*, b*, c* order, but also reference max eigenvalue
     */
    class AnisoDirection  {
    public:
		AnisoDirection() {}
		~AnisoDirection() {}
        //! construct from UAO
        explicit AnisoDirection(const clipper::U_aniso_orth& uao);
		//! compute eigenvectors and eigenvalues
		const std::vector<clipper::ftype>& operator()(const clipper::U_aniso_orth& uao);
        //! return eigenvalues, sorted closest to a*,b*,c*
        const std::vector<clipper::ftype>& eigenValues() const 
        { return _eigenvalues; }
        //! return eigenvectors, sorted closest to a*,b*,c*
        const std::vector<clipper::Vec3<clipper::ftype> >& eigenVectors() const 
        { return _eigenvectors; }
        //! return max eigenvalue
		clipper::ftype max() const { return _max; }
        
    private:
        const clipper::U_aniso_orth* _uao;       //!< reference setup UAO
        std::vector<clipper::ftype> _eigenvalues; //!< sorted eigenvalues
        std::vector<clipper::Vec3<clipper::ftype> > _eigenvectors; //!< sorted eigenvalues
		clipper::ftype _max;               //!< max eigenvalue
    };
	
	//AnisoAnalysis--------------------------------------------------
	//! Wrapper for Anisotropy of data analysis
	/*! Wrap AnisoCorr, AnisoDirection and YorgoModis
	 */
	class AnisoAnalysis {
	public:
		//! constructor
		AnisoAnalysis() { }
		template<class T, template<class> class D> AnisoAnalysis(const clipper::HKL_data< D<T> >&);
		template<class T, template<class> class D, class TT> AnisoAnalysis(const clipper::HKL_data< D<T> >&, const clipper::Range<TT>&);
		//! destructor
		~AnisoAnalysis() { }
		//! anisotropy allowed by symmetry restraints
		bool allowed_by_symmetry() const { return _anisobysymm; }
		static bool allowed_by_symmetry(const clipper::Spacegroup&, const clipper::Cell&);
		//! return anisotropic temp
		const clipper::U_aniso_orth& u_aniso_orth_F() const { return _u_aniso_orth_f; }
		const clipper::U_aniso_orth& u_aniso_orth_I() const { return _u_aniso_orth_i; }
        //! return anisotropy correction
        const clipper::U_aniso_orth& u_aniso_orth_corr_F() const { return _u_aniso_orth_corr_f; }
		const clipper::U_aniso_orth& u_aniso_orth_corr_I() const { return _u_aniso_orth_corr_i; }
		//is anisotropic
		bool is_anisotropic() const;
		//output to std::out
		void output() const;
		//return xml
		std::stringstream& xml_output(std::stringstream&) const;
		//return data type (intensity or amplitudes)
		bool is_intensity() const { return _base->type() == "I_sigI_ano" || _base->type() == "J_sigJ_ano" || _base->type() == "I_sigI"; }
		bool is_anomalous() const { return _base->type() == "I_sigI_ano" || _base->type() == "J_sigJ_ano" || _base-> type() == "F_sigF_ano"; }
	private:
		const clipper::HKL_data_base *_base; //!< pointer to parent data
		bool _anisobysymm; //!< cache response of allowed_by_symm
		clipper::U_aniso_orth _u_aniso_orth_corr_f; //!< computed anisotropy correction for amplitudes
		clipper::U_aniso_orth _u_aniso_orth_corr_i; //!< computed anisotropy correction for intensities
		clipper::U_aniso_orth _u_aniso_orth_f; //!< computed anisotropic U for amplitudes
		clipper::U_aniso_orth _u_aniso_orth_i; //!< computed anisotropic U for intensities
		clipper::Range<clipper::ftype> _range; //!< range used in calculations
		
		AnisoDirection _dm; //!<
		YorgoModis _ym; //!< Yorgo Modis plot
		
		template<class T, template<class> class D, class TT> void init(const clipper::HKL_data< D<T> >&, const clipper::Range<TT>&);
	};
	
	//definitions-----------------------------------------------------

	template<class T, template<class> class D> AnisoAnalysis::AnisoAnalysis(const clipper::HKL_data<D<T> >& hkldata)
	{
		clipper::Range<clipper::ftype> range(hkldata.invresolsq_range() );
		init(hkldata,range);
	}
	
	template<class T, template<class> class D, class TT> AnisoAnalysis::AnisoAnalysis(const clipper::HKL_data<D<T> >& hkldata, const clipper::Range<TT>& range)
	{
		init(hkldata,range);
	}
	
	template<class T, template<class> class D, class TT> void AnisoAnalysis::init(const clipper::HKL_data<D<T> >& hkldata, const clipper::Range<TT>& range)
	{
		_anisobysymm = allowed_by_symmetry(hkldata.hkl_info().spacegroup(), hkldata.hkl_info().cell() );
		_base = &hkldata;
		_range = ((range.min() > range.max() ) ? hkldata.invresolsq_range() : range );
		
		if (_anisobysymm) {
            clipper::HKL_data<D<T> > ianiso(hkldata.hkl_info() );
            
            for ( clipper::HKL_data_base::HKL_reference_index ih = ianiso.first(); !ih.last(); ih.next() ) {
                if (_range.contains(ih.invresolsq() ) ) {
                    ianiso[ih] = D<T>( hkldata[ih.hkl() ] );
                }
            }
            
            //corrections
			try { 
				AnisoCorr<Iscale_logLikeAniso<T>, D<T>, T > llscl(ianiso, false, false, _range);
				_u_aniso_orth_corr_i = -(llscl.u_aniso_orth(Scaling::I) );
				_u_aniso_orth_corr_f = -(llscl.u_aniso_orth(Scaling::F) );
			} catch (clipper::Message_fatal) {
                //clipper::Message_(1, "Anisotropy anlysis failed.");
                throw;
			}
			
            //reset for full temperature term
            for ( clipper::HKL_data_base::HKL_reference_index ih = ianiso.first(); !ih.last(); ih.next() ) {
                if (_range.contains(ih.invresolsq() ) ) {
                    ianiso[ih] = D<T>( hkldata[ih.hkl() ] );
                }
            }
            try {
				AnisoCorr<Iscale_logLikeAniso<T>, D<T>, T > llscl(ianiso, true, false, _range);
				_u_aniso_orth_i = (llscl.u_aniso_orth(Scaling::I) );
				_u_aniso_orth_f = (llscl.u_aniso_orth(Scaling::F) );
			} catch (clipper::Message_fatal) {
                //clipper::Message_(1, "Anisotropy anlysis failed.");
                throw;
			}
            
			_dm(_u_aniso_orth_i);

			//YorgoModis plot
			_ym(hkldata,_u_aniso_orth_i);
		}
	}
	
	template<class T, template<class> class D, class TT> YorgoModis::YorgoModis(const clipper::HKL_data<D<T> >& hkldata, const clipper::Range<TT>& range,clipper::U_aniso_orth uao)
	{
		this->operator()(hkldata ,range, uao );
	}
	
	template<class T, template<class> class D> YorgoModis::YorgoModis(const clipper::HKL_data<D<T> >& hkldata,clipper::U_aniso_orth uao)
	{
		clipper::Range<clipper::ftype> range;
		this->operator()(hkldata ,range, uao );
	}
	
	template<class T, template<class> class D> void YorgoModis::operator()(const clipper::HKL_data<D<T> >& hkldata, clipper::U_aniso_orth uao)
	{
		clipper::Range<clipper::ftype> range;
		this->operator()(hkldata ,range, uao );
	}
	
	template<class T, template<class> class D, class TT> void YorgoModis::operator()(const clipper::HKL_data<D<T> >& hkldata, const clipper::Range<TT>& range, clipper::U_aniso_orth uao )
	{
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		_base = &hkldata;
		_reso = ((range.min() > range.max() ) ? hkldata.invresolsq_range() : range );
		_uao = uao;
		
        //get uao eigenvectors
		AnisoDirection direct(uao);
		for (int i = 0 ; i !=3 ; ++i)
			for (int j = 0 ; j != 3 ; ++j)
                _e123(i,j) = (direct.eigenVectors())[i][j];
		
		clipper::Cell cell = hkldata.hkl_info().cell();
		clipper::Spacegroup spg = hkldata.hkl_info().spacegroup();
		
		int _nzerosigma = 0;
		float cone = 30.0; //hardwired for now
		float ang;
		float cosang;
		
		//ctruncate::ResolutionRange_ordinal s_ord;
		_s_ord.init( hkldata.hkl_info(), _reso, 1.0 );
		
		int nbins(0);
		int Nreflections(0);
		for ( HRI ih = hkldata.first(); !ih.last(); ih.next() ) {
			if (_reso.contains(ih.invresolsq() ) ) ++Nreflections;
		}
		
		int nspg = spg.num_primitive_symops();
		const int nreflns(1000);
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
				
		_b_reso.resize(nbins,0.0);
		_somov.resize(nbins,0.0); _somsdov.resize(nbins,0.0); _numov.resize(nbins,0);
		_somdir.resize(3*nbins,0.0); _somsddir.resize(3*nbins,0.0); _numdir.resize(3*nbins,0.0);
        _b_range.resize(nbins);
		
		//std::fill(_b_reso.begin(), _b_reso.end(), 0.0); std::fill(_somov.begin(), _somov.end(), 0.0); std::fill(_somsdov.begin(), _somsdov.end(), 0.0);
		//std::fill(_numov.begin(), _numov.end(), 0.0); std::fill(_somdir.begin(), _somdir.end(), 0.0); std::fill(_somsddir.begin(), _somsddir.end(), 0.0);
		//std::fill(_numdir.begin(), _numdir.end(), 0.0);
		
		std::vector<clipper::ftype> nr(nbins,0.0);
		for ( HRI ih = hkldata.first(); !ih.last(); ih.next() ) {
			clipper::ftype s = ih.invresolsq();
			if (_reso.contains(s) ) {
				clipper::ftype epsiln = (is_intensity() ) ? 1.0f/ih.hkl_class().epsilonc() : std::sqrt(1.0f/ih.hkl_class().epsilonc() );
				int bin = clipper::Util::bound( 0,clipper::Util::intf( clipper::ftype(nbins) * _s_ord.ordinal( s ) ), nbins-1 );
				_b_reso[bin] += nspg*epsiln*s;
				nr[bin] += nspg*epsiln;
                _b_range[bin].include(s);
				if ( !hkldata[ih].missing() ) {
					
					clipper::ftype val=obs(hkldata[ih])*epsiln;
					clipper::ftype sig=sigobs(hkldata[ih])*epsiln;
					if ( sig > 0.0 ) {
						for ( int jsym = 0; jsym != nspg ; ++jsym ) {
							for (int friedal = 0 ; friedal != 2 ; ++friedal) {
								clipper::HKL ri = int(std::pow( -1.0f, float(friedal) ))*ih.hkl();
								clipper::HKL rj = ri.transform( spg.primitive_symop( jsym ) );
								
								clipper::Vec3<clipper::ftype> hc = _e123*clipper::Vec3<clipper::ftype>(rj.coord_reci_orth(cell) );  //transpose into eigenspace
								
								_somov[bin] += val;
								_somsdov[bin] += val/sig;
								_numov[bin] += 1.0;
								for (int j=0;j!=3;++j) {
									int jn = j*nbins+bin;
									cosang = fabs( hc[j] )/sqrt(ih.invresolsq());
									// cosang can stray just past 1.0
									cosang = std::min(cosang, 1.0f);
									ang = acos(cosang);
									if ( ang < clipper::Util::d2rad(cone) ) {
										_somdir[jn] += val;
										_somsddir[jn] += val/sig;
										_numdir[jn] += 1.0;
									}
								}
							}
						}
					}
				}
			} 
		}
		
		for (int i=0;i != nbins; ++i) {
			for (int j=0;j!=3;++j) {
				int jn = j*nbins+i;
				if (_numdir[jn] < 0.5) {
					_somdir[jn] = 0.0;
					_somsddir[jn] = 0.0;
				} else {
					_somdir[jn] /= _numdir[jn];
					_somsddir[jn] /= _numdir[jn];
				}
			}
			if (_numov[i] < 0.5) {
				_somov[i] = 0.0;
				_somsdov[i] = 0.0;
			} else {
				_somov[i] /= _numov[i];
				_somsdov[i] /= _numov[i];
			}
			_b_reso[i] /= nr[i];
		}
		return;
	}
    
    //Directional NonAnomAnalysis-------------------------------------------------
    //! analyse data for directional completeness
    class HKLStats_d_completeness : public ResolStats_base {
    public:
        HKLStats_d_completeness(const clipper::HKL_data_base& hkldata, clipper::U_aniso_orth& uao, clipper::ftype val=-99.0);
        HKLStats_d_completeness(const ResolStats_base&, clipper::U_aniso_orth& uao, clipper::ftype val=-99.0);
        HKLStats_d_completeness(clipper::ftype val=-99.0) { _val=val; }
        ~HKLStats_d_completeness() {}
        clipper::ftype operator[](int index) const;
        clipper::ftype operator()(const int, const int) const;
        clipper::ftype IoversigI() { return _val; }
        clipper::Vec3<clipper::ftype> direction(const int) const;
        HKLStats_d_completeness& operator=(const HKLStats_d_completeness&);
    private:
        clipper::Mat33<clipper::ftype> _e123; //!<directions of eigevectors
        std::vector<clipper::ftype> _d1_completeness;  //!< completeness stats
        std::vector<clipper::ftype> _d2_completeness;  //!< completeness stats
        std::vector<clipper::ftype> _d3_completeness;  //!< completeness stats
        clipper::ftype _val;                //!< reference value
        
        void calc(const clipper::HKL_data_base& hkldata);
    };

	
}


#endif
