
#ifndef CTRUNCATE_ISCALE
#define CTRUNCATE_ISCALE


#include "clipper/contrib/function_object_bases.h"
#include "clipper/core/resol_targetfn.h"


namespace ctruncate {
	
    class Scaling {
    public:
        enum TYPE { F, I };
        enum MODE { NORMAL, SHARPEN, UNSHARPEN };
    };
    
	template<class T> class Iscale_aniso_base : public clipper::SFscale_base<T> {
	public:
		//enum TYPE { F, I };  //!< type for returning U_aniso_orth
		//enum MODE { NORMAL, SHARPEN, UNSHARPEN };  //!< mode for scaling
        virtual bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& ) = 0;
        virtual bool operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& ) = 0;
        virtual bool operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >&, 
								 const clipper::HKL_data<clipper::datatypes::I_sigI<T> >&,
								 clipper::Range<clipper::ftype>) = 0;
        virtual bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >&,
								 const clipper::HKL_data<clipper::datatypes::F_sigF<T> >&,
								 clipper::Range<clipper::ftype> ) = 0;
		const clipper::U_aniso_orth& u_aniso_orth( Scaling::TYPE t ) const;
		const clipper::U_aniso_orth& u_aniso_orth() const { return u_i; }  //!< \deprecated
		const clipper::ftype kscale() const { return bscale_; }
		//protected:
		const T&    obs( const clipper::datatypes::F_sigF<T>& f ) { return f.f(); }
		const T&    obs( const clipper::datatypes::I_sigI<T>& f ) { return f.I(); }
		const T& sigobs( const clipper::datatypes::F_sigF<T>& f ) { return f.sigf(); }
		const T& sigobs( const clipper::datatypes::I_sigI<T>& f ) { return f.sigI(); }
		clipper::U_aniso_orth u_i, u_f;
		clipper::ftype bscale_;
		clipper::ftype nsig_;
        Scaling::MODE mode_;
	};
	
	
    //! abstract base class for restraints in fit based on global parameters
    /*! return addition to gradients and hessian based on restraints wrt parameters
	 */
    
    class RestraintFn_base
    {
    public:
        class Aderiv
        {
        public:
            clipper::ftype f;
            std::vector<clipper::ftype> df;
            clipper::Matrix<> df2;
            Aderiv() {}
            Aderiv(const int& np) : df(np,0.0), df2(np,np,0.0) {}
        };
        //! null constructo
        RestraintFn_base() {}
        //! constructor with number of parameters
        RestraintFn_base( const int np) : np_(np), aderiv_(np) {}
        //! the value of the global restraint function
        virtual clipper::ftype f( const clipper::Cell& cell, const std::vector<clipper::ftype>& params ) const { return aderiv( cell, params ).f; }
        //! compute the adustment due to the restraint
        virtual const Aderiv& aderiv(const clipper::Cell& cell, const std::vector<clipper::ftype>&  params) const = 0;
        //! provide result for derived classes
        Aderiv& result() const { return aderiv_; }
        //! return number of parameters
        const int& num_params() const { return np_; }
        //! destructor
        virtual ~RestraintFn_base() {}
        
    private:
        int np_;  //!< number of parameters
        mutable Aderiv aderiv_; //<!store state
    };
    
    //! spherical restraint on atomic displacement parameter
    /*! This restrains the the parameter to be spherical
     
     Murshudov et al. Acta D55 (1999) 247
     
	 d^2 = sum(uii-Uiso)^2 + sum(uij^2)
     
     used in exponental form of PDF
     */
    class RestraintFn_sphericalU : public RestraintFn_base
    {
    public:
        //! constructor taking weight for restraint
        RestraintFn_sphericalU(clipper::ftype f=10.0, clipper::ftype n=1.0) : sigma_(f), nobs_(n), RestraintFn_base(7) { }
        //! the value of the global restraint function
        clipper::ftype f( const clipper::Cell& cell, const std::vector<clipper::ftype>& params ) const { return f_s(params); }
        //! compute value of restaint without doing full aderiv
        clipper::ftype f_s( const std::vector<clipper::ftype>& params ) const;
        //! return adjustments to derivatives and hessian
        const Aderiv& aderiv(const clipper::Cell& cell, const std::vector<clipper::ftype>&  params) const;
        
    private:
        clipper::ftype sigma_; //!< scale factor for restraint
        clipper::ftype nobs_;  //!< number of observations to scale sigma
		
    };
    
	//! 2nd order resolution function evaluator
    /*! This is an automatic evaluator for arbitrary functions of HKL,
     most commonly used for evaluating a function of resolution (such a
     mean F^2 or sigmaa), although more general tasks including local
     scaling of reflections and anisotropic functions can also be
     handled. This form is for target functions which approach zero
     quadratically, e.g. least-squares targets.
     
     Symmetry is used as an additional constraint.
     
     Allow additional restraints function.
     
     \note This version implements a minimiser which uses both
     Newton-Raphson and gradient steps depending on the situation. It
     can be used for non-quadratic targets or non-linear basis
     functions.
     
     To evaluate a resolution function, this class must be provided
     with two objects:
     - The basis function (and gradients), which describes the value
     of the function for any reflection given a set of paramters.
     - The target function (and derivatives), which is used to determine
     the values of the basis function parameters.
     */
    class ResolutionFn_nonlinear : public clipper::ResolutionFn
    {
    public:
        //! constructor: need reflections, basis fn, target fn and restraints. No symops.
        ResolutionFn_nonlinear( const clipper::HKL_info& hkl_info, 
                                    const clipper::BasisFn_base& basisfn, 
                                    const clipper::TargetFn_base& targetfn, 
                                    const std::vector<clipper::ftype>& params, 
			        	const std::vector<bool>& mask,
                                    const clipper::ftype damp = 0.0, 
                                    const bool debug = false ) throw(clipper::Message_fatal);
		
    };
	
    //! 2nd order resolution function evaluator
    /*! This is an automatic evaluator for arbitrary functions of HKL,
     most commonly used for evaluating a function of resolution (such a
     mean F^2 or sigmaa), although more general tasks including local
     scaling of reflections and anisotropic functions can also be
     handled. This form is for target functions which approach zero
     quadratically, e.g. least-squares targets.
     
     Symmetry is used as an additional constraint.
     
     Allow additional restraints function.
     
     \note This version implements a minimiser which uses both
     Newton-Raphson and gradient steps depending on the situation. It
     can be used for non-quadratic targets or non-linear basis
     functions.
     
     To evaluate a resolution function, this class must be provided
     with two objects:
     - The basis function (and gradients), which describes the value
     of the function for any reflection given a set of paramters.
     - The target function (and derivatives), which is used to determine
     the values of the basis function parameters.
     */
    class ResolutionFn_nonlinear_rest : public clipper::ResolutionFn
    {
    public:
        //! constructor: need reflections, basis fn, target fn and restraints. No symops.
        ResolutionFn_nonlinear_rest( const clipper::HKL_info& hkl_info, 
                                    const clipper::BasisFn_base& basisfn, 
                                    const clipper::TargetFn_base& targetfn, 
                                    const std::vector<clipper::ftype>& params, 
					const std::vector<bool>& mask,
                                    const ctruncate::RestraintFn_base& restraint,
                                    const clipper::ftype damp = 0.0, 
                                    const bool debug = false ) throw(clipper::Message_fatal);
		
    protected:
        const RestraintFn_base* restfn_;   //!< retraints function        
		
    };
    
	
	//! Structure factor anisotropic scaling
	/*! Perform structure factor anisotropic scaling, observed to calculated,
	 calculated to observed, or observed against itself using log target
	 \ingroup g_funcobj */
	template<class T> class Iscale_logAniso : public ctruncate::Iscale_aniso_base<T> {
	public:
		Iscale_logAniso( clipper::ftype nsig = 3.0, Scaling::MODE mode = Scaling::NORMAL )
        { Iscale_aniso_base<T>::nsig_ = nsig; Iscale_aniso_base<T>::mode_ = mode; }
		//! Scale Fc to Fo
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_phi<T> >& fo, 
						 const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fc );
		//! Scale Fc to Fo
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_phi<T> >& fo, 
						 const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fc,
						 clipper::Range<clipper::ftype> resfilter );
        //! Scale Fo to Fc
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, 
						 const clipper::HKL_data<clipper::datatypes::F_phi<T> >& fc );
		//! Scale Fo to Fc
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, 
						 const clipper::HKL_data<clipper::datatypes::F_phi<T> >& fc,
						 clipper::Range<clipper::ftype> resfilter );		
		//! Scale Fc to Fo
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fc, 
						 const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo,
						 clipper::Range<clipper::ftype> resfilter=clipper::Range<clipper::ftype>() );
		//! Scale Fo to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo );
		//! Scale Fo to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, 
						 clipper::Range<clipper::ftype>& resfilter,
						 const int npar_scl );
		//! Scale Io to Ic (psuedo I), or vice versa
		bool operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io, 
						 const clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Ic,
						 clipper::Range<clipper::ftype> resfilter=clipper::Range<clipper::ftype>() );
		//! Scale Io to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io );
		//! Scale Io to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io, 
						 clipper::Range<clipper::ftype>& resfilter, const int npar_scl );
		//! Primitive scaling functions for F and I
		template<class D, class T1, class T2, class S>
		bool scale( clipper::HKL_data<D>& fo, clipper::Range<clipper::ftype>& resfilter, const int npar_scl );
	};
	
    //! Structure factor anisotropic scaling
	/*! Perform structure factor anisotropic scaling, observed to calculated,
	 calculated to observed, or observed against itself using standard target
	 \ingroup g_funcobj */
	template<class T> class Iscale_aniso : public Iscale_aniso_base<T> {
	public:
		Iscale_aniso( clipper::ftype nsig = 3.0, Scaling::MODE mode = Scaling::NORMAL )
		{ Iscale_aniso_base<T>::nsig_=nsig; Iscale_aniso_base<T>::mode_=mode; }
		//! Scale Fo to Fc
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_phi<T> >& fo, 
						 const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fc );
		//! Scale Fc to Fo
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_phi<T> >& fo, 
						 const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fc,
						 clipper::Range<clipper::ftype> resfilter );
        //! Scale Fo to Fc
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, 
						 const clipper::HKL_data<clipper::datatypes::F_phi<T> >& fc );
		//! Scale Fo to Fc
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, 
						 const clipper::HKL_data<clipper::datatypes::F_phi<T> >& fc,
						 clipper::Range<clipper::ftype> resfilter );		
		//! Scale Fc to Fo
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fc, 
						 const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo,
						 clipper::Range<clipper::ftype> resfilter=clipper::Range<clipper::ftype>() );
		//! Scale Fo to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo );
		//! Scale Fo to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, 
						 clipper::Range<clipper::ftype> resfilter,
						 const int npar_scl );
		//! Scale Io to Ic (psuedo I), or vice versa
		bool operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io, 
						 const clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Ic,
						 clipper::Range<clipper::ftype> resfilter=clipper::Range<clipper::ftype>() );
		//! Scale Io to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io );
		//! Scale Io to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io, 
						 clipper::Range<clipper::ftype>,
						 const int npar_scl );
		//! Primitive scaling functions for F and I
		template<class D, class T1, class T2, class S>
		bool scale( clipper::HKL_data<D>& fo, clipper::Range<clipper::ftype>& resfilter, const int npar_scl );
		//! return aniso correction on F or I
	};

    //! Structure factor anisotropic scaling
	/*! Perform structure factor anisotropic scaling, observed to calculated,
	 calculated to observed, or observed against itself using loglikelihood target
	 \ingroup g_funcobj */
	template<class T> class Iscale_logLikeAniso : public ctruncate::Iscale_aniso_base<T> {
	public:
		Iscale_logLikeAniso( clipper::ftype nsig = 0.0, Scaling::MODE mode = Scaling::NORMAL )
        { Iscale_aniso_base<T>::nsig_ = nsig; Iscale_aniso_base<T>::mode_ = mode; }
		//! Scale Fc to Fo
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_phi<T> >& fo, 
						 const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fc );
		//! Scale Fc to Fo with resolution range
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_phi<T> >& fo, 
						 const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fc,
						 clipper::Range<clipper::ftype> resfilter );
		//! Scale Fc to Fo with resolution range and mask
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_phi<T> >& fo, 
						 const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fc,
						 clipper::Range<clipper::ftype> resfilter, std::vector<bool> mask );		
        //! Scale Fo to Fc
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, 
						 const clipper::HKL_data<clipper::datatypes::F_phi<T> >& fc );
		//! Scale Fo to Fc with resolution range
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, 
						 const clipper::HKL_data<clipper::datatypes::F_phi<T> >& fc,
						 clipper::Range<clipper::ftype> resfilter );
		//! Scale Fo to Fc with resolution range
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, 
						 const clipper::HKL_data<clipper::datatypes::F_phi<T> >& fc,
						 clipper::Range<clipper::ftype> resfilter, std::vector<bool> mask );
		//! Scale Fo to Fo with resolution range (optional)
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fc, 
						 const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo,
						 clipper::Range<clipper::ftype> resfilter=clipper::Range<clipper::ftype>() );
		//! Scale Fo to Fo with resolution range and mask
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fc, 
						 const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo,
						 clipper::Range<clipper::ftype> resfilter, std::vector<bool> mask );
		//! Scale Fo to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo );
		//! Scale Fo to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, 
						 clipper::Range<clipper::ftype> resfilter,
						 const int npar_scl );
		//! Scale Io to Ic (psuedo I), or vice versa, with resolution range (optional)
		bool operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io, 
						 const clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Ic,
						 clipper::Range<clipper::ftype> resfilter=clipper::Range<clipper::ftype>() );
		//! Scale Io to Ic (psuedo I), or vice versa, with resolution range and mask
		bool operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io, 
						 const clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Ic,
						 clipper::Range<clipper::ftype> resfilter, std::vector<bool> mask );
		//! Scale Io to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io );
		//! Scale Io to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io, 
						 clipper::Range<clipper::ftype> resfilter,
						 const int npar_scl );
		//! Primitive scaling functions for F and I
		template<class D, class T1, class T2, class R2, class S>
		bool scale( clipper::HKL_data<D>& fo, clipper::Range<clipper::ftype>& resfilter, const int npar_scl );
        template<class D1, class D2, class S>
        bool prescale( const clipper::HKL_data<D1>& fo, const clipper::HKL_data<D2>& fc, clipper::ftype& A, clipper::ftype& B, clipper::Range<clipper::ftype>, const int npar_scl=100 );
	};
    
    //! Structure factor anisotropic scaling
	/*! Perform structure factor anisotropic scaling, observed to calculated,
	 calculated to observed, or observed against itself using standard target.
     Use Iscale_aniso plus wilson scaling to get approximation to anisotropic scaling.
	 \ingroup g_funcobj */
	template<class T> class Iscale_wilsonAniso : public Iscale_aniso_base<T> {
	public:
		Iscale_wilsonAniso( clipper::ftype nsig = 0.0, Scaling::MODE mode = Scaling::NORMAL )
		{ Iscale_aniso_base<T>::nsig_=nsig; Iscale_aniso_base<T>::mode_=mode; }
		//! Scale Fc to Fo
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_phi<T> >& fo, 
						 const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fc );
		//! Scale Fc to Fo
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_phi<T> >& fo, 
						 const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fc,
						 clipper::Range<clipper::ftype> resfilter );
        //! Scale Fo to Fc
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, 
						 const clipper::HKL_data<clipper::datatypes::F_phi<T> >& fc );
		//! Scale Fo to Fc
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, 
						 const clipper::HKL_data<clipper::datatypes::F_phi<T> >& fc,
						 clipper::Range<clipper::ftype> resfilter);		
        //! Scale Fo to Fc
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, 
						 const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fc,
						 clipper::Range<clipper::ftype> resfilter=clipper::Range<clipper::ftype>());
		//! Scale Fo to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo );
		//! Scale Fo to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, 
						 clipper::Range<clipper::ftype>& resfilter, const int npar_scl );
		//! Scale Io to Ic (psuedo I), or vice versa
		bool operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io, 
						 const clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Ic,
						 clipper::Range<clipper::ftype> resfilter=clipper::Range<clipper::ftype>() );
		//! Scale Io to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io );
		//! Scale Io to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io, 
						 clipper::Range<clipper::ftype>& resfilter, const int npar_scl );
		//! Primitive scaling functions for F and I
		template<class D, class T1, class T2, class S>
		bool scale( clipper::HKL_data<D>& fo, clipper::Range<clipper::ftype>& resfilter, const int npar_scl );
		//! return aniso correction on F or I
	};
	
	//! Structure factor isotropic scaling
	/*! Perform structure factor isotropic scaling, observed to calculated,
	 calculated to observed, or observed against itself using loglikelihood target
	 \ingroup g_funcobj */
	template<class T> class Iscale_logLikeIso {
	public:
		Iscale_logLikeIso( clipper::ftype nsig = 0.0, Scaling::MODE mode = Scaling::NORMAL )
        { nsig_ = nsig; mode_ = mode; }
		//! Scale Fc to Fo
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_phi<T> >& fo, 
						 const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fc, 
						 clipper::Range<clipper::ftype> resfilter=clipper::Range<clipper::ftype>(),
						 std::vector<bool> mask=std::vector<bool>(2,false) );
        //! Scale Fo to Fc
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, 
						 const clipper::HKL_data<clipper::datatypes::F_phi<T> >& fc,
						 clipper::Range<clipper::ftype> resfilter=clipper::Range<clipper::ftype>(),
						 std::vector<bool> mask=std::vector<bool>(2,false) );
		//! Scale Fo to Fo
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fc, 
						 const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo,
						 clipper::Range<clipper::ftype> resfilter=clipper::Range<clipper::ftype>(),
						 std::vector<bool> mask=std::vector<bool>(2,false) );
		//! Scale Io to Ic (psuedo I), or vice versa
		bool operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io, 
						 const clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Ic,
						 clipper::Range<clipper::ftype> resfilter=clipper::Range<clipper::ftype>(),
						 std::vector<bool> mask=std::vector<bool>(2,false) );
        template<class D1, class D2, class S>
        bool prescale( const clipper::HKL_data<D1>& fo, const clipper::HKL_data<D2>& fc, 
					  clipper::ftype& A, clipper::ftype& B, 
					  clipper::Range<clipper::ftype>& resfilter, 
					  const int npar_scl=60 );
		
		const clipper::U_aniso_orth& u_aniso_orth( Scaling::TYPE t ) const;
		const clipper::U_aniso_orth& u_aniso_orth() const { return u_i; }  //!< \deprecated
		const clipper::ftype kscale() const { return bscale_; }
		
	protected:
		const T&    obs( const clipper::datatypes::F_sigF<T>& f ) { return f.f(); }
		const T&    obs( const clipper::datatypes::I_sigI<T>& f ) { return f.I(); }
		const T& sigobs( const clipper::datatypes::F_sigF<T>& f ) { return f.sigf(); }
		const T& sigobs( const clipper::datatypes::I_sigI<T>& f ) { return f.sigI(); }
		clipper::U_aniso_orth u_i, u_f;
		clipper::ftype bscale_;
		clipper::ftype nsig_;
        Scaling::MODE mode_;
	};
	
} // namespace ctruncate

#endif
