/*!  \file sfscale.h
  Header file for structure factor anisotropic scaling object
*/
//C Copyright (C) 2000-2006 Kevin Cowtan and University of York
//L
//L  This library is free software and is distributed under the terms
//L  and conditions of version 2.1 of the GNU Lesser General Public
//L  Licence (LGPL) with the following additional clause:
//L
//L     `You may also combine or link a "work that uses the Library" to
//L     produce a work containing portions of the Library, and distribute
//L     that work under terms of your choice, provided that you give
//L     prominent notice with each copy of the work that the specified
//L     version of the Library is used in it, and that you include or
//L     provide public access to the complete corresponding
//L     machine-readable source code for the Library including whatever
//L     changes were used in the work. (i.e. If you make changes to the
//L     Library you must distribute those, but you do not need to
//L     distribute source or object code to those portions of the work
//L     not covered by this licence.)'
//L
//L  Note that this clause grants an additional right and does not impose
//L  any additional restriction, and so does not affect compatibility
//L  with the GNU General Public Licence (GPL). If you wish to negotiate
//L  other terms, please contact the maintainer.
//L
//L  You can redistribute it and/or modify the library under the terms of
//L  the GNU Lesser General Public License as published by the Free Software
//L  Foundation; either version 2.1 of the License, or (at your option) any
//L  later version.
//L
//L  This library is distributed in the hope that it will be useful, but
//L  WITHOUT ANY WARRANTY; without even the implied warranty of
//L  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//L  Lesser General Public License for more details.
//L
//L  You should have received a copy of the CCP4 licence and/or GNU
//L  Lesser General Public License along with this library; if not, write
//L  to the CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//L  The GNU Lesser General Public can also be obtained by writing to the
//L  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//L  MA 02111-1307 USA


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
        virtual bool operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >&, const clipper::HKL_data<clipper::datatypes::I_sigI<T> >& ) = 0;
        virtual bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >&, const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& ) = 0;
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
    class ResolutionFn_nonlinear_rest : public clipper::ResolutionFn
    {
    public:
        //! constructor: need reflections, basis fn, target fn and restraints. No symops.
        ResolutionFn_nonlinear_rest( const clipper::HKL_info& hkl_info, 
                                    const clipper::BasisFn_base& basisfn, 
                                    const clipper::TargetFn_base& targetfn, 
                                    const std::vector<clipper::ftype>& params, 
                                    const ctruncate::RestraintFn_base& restraint,
                                    const clipper::ftype damp = 0.0, 
                                    const bool debug = false );
                    
    protected:
        const RestraintFn_base* restfn_;   //!< retraints function        

    };
    

	//! Structure factor anisotropic scaling
	/*! Perform structure factor anisotropic scaling, observed to calculated,
	 calculated to observed, or observed against itself using log target
	 \ingroup g_funcobj */
	template<class T> class Iscale_loganiso : public ctruncate::Iscale_aniso_base<T> {
	public:
		Iscale_loganiso( clipper::ftype nsig = 3.0, Scaling::MODE mode = Scaling::NORMAL )
        { Iscale_aniso_base<T>::nsig_ = nsig; Iscale_aniso_base<T>::mode_ = mode; }
		//! Scale Fo to Fc
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, const clipper::HKL_data<clipper::datatypes::F_phi<T> >& fc );
        //! Scale Fo to Fc
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fc );
		//! Scale Fc to Fo
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_phi<T> >& fc, const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo );
		//! Scale Fo to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo );
		//! Scale Fo to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, const clipper::ftype resfilter, const int npar_scl );
		//! Scale Io to Ic (psuedo I), or vice versa
		bool operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io, const clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Ic );
		//! Scale Io to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io );
		//! Scale Io to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io, const clipper::ftype resfilter, const int npar_scl );
		//! Primitive scaling functions for F and I
		template<class D, class T1, class T2, class S>
		bool scale( clipper::HKL_data<D>& fo, const clipper::ftype resfilter, const int npar_scl );
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
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, const clipper::HKL_data<clipper::datatypes::F_phi<T> >& fc );
        //! Scale Fo to Fc
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fc );
		//! Scale Fc to Fo
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_phi<T> >& fc, const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo );
		//! Scale Fo to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo );
		//! Scale Fo to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, const clipper::ftype resfilter, const int npar_scl );
		//! Scale Io to Ic (psuedo I), or vice versa
		bool operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io, const clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Ic );
		//! Scale Io to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io );
		//! Scale Io to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io, const clipper::ftype resfilter, const int npar_scl );
		//! Primitive scaling functions for F and I
		template<class D, class T1, class T2, class S>
		bool scale( clipper::HKL_data<D>& fo, const clipper::ftype resfilter, const int npar_scl );
		//! return aniso correction on F or I
	};

    //! Structure factor anisotropic scaling
	/*! Perform structure factor anisotropic scaling, observed to calculated,
	 calculated to observed, or observed against itself using loglikelihood target
	 \ingroup g_funcobj */
	template<class T> class Iscale_llaniso : public ctruncate::Iscale_aniso_base<T> {
	public:
		Iscale_llaniso( clipper::ftype nsig = 0.0, Scaling::MODE mode = Scaling::NORMAL )
        { Iscale_aniso_base<T>::nsig_ = nsig; Iscale_aniso_base<T>::mode_ = mode; }
		//! Scale Fc to Fo
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_phi<T> >& fo, const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fc );
        //! Scale Fo to Fc
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, const clipper::HKL_data<clipper::datatypes::F_phi<T> >& fc );
		//! Scale Fo to Fo
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fc, const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo );
		//! Scale Fo to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo );
		//! Scale Fo to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, const clipper::ftype resfilter, const int npar_scl );
		//! Scale Io to Ic (psuedo I), or vice versa
		bool operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io, const clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Ic );
		//! Scale Io to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io );
		//! Scale Io to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io, const clipper::ftype resfilter, const int npar_scl );
		//! Primitive scaling functions for F and I
		template<class D, class T1, class T2, class R2, class S>
		bool scale( clipper::HKL_data<D>& fo, const clipper::ftype resfilter, const int npar_scl );
        template<class D, class S>
        bool prescale( const clipper::HKL_data<D>& fo, const clipper::HKL_data<D>& fc, clipper::ftype& A, clipper::ftype& B, const clipper::ftype resfilter, const int npar_scl=60 );
	};
    
    //! Structure factor anisotropic scaling
	/*! Perform structure factor anisotropic scaling, observed to calculated,
	 calculated to observed, or observed against itself using standard target.
     Use Iscale_aniso plus wilson scaling to get approximation to anisotropic scaling.
	 \ingroup g_funcobj */
	template<class T> class Iscale_wilson_aniso : public Iscale_aniso_base<T> {
	public:
		Iscale_wilson_aniso( clipper::ftype nsig = 0.0, Scaling::MODE mode = Scaling::NORMAL )
		{ Iscale_aniso_base<T>::nsig_=nsig; Iscale_aniso_base<T>::mode_=mode; }
		//! Scale Fo to Fc
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, const clipper::HKL_data<clipper::datatypes::F_phi<T> >& fc );
		//! Scale Fc to Fo
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_phi<T> >& fc, const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo );
        //! Scale Fo to Fc
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fc );
		//! Scale Fo to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo );
		//! Scale Fo to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, const clipper::ftype resfilter, const int npar_scl );
		//! Scale Io to Ic (psuedo I), or vice versa
		bool operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io, const clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Ic );
		//! Scale Io to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io );
		//! Scale Io to isotropic (approximate)
		bool operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io, const clipper::ftype resfilter, const int npar_scl );
		//! Primitive scaling functions for F and I
		template<class D, class T1, class T2, class S>
		bool scale( clipper::HKL_data<D>& fo, const clipper::ftype resfilter, const int npar_scl );
		//! return aniso correction on F or I
	};

} // namespace ctruncate

#endif
