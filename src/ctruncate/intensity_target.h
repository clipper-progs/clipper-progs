/*! \file lib/resol_targetfn.h
 Header file for resolution function generator
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


#ifndef INTENSITY_TARGET
#define INTENSITY_TARGET

#include "clipper/clipper.h"
#include "ctruncate_utils.h"

//! simple mean |I|<sup>n</sup> target
/*! This class implements the target function for calculating mean
 |I|<sup>n</sup> as a function of position in reciprocal space. It
 includes the appropriate multiplicity correction, and so can be
 applied to any type with an 'I' member with the same dimensions as
 an |I| or |Z| (or an uncorrected |I|).
 
 recoding of TargetFn_meanInth, with additional ice rings term
*/
template<class T> class TargetFn_meanInth : public clipper::TargetFn_base
{
public:
    //! constructor: takes the datalist against which to calc target, and power
	TargetFn_meanInth( const clipper::HKL_data<T>& hkl_data_, const clipper::ftype& n ) : 
	hkl_data(&hkl_data_), power(n) {}
	//! constructor: takes the datalist against which to calc target, and power, plus rings for exclusion
	TargetFn_meanInth( const clipper::HKL_data<T>& hkl_data_, const clipper::ftype& n, const ctruncate::Rings& rings_ ) :
	hkl_data(&hkl_data_), power(n), rings(rings_) {}
	//! return the value and derivatives of the target function
	Rderiv rderiv( const clipper::HKL_info::HKL_reference_index& ih, const clipper::ftype& fh ) const;
    //! the type of the function: optionally used to improve convergence
	FNtype type() const { return QUADRATIC; }
    
protected:
    template <class D > const D&    obs( const clipper::datatypes::F_sigF<D>& f ) const { return f.f(); }
    template <class D > const D&    obs( const clipper::datatypes::F_phi<D>& f ) const { return f.f(); }
    template <class D > const D&    obs( const clipper::datatypes::I_sigI<D>& f ) const { return f.I(); }
    template <class D > const D& sigobs( const clipper::datatypes::F_sigF<D>& f ) { return f.sigf(); }
    template <class D > const D& sigobs( const clipper::datatypes::F_phi<D>& f ) { return 0.0; }
    template <class D > const D& sigobs( const clipper::datatypes::I_sigI<D>& f ) { return f.sigI(); }
    template <class D > int   type( const clipper::datatypes::F_sigF<D>& f ) const { return 1; }
    template <class D > int   type( const clipper::datatypes::F_phi<D>& f ) const { return 1; }
    template <class D > int   type( const clipper::datatypes::I_sigI<D>& f ) const { return 2; }

private:
	clipper::ftype power;
    const clipper::HKL_data<T>* hkl_data;
	const ctruncate::Rings rings;
};

//! log likelihood |I| scaling target
/*! This class implements the target function for calculating the
 scale factor to scale the wilson log likelihood functions
 
 -LL = -1/2.ln(2pi.[E.S+s**2)I) + I/(2(E.S+s**2))   <-- centric
            -ln(E.S+S**2)       + I/(E.S+s**2)      <-- acentric
 
 can be used to scale against Best, scaled Best etc.
 */
template<class T1, class T2> class TargetFn_scaleLogLikeI1I2 : public clipper::TargetFn_base
{
public:
    //! constructor: takes the datalist against which to calc target
    TargetFn_scaleLogLikeI1I2( const clipper::HKL_data<T1>& hkl_data1_, const clipper::HKL_data<T2>& hkl_data2_) :
    hkl_data1(&hkl_data1_), hkl_data2(&hkl_data2_)  {}
    //! return the value and derivatives of the target function
    Rderiv rderiv( const clipper::HKL_info::HKL_reference_index& ih, const clipper::ftype& fh ) const;
    //! the type of the function: optionally used to improve convergence
    FNtype type() const { return GENERAL; }
private:
    const clipper::HKL_data<T1>* hkl_data1;  // experimental
    const clipper::HKL_data<T2>* hkl_data2;  // zero order reference (typically Best)
};

template<class T1, class T2> class TargetFn_scaleLogLikeF1F2 : public clipper::TargetFn_base
{
public:
    //! constructor: takes the datalist against which to calc target
    TargetFn_scaleLogLikeF1F2( const clipper::HKL_data<T1>& hkl_data1_, const clipper::HKL_data<T2>& hkl_data2_) :
    hkl_data1(&hkl_data1_), hkl_data2(&hkl_data2_)  {}
    //! return the value and derivatives of the target function
    Rderiv rderiv( const clipper::HKL_info::HKL_reference_index& ih, const clipper::ftype& fh ) const;
    //! the type of the function: optionally used to improve convergence
    FNtype type() const { return GENERAL; }
private:
    const clipper::HKL_data<T1>* hkl_data1;  // experimental
    const clipper::HKL_data<T2>* hkl_data2;  // zero order reference (typically Best)
};

// mean I calculation

template<class T> clipper::TargetFn_base::Rderiv 
TargetFn_meanInth<T>::rderiv( const clipper::HKL_info::HKL_reference_index& ih, const clipper::ftype& fh ) const
{
    // it's really this bit that does the work
	Rderiv result;
    const clipper::HKL_data<T>& data = *hkl_data;
	if ( !data[ih].missing() && (rings.InRing(ih.hkl().invresolsq(data.base_cell() ) ) == -1 ) ) {
		clipper::ftype d = fh - pow( obs(data[ih]) / pow(ih.hkl_class().epsilon(), clipper::ftype(type(data[ih])/2.0) ), power );
		result.r = d * d;
		result.dr = 2.0 * d;
		result.dr2 = 2.0;
    } else {
		result.r = result.dr = result.dr2 = 0.0;
    }
    return result;
}

template<class T> class TargetFn_meanSignth : public clipper::TargetFn_base
{
public:
    //! constructor: takes the datalist against which to calc target, and power
	TargetFn_meanSignth( const clipper::HKL_data<T>& hkl_data_, const clipper::ftype& n ) : 
	hkl_data(&hkl_data_), power(n) {}
	//! constructor: takes the datalist against which to calc target, and power, plus rings for exclusion
	TargetFn_meanSignth( const clipper::HKL_data<T>& hkl_data_, const clipper::ftype& n, const ctruncate::Rings& rings_ ) :
	hkl_data(&hkl_data_), power(n), rings(rings_) {}
	//! return the value and derivatives of the target function
	Rderiv rderiv( const clipper::HKL_info::HKL_reference_index& ih, const clipper::ftype& fh ) const;
    //! the type of the function: optionally used to improve convergence
	FNtype type() const { return QUADRATIC; }
private:
	clipper::ftype power;
    const clipper::HKL_data<T>* hkl_data;
	const ctruncate::Rings rings;
};


template<class T> clipper::TargetFn_base::Rderiv 
TargetFn_meanSignth<T>::rderiv( const clipper::HKL_info::HKL_reference_index& ih, const clipper::ftype& fh ) const
{
    // it's really this bit that does the work
	Rderiv result;
    const clipper::HKL_data<T>& data = *hkl_data;
    if ( !data[ih].missing() && (rings.InRing(ih.hkl().invresolsq(data.base_cell() ) ) == -1 ) ) {
		clipper::ftype d = fh - pow( clipper::ftype(data[ih].sigI()) / pow(ih.hkl_class().epsilon(), clipper::ftype(type(data[ih])/2.0) ), power );
		result.r = d * d;
		result.dr = 2.0 * d;
		result.dr2 = 2.0;
    } else {
		result.r = result.dr = result.dr2 = 0.0;
    }
    return result;
}

// Log Likelihood I1-I2 scaling

template<class T1, class T2> clipper::TargetFn_base::Rderiv TargetFn_scaleLogLikeI1I2<T1,T2>::rderiv( const clipper::HKL_info::HKL_reference_index& ih, const clipper::ftype& fh ) const
{
    Rderiv result;
    result.r = result.dr = result.dr2 = 0.0;
    const T1& ft1 = (*hkl_data1)[ih];  //reference which is scaled to observed
    const T2& ft2 = (*hkl_data2)[ih];  //observed
    if ( !ft1.missing() && !ft2.missing() )
        if ( ft1.I() > 1.0e-6 && ft2.I() > 1.0e-6 ) {
            const clipper::ftype eps = ih.hkl_class().epsilonc();
            const clipper::ftype centfac = ih.hkl_class().centric() ? 0.5 : 1.0;
            const clipper::ftype f1 = ft1.I();
            const clipper::ftype f2 = ft2.I();
            const clipper::ftype s = std::fabs(ft2.sigI());
            //const clipper::ftype s = 0.0;
            const clipper::ftype ref = eps*fh*f1 + s;
            const clipper::ftype rel_ref = ref/(eps*f1);
            result.r   = centfac*(std::log(ref)+f2/ref);  // likelihood value
            result.dr  = centfac/rel_ref*(1.0-f2/ref); // first derivative wrt Sigma
            result.dr2 = centfac/(rel_ref*rel_ref)*(2.0*f2/ref-1.0); // second derivative wrt Sigma
        }
    return result;
}

template<class T1, class T2> clipper::TargetFn_base::Rderiv TargetFn_scaleLogLikeF1F2<T1,T2>::rderiv( const clipper::HKL_info::HKL_reference_index& ih, const clipper::ftype& fh ) const
{
    //fh is scaled value of reference in this case
    Rderiv result;
    result.r = result.dr = result.dr2 = 0.0;
    const T1& ft1 = (*hkl_data1)[ih];  //reference which is scaled to observed
    const T2& ft2 = (*hkl_data2)[ih];  //observed
    if ( !ft1.missing() && !ft2.missing() )
        if ( ft1.f() > 1.0e-6 && ft2.f() > 1.0e-6 ) {
            const clipper::ftype eps = ih.hkl_class().epsilonc();
            const clipper::ftype centfac = ih.hkl_class().centric() ? 0.5 : 1.0;
            const clipper::ftype f1 = ft1.f();
            const clipper::ftype f2 = ft2.f();
            const clipper::ftype s = std::pow(clipper::ftype(ft2.sigf()),clipper::ftype(2.0));
            const clipper::ftype ref = eps*fh*f1*f1 + s;
            const clipper::ftype rel_ref = ref/(eps*f1*f1);
            result.r   = centfac*(std::log(ref)+f2*f2/ref);  // likelihood value
            result.dr  = centfac/rel_ref*(1.0-f2*f2/ref); // first derivative wrt Sigma
            result.dr2 = centfac/(rel_ref*rel_ref)*(2.0*f2*f2/ref-1.0); // second derivative wrt Sigma
        }
    return result;
}

#endif
