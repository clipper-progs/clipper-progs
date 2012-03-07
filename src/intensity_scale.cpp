/* Iscale.cpp: structure factor anisotropic scaling implementation */
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


#include "intensity_scale.h"
#include "intensity_target.h"
#include "clipper/core/resol_targetfn.h"
#include "clipper/core/hkl_compute.h"
#include <limits>


namespace ctruncate {
	
    template<class T> const clipper::U_aniso_orth& Iscale_aniso_base<T>::u_aniso_orth( Scaling::TYPE t ) const
	{
		if ( t == Scaling::I ) return this->u_i;
		else          return this->u_f;
	}
    
    /*------------------log aniso scaling ---------------------------------*/
    
	template<class T> bool Iscale_loganiso<T>::operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, const clipper::HKL_data<clipper::datatypes::F_phi<T> >& fc )
	{
		typedef clipper::HKL_info::HKL_reference_index HRI;
		// expand to P1 in order to preserve symmetry
		const clipper::HKL_info& hkls = fo.hkl_info();
		clipper::Spacegroup spgrp1( clipper::Spacegroup::P1 );
		clipper::HKL_info hkl1( spgrp1, hkls.cell(), hkls.resolution(), true );
		clipper::HKL_data<clipper::datatypes::F_sigF<T> > fo1( hkl1 );
		clipper::HKL_data<clipper::datatypes::F_phi<T> >  fc1( hkl1 );
		for ( HRI ih = hkl1.first(); !ih.last(); ih.next() ) {
			clipper::datatypes::F_sigF<T> f = fo[ih.hkl()];
			if ( f.f() >= this->nsig_ * f.sigf() ) {
				fo1[ih] = f;
				fc1[ih] = fc[ih.hkl()];
			}
		}
		// do the aniso scaling
		std::vector<double> param( 7, 0.0 );
		clipper::BasisFn_log_aniso_gaussian bfn;
		clipper::TargetFn_scaleLogF1F2<clipper::datatypes::F_sigF<T>,clipper::datatypes::F_phi<T> >
		tfn( fo1, fc1 );
		clipper::ResolutionFn rfn( hkl1, bfn, tfn, param );
		for ( HRI ih = hkls.first(); !ih.last(); ih.next() )
			if ( !fo[ih].missing() )
				fo[ih].scale( exp( 0.5*bfn.f(ih.hkl(),hkls.cell(),rfn.params()) ) );
		this->u_i = bfn.u_aniso_orth( rfn.params() );
		this->u_f = 0.5 * this->u_i;
		this->bscale_ = bfn.scale(rfn.params());
		return true;
	}
    
	template<class T> bool Iscale_loganiso<T>::operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io, const clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Ic )
	{
		typedef clipper::HKL_info::HKL_reference_index HRI;
		// expand to P1 in order to preserve symmetry
		const clipper::HKL_info& hkls = Io.hkl_info();
		clipper::Spacegroup spgrp1( clipper::Spacegroup::P1 );
		clipper::HKL_info hkl1( spgrp1, hkls.cell(), hkls.resolution(), true );
		clipper::HKL_data<clipper::datatypes::I_sigI<T> > Io1( hkl1 );
		clipper::HKL_data<clipper::datatypes::I_sigI<T> > Ic1( hkl1 );
		for ( HRI ih = hkl1.first(); !ih.last(); ih.next() ) {
			clipper::datatypes::I_sigI<T> I = Io[ih.hkl()];
			if ( I.I() >= this->nsig_ * I.sigI() ) {
				Io1[ih] = I;
				Ic1[ih] = Ic[ih.hkl()];
			}
		}
		// do the aniso scaling
		std::vector<double> param( 7, 0.0 );
		clipper::BasisFn_log_aniso_gaussian bfn;
        clipper::TargetFn_scaleLogI1I2<clipper::datatypes::I_sigI<T>,clipper::datatypes::I_sigI<T> >
        			tfn( Io1, Ic1 );
		clipper::ResolutionFn rfn( hkl1, bfn, tfn, param );

		this->u_i = bfn.u_aniso_orth( rfn.params() );
		this->u_f = 0.5 * this->u_i;
		this->bscale_ = bfn.scale( rfn.params() );
		
		clipper::datatypes::Compute_scale_u_aniso<clipper::datatypes::I_sigI<T> > compute_s(std::sqrt(this->bscale_),-(this->u_f )); //always take u_f
		Io.compute( Io, compute_s );
		return true;
	}
		
	template<class T> bool Iscale_loganiso<T>::operator() ( clipper::HKL_data<clipper::datatypes::F_phi<T> >& fc, const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo )
	{
		typedef clipper::HKL_info::HKL_reference_index HRI;
		// expand to P1 in order to preserve symmetry
		const clipper::HKL_info& hkls = fo.hkl_info();
		clipper::Spacegroup spgrp1( clipper::Spacegroup::P1 );
		clipper::HKL_info hkl1( spgrp1, hkls.cell(), hkls.resolution(), true );
		clipper::HKL_data<clipper::datatypes::F_sigF<T> > fo1( hkl1 );
		clipper::HKL_data<clipper::datatypes::F_phi<T> >  fc1( hkl1 );
		for ( HRI ih = hkl1.first(); !ih.last(); ih.next() ) {
			clipper::datatypes::F_sigF<T> f = fo[ih.hkl()];
			if ( f.f() >= this->nsig_ * f.sigf() ) {
				fo1[ih] = f;
				fc1[ih] = fc[ih.hkl()];
			}
		}
		// do the aniso scaling
		std::vector<double> param( 7, 0.0 );
		clipper::BasisFn_log_aniso_gaussian bfn;
		clipper::TargetFn_scaleLogF1F2<clipper::datatypes::F_phi<T>,clipper::datatypes::F_sigF<T> >
		tfn( fc1, fo1 );
		clipper::ResolutionFn rfn( hkl1, bfn, tfn, param );
		for ( HRI ih = hkls.first(); !ih.last(); ih.next() )
			if ( !fc[ih].missing() )
				fc[ih].scale( exp( 0.5*bfn.f(ih.hkl(),hkls.cell(),rfn.params()) ) );
		this->u_i = bfn.u_aniso_orth( rfn.params() );
		this->u_f = 0.5 * this->u_i;
		this->bscale_ = bfn.scale( rfn.params() );
		return true;
	}
    
    template<class T> bool Iscale_loganiso<T>::operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fc, const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo )
	{
		typedef clipper::HKL_info::HKL_reference_index HRI;
		// expand to P1 in order to preserve symmetry
		const clipper::HKL_info& hkls = fo.hkl_info();
		clipper::Spacegroup spgrp1( clipper::Spacegroup::P1 );
		clipper::HKL_info hkl1( spgrp1, hkls.cell(), hkls.resolution(), true );
		clipper::HKL_data<clipper::datatypes::F_sigF<T> > fo1( hkl1 );
		clipper::HKL_data<clipper::datatypes::F_sigF<T> > fc1( hkl1 );
		for ( HRI ih = hkl1.first(); !ih.last(); ih.next() ) {
			clipper::datatypes::F_sigF<T> f = fo[ih.hkl()];
			if ( f.f() >= this->nsig_ * f.sigf() ) {
				fo1[ih] = f;
				fc1[ih] = fc[ih.hkl()];
			}
		}
		// do the aniso scaling
		std::vector<double> param( 7, 0.0 );
		clipper::BasisFn_log_aniso_gaussian bfn;
		clipper::TargetFn_scaleLogF1F2<clipper::datatypes::F_sigF<T>,clipper::datatypes::F_sigF<T> >
		tfn( fc1, fo1 );
		clipper::ResolutionFn rfn( hkl1, bfn, tfn, param );
		for ( HRI ih = hkls.first(); !ih.last(); ih.next() )
			if ( !fc[ih].missing() )
				fc[ih].scale( exp( 0.5*bfn.f(ih.hkl(),hkls.cell(),rfn.params()) ) );
		this->u_i = bfn.u_aniso_orth( rfn.params() );
		this->u_f = 0.5 * this->u_i;
		this->bscale_ = bfn.scale( rfn.params() );
		return true;
	}

	
	template<class T> bool Iscale_loganiso<T>::operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo )
	{
		typedef clipper::datatypes::F_sigF<T> DATA;
		typedef clipper::TargetFn_scaleF1F2<DATA,DATA> TGT1;
		typedef clipper::TargetFn_scaleLogF1F2<DATA,DATA> TGT2;
		typedef clipper::BasisFn_spline SCALETYPE;
		return scale<DATA,TGT1,TGT2,SCALETYPE>( fo, -1.0, 12 );
	}
	
	template<class T> bool Iscale_loganiso<T>::operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io )
	{
		typedef clipper::datatypes::I_sigI<T> DATA;
		typedef clipper::TargetFn_scaleI1I2<DATA,DATA> TGT1;
		typedef clipper::TargetFn_scaleLogI1I2<DATA,DATA> TGT2;
		typedef clipper::BasisFn_spline SCALETYPE;		
		return scale<DATA,TGT1,TGT2,SCALETYPE>( Io, -1.0, 12 );
	}
	
	template<class T> bool Iscale_loganiso<T>::operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, clipper::ftype resfilter, const int npar_scl )
	{
		typedef clipper::datatypes::F_sigF<T> DATA;
		typedef clipper::TargetFn_scaleF1F2<DATA,DATA> TGT1;
		typedef clipper::TargetFn_scaleLogF1F2<DATA,DATA> TGT2;
		typedef clipper::BasisFn_spline SCALETYPE;
		return scale<DATA,TGT1,TGT2,SCALETYPE>( fo, resfilter, npar_scl );
	}
	
	
	template<class T> bool Iscale_loganiso<T>::operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io,clipper:: ftype resfilter, const int npar_scl )
	{
		typedef clipper::datatypes::I_sigI<T> DATA;
		typedef clipper::TargetFn_scaleI1I2<DATA,DATA> TGT1;
		typedef clipper::TargetFn_scaleLogI1I2<DATA,DATA> TGT2;
		typedef clipper::BasisFn_spline SCALETYPE;
		return scale<DATA,TGT1,TGT2,SCALETYPE>( Io, resfilter, npar_scl );
	}
	
	
	template<class T> template<class D, class T1, class T2, class S> bool Iscale_loganiso<T>::scale( clipper::HKL_data<D>& fo, const clipper::ftype resfilter, const int npar_scl )
	{
		typedef clipper::HKL_info::HKL_reference_index HRI;
		// expand to P1 in order to preserve symmetry
		const clipper::HKL_info& hkls = fo.hkl_info();
		clipper::Spacegroup spgrp1( clipper::Spacegroup::P1 );
		clipper::HKL_info hkl1( spgrp1, hkls.cell(), hkls.resolution(), true );
		clipper::HKL_data<D> fo1( hkl1 ), fs1( hkl1 ), fc1( hkl1 );
		for ( HRI ih = hkl1.first(); !ih.last(); ih.next() ) {
			D f = fo[ih.hkl()];
			if ( obs(f) >= this->nsig_ * sigobs(f) ) fo1[ih] = f;
		}
		
		// calc target values
		fc1 = D( 1.0, 1.0 );
		for ( HRI ih = fc1.first(); !ih.last(); ih.next() )
			fc1[ih].scale( sqrt( ih.hkl_class().epsilon() ) );
		
		// do correction
		this->u_i = this->u_f = clipper::U_aniso_orth( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
		std::vector<clipper::ftype> param( 7, 0.0 );
		clipper::BasisFn_log_aniso_gaussian bfn;
		clipper::ftype dp;
		for ( int c = 0; c < 2; c++ ) {
			// remove current anistropy estimate
			clipper::datatypes::Compute_scale_u_aniso<D> compute_s(1.0,-(this->u_f) );
			fs1.compute( fo1, compute_s );
			
			// and calc E scale
			S basissc( fs1, npar_scl, 1.0 );
			std::vector<clipper::ftype> param_fo( basissc.num_params(), 1.0 );
			T1 tfn_fo( fs1, fc1 );
			clipper::ResolutionFn rfn_fo( hkl1, basissc, tfn_fo, param_fo );
			param_fo = rfn_fo.params();
			
			// prescale F to E-like scale
			for ( HRI ih = fs1.first(); !ih.last(); ih.next() ) {
				fs1[ih] = fo1[ih];
				fs1[ih].scale( sqrt( basissc.f_s( ih.invresolsq(), param_fo ) ) );
			}
						
			// do the aniso scaling
			T2 tfn( fs1, fc1 );
			param = std::vector<clipper::ftype>( 7, 0.0 );
			clipper::ResolutionFn rfn( hkl1, bfn, tfn, param );
			param = rfn.params();
			// set trace to zero (i.e. no isotropic correction)
			dp = (param[1]+param[2]+param[3])/3.0;
			param[1] -= dp; param[2] -= dp; param[3] -= dp;
			this->u_i = bfn.u_aniso_orth( param );
			this->u_f = 0.5 * this->u_i;
			//std::cout << c << " | " << param[1] << " " << param[2] << " " << param[3] << " " << param[4] << " " << param[5] << " " << param[6] << "\n"; std::cout << " DP " << dp << "\n";
		}
		
		// sharpen or smooth as required
		clipper::Matrix<clipper::ftype> m(3,3);
		m(0,0)=       param[1]; m(1,1)=       param[2]; m(2,2)=       param[3];
		m(0,1)=m(1,0)=param[4]; m(0,2)=m(2,0)=param[5]; m(1,2)=m(2,1)=param[6];
		std::vector<clipper::ftype> ev = m.eigen();
		//std::cout << "EIGEN " << param[1] << " " << param[2] << " " << param[3] << " " << ev[0] << " " << ev[1] << " " << ev[2] << std::endl;
		dp = 0.0;
		if ( this->mode_ == Scaling::SHARPEN   ) dp = ev[2];
		if ( this->mode_ == Scaling::UNSHARPEN ) dp = ev[0];
		param[1] -= dp; param[2] -= dp; param[3] -= dp;
		this->u_i = bfn.u_aniso_orth( param );
		this->u_f = 0.5 * this->u_i;
		this->bscale_ = bfn.scale(param);
		// store the results
		clipper::datatypes::Compute_scale_u_aniso<D> compute_s(1.0,-(this->u_f));
		fo.compute( fo, compute_s );
		
		return true;
	}
	
    /*------------------aniso scaling ---------------------------------*/
    
    template<class T> bool Iscale_aniso<T>::operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, const clipper::HKL_data<clipper::datatypes::F_phi<T> >& fc )
	{
		typedef clipper::HKL_info::HKL_reference_index HRI;
		// expand to P1 in order to preserve symmetry
		const clipper::HKL_info& hkls = fo.hkl_info();
		clipper::Spacegroup spgrp1( clipper::Spacegroup::P1 );
		clipper::HKL_info hkl1( spgrp1, hkls.cell(), hkls.resolution(), true );
		clipper::HKL_data<clipper::datatypes::F_sigF<T> > fo1( hkl1 );
		clipper::HKL_data<clipper::datatypes::F_phi<T> >  fc1( hkl1 );
		for ( HRI ih = hkl1.first(); !ih.last(); ih.next() ) {
			clipper::datatypes::F_sigF<T> f = fo[ih.hkl()];
			if ( f.f() >= this->nsig_ * f.sigf() ) {
				fo1[ih] = f;
				fc1[ih] = fc[ih.hkl()];
			}
		}
		// do the aniso scaling
		std::vector<double> param( 7, 0.0 );
		clipper::BasisFn_aniso_gaussian bfn;
		clipper::TargetFn_scaleF1F2<clipper::datatypes::F_sigF<T>,clipper::datatypes::F_phi<T> >
		tfn( fo1, fc1 );
		clipper::ResolutionFn_nonlinear rfn( hkl1, bfn, tfn, param );
		for ( HRI ih = hkls.first(); !ih.last(); ih.next() )
			if ( !fo[ih].missing() )
				fo[ih].scale( exp( 0.5*bfn.f(ih.hkl(),hkls.cell(),rfn.params()) ) );
		this->u_i = bfn.u_aniso_orth( rfn.params() );
		this->u_f = 0.5 * this->u_i;
		this->bscale_ = bfn.scale(rfn.params());
		return true;
	}
	
    template<class T> bool Iscale_aniso<T>::operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fc )
	{
		typedef clipper::HKL_info::HKL_reference_index HRI;
		// expand to P1 in order to preserve symmetry
		const clipper::HKL_info& hkls = fo.hkl_info();
		clipper::Spacegroup spgrp1( clipper::Spacegroup::P1 );
		clipper::HKL_info hkl1( spgrp1, hkls.cell(), hkls.resolution(), true );
		clipper::HKL_data<clipper::datatypes::F_sigF<T> > fo1( hkl1 );
		clipper::HKL_data<clipper::datatypes::F_sigF<T> >  fc1( hkl1 );
		for ( HRI ih = hkl1.first(); !ih.last(); ih.next() ) {
			clipper::datatypes::F_sigF<T> f = fo[ih.hkl()];
			if ( f.f() >= this->nsig_ * f.sigf() ) {
				fo1[ih] = f;
				fc1[ih] = fc[ih.hkl()];
			}
		}
		// do the aniso scaling
		std::vector<double> param( 7, 0.0 );
		clipper::BasisFn_aniso_gaussian bfn;
		clipper::TargetFn_scaleF1F2<clipper::datatypes::F_sigF<T>,clipper::datatypes::F_sigF<T> >
		tfn( fo1, fc1 );
		clipper::ResolutionFn_nonlinear rfn( hkl1, bfn, tfn, param );
		for ( HRI ih = hkls.first(); !ih.last(); ih.next() )
			if ( !fo[ih].missing() )
				fo[ih].scale( exp( 0.5*bfn.f(ih.hkl(),hkls.cell(),rfn.params()) ) );
		this->u_i = bfn.u_aniso_orth( rfn.params() );
		this->u_f = 0.5 * this->u_i;
		this->bscale_ = bfn.scale(rfn.params());
		return true;
	}

	template<class T> bool Iscale_aniso<T>::operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io, const clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Ic )
	{
		typedef clipper::HKL_info::HKL_reference_index HRI;
		// expand to P1 in order to preserve symmetry
		const clipper::HKL_info& hkls = Io.hkl_info();
		clipper::Spacegroup spgrp1( clipper::Spacegroup::P1 );
		clipper::HKL_info hkl1( spgrp1, hkls.cell(), hkls.resolution(), true );
		clipper::HKL_data<clipper::datatypes::I_sigI<T> > Io1( hkl1 );
		clipper::HKL_data<clipper::datatypes::I_sigI<T> > Ic1( hkl1 );
		for ( HRI ih = hkl1.first(); !ih.last(); ih.next() ) {
			clipper::datatypes::I_sigI<T> I = Io[ih.hkl()];
			if ( I.I() >= this->nsig_ * I.sigI() ) {
				Io1[ih] = I;
				Ic1[ih] = Ic[ih.hkl()];
			}
		}
		// do the aniso scaling
		std::vector<double> param( 7, 0.0 );
        clipper::BasisFn_aniso_gaussian bfn;
        clipper::TargetFn_scaleI1I2<clipper::datatypes::I_sigI<T>,clipper::datatypes::I_sigI<T> >
        tfn( Io1, Ic1 );
		clipper::ResolutionFn_nonlinear rfn( hkl1, bfn, tfn, param );
        
		this->u_i = bfn.u_aniso_orth( rfn.params() );
		this->u_f = 0.5 * this->u_i;
		this->bscale_ = bfn.scale( rfn.params() );
		
		clipper::datatypes::Compute_scale_u_aniso<clipper::datatypes::I_sigI<T> > compute_s(std::sqrt(this->bscale_),-(this->u_f) ); //always take u_f
		Io.compute( Io, compute_s );
		return true;
	}
    
	template<class T> bool Iscale_aniso<T>::operator() ( clipper::HKL_data<clipper::datatypes::F_phi<T> >& fc, const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo )
	{
		typedef clipper::HKL_info::HKL_reference_index HRI;
		// expand to P1 in order to preserve symmetry
		const clipper::HKL_info& hkls = fo.hkl_info();
		clipper::Spacegroup spgrp1( clipper::Spacegroup::P1 );
		clipper::HKL_info hkl1( spgrp1, hkls.cell(), hkls.resolution(), true );
		clipper::HKL_data<clipper::datatypes::F_sigF<T> > fo1( hkl1 );
		clipper::HKL_data<clipper::datatypes::F_phi<T> >  fc1( hkl1 );
		for ( HRI ih = hkl1.first(); !ih.last(); ih.next() ) {
			clipper::datatypes::F_sigF<T> f = fo[ih.hkl()];
			if ( f.f() >= this->nsig_ * f.sigf() ) {
				fo1[ih] = f;
				fc1[ih] = fc[ih.hkl()];
			}
		}
		// do the aniso scaling
		std::vector<double> param( 7, 0.0 );
		clipper::BasisFn_aniso_gaussian bfn;
		clipper::TargetFn_scaleF1F2<clipper::datatypes::F_phi<T>,clipper::datatypes::F_sigF<T> >
		tfn( fc1, fo1 );
		clipper::ResolutionFn_nonlinear rfn( hkl1, bfn, tfn, param );
		for ( HRI ih = hkls.first(); !ih.last(); ih.next() )
			if ( !fc[ih].missing() )
				fc[ih].scale( exp( 0.5*bfn.f(ih.hkl(),hkls.cell(),rfn.params()) ) );
		this->u_i = bfn.u_aniso_orth( rfn.params() );
		this->u_f = 0.5 * this->u_i;
		this->bscale_ = bfn.scale( rfn.params() );
		return true;
	}
	
	template<class T> bool Iscale_aniso<T>::operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo )
	{
		typedef clipper::datatypes::F_sigF<T> DATA;
		typedef clipper::TargetFn_scaleF1F2<DATA,DATA> TGT1;
		typedef clipper::TargetFn_scaleF1F2<DATA,DATA> TGT2;
		typedef clipper::BasisFn_spline SCALETYPE;
		return scale<DATA,TGT1,TGT2,SCALETYPE>( fo, -1.0, 12 );
	}
	
	template<class T> bool Iscale_aniso<T>::operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io )
	{
		typedef clipper::datatypes::I_sigI<T> DATA;
		typedef clipper::TargetFn_scaleI1I2<DATA,DATA> TGT1;
		typedef clipper::TargetFn_scaleI1I2<DATA,DATA> TGT2;
		typedef clipper::BasisFn_spline SCALETYPE;		
		return scale<DATA,TGT1,TGT2,SCALETYPE>( Io, -1.0, 12 );
	}
	
	template<class T> bool Iscale_aniso<T>::operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, clipper::ftype resfilter, const int npar_scl )
	{
		typedef clipper::datatypes::F_sigF<T> DATA;
		typedef clipper::TargetFn_scaleF1F2<DATA,DATA> TGT1;
		typedef clipper::TargetFn_scaleF1F2<DATA,DATA> TGT2;
		typedef clipper::BasisFn_spline SCALETYPE;
		return scale<DATA,TGT1,TGT2,SCALETYPE>( fo, resfilter, npar_scl );
	}
	
	
	template<class T> bool Iscale_aniso<T>::operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io,clipper:: ftype resfilter, const int npar_scl )
	{
		typedef clipper::datatypes::I_sigI<T> DATA;
		typedef clipper::TargetFn_scaleI1I2<DATA,DATA> TGT1;
		typedef clipper::TargetFn_scaleI1I2<DATA,DATA> TGT2;
		typedef clipper::BasisFn_spline SCALETYPE;
		return scale<DATA,TGT1,TGT2,SCALETYPE>( Io, resfilter, npar_scl );
	}
	
	
	template<class T> template<class D, class T1, class T2, class S> bool Iscale_aniso<T>::scale( clipper::HKL_data<D>& fo, const clipper::ftype resfilter, const int npar_scl )
	{
		typedef clipper::HKL_info::HKL_reference_index HRI;
		// expand to P1 in order to preserve symmetry
		const clipper::HKL_info& hkls = fo.hkl_info();
		clipper::Spacegroup spgrp1( clipper::Spacegroup::P1 );
		clipper::HKL_info hkl1( spgrp1, hkls.cell(), hkls.resolution(), true );
		clipper::HKL_data<D> fo1( hkl1 ), fs1( hkl1 ), fc1( hkl1 );
		for ( HRI ih = hkl1.first(); !ih.last(); ih.next() ) {
			D f = fo[ih.hkl()];
			if ( obs(f) >= this->nsig_ * sigobs(f) ) fo1[ih] = f;
		}
		
		// calc target values
		fc1 = D( 1.0, 1.0 );
		for ( HRI ih = fc1.first(); !ih.last(); ih.next() )
			fc1[ih].scale( sqrt( ih.hkl_class().epsilon() ) );
		
		// do correction
		this->u_i = this->u_f = clipper::U_aniso_orth( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
		std::vector<clipper::ftype> param( 7, 0.0 );
		clipper::BasisFn_aniso_gaussian bfn;
		clipper::ftype dp;
		for ( int c = 0; c < 2; c++ ) {
			// remove current anistropy estimate
			clipper::datatypes::Compute_scale_u_aniso<D> compute_s(1.0,-(this->u_f) );
			fs1.compute( fo1, compute_s );
						
			// and calc E scale
			S basissc( fs1, npar_scl, 1.0 );
			std::vector<clipper::ftype> param_fo( basissc.num_params(), 1.0 );
			T1 tfn_fo( fs1, fc1 );
			clipper::ResolutionFn rfn_fo( hkl1, basissc, tfn_fo, param_fo );
			param_fo = rfn_fo.params();
			
			// prescale F to E-like scale
			for ( HRI ih = fs1.first(); !ih.last(); ih.next() ) {
				fs1[ih] = fo1[ih];
				fs1[ih].scale( sqrt( basissc.f_s( ih.invresolsq(), param_fo ) ) );
			}
						
			// do the aniso scaling
			T2 tfn( fs1, fc1 );
			param = std::vector<clipper::ftype>( 7, 0.0 );
			clipper::ResolutionFn_nonlinear rfn( hkl1, bfn, tfn, param );
			param = rfn.params();
			// set trace to zero (i.e. no isotropic correction)
			dp = (param[1]+param[2]+param[3])/3.0;
			param[1] -= dp; param[2] -= dp; param[3] -= dp;
			this->u_i = bfn.u_aniso_orth( param );
			this->u_f = 0.5 * this->u_i;
		}
		
		// sharpen or smooth as required
		clipper::Matrix<clipper::ftype> m(3,3);
		m(0,0)=       param[1]; m(1,1)=       param[2]; m(2,2)=       param[3];
		m(0,1)=m(1,0)=param[4]; m(0,2)=m(2,0)=param[5]; m(1,2)=m(2,1)=param[6];
		std::vector<clipper::ftype> ev = m.eigen();
		dp = 0.0;
		if ( this->mode_ == Scaling::SHARPEN   ) dp = ev[2];
		if ( this->mode_ == Scaling::UNSHARPEN ) dp = ev[0];
		param[1] -= dp; param[2] -= dp; param[3] -= dp;
		this->u_i = bfn.u_aniso_orth( param );
		this->u_f = 0.5 * this->u_i;
		this->bscale_ = bfn.scale(param);
		// store the results
		clipper::datatypes::Compute_scale_u_aniso<D> compute_s(1.0,-(this->u_f) );
		fo.compute( fo, compute_s );
		
		return true;
    }
    
    /*------------------log likelihood ---------------------------------*/
    
    
    template<class T> template<class D, class S> bool Iscale_llaniso<T>::prescale(const clipper::HKL_data<D>& f_ref, const clipper::HKL_data<D>& f_target, clipper::ftype& a, clipper::ftype& b, clipper::ftype resfilter, const int nbins)
    {
        const clipper::HKL_info& hkls = f_ref.hkl_info();
        
        std::vector<clipper::ftype> params_init1( nbins, 1.0 );
        S basis_fo_beta1( f_ref, nbins, 1.0 );
        TargetFn_meanInth<D > target_fo_beta1( f_ref, 1.0);
        clipper::ResolutionFn_nonlinear betaplot1( hkls, basis_fo_beta1, target_fo_beta1, params_init1 ); 
        
        std::vector<clipper::ftype> params_init2( nbins, 1.0 );
        S basis_fo_beta2( f_target, nbins, 1.0 );
        TargetFn_meanInth<D > target_fo_beta2( f_target, 1.0);
        clipper::ResolutionFn_nonlinear betaplot2( hkls, basis_fo_beta2, target_fo_beta2, params_init2 );
        
        int nobs = f_target.num_obs();
        std::vector<clipper::ftype> xi, yi, wi; 
        xi.reserve(nobs) ; yi.reserve(nobs) ; wi.reserve(nobs);
        
        for ( clipper::HKL_info::HKL_reference_index ih = f_target.first(); !ih.last(); ih.next() ) {
            if ( !f_target[ih].missing() && betaplot1.f(ih) > 0.0 && betaplot2.f(ih) > 0.0 ) {
                clipper::ftype lnS = log(betaplot2.f(ih)) - log(betaplot1.f(ih));
                clipper::ftype res = ih.invresolsq();
                clipper::ftype eps = ih.hkl_class().epsilon();
                xi.push_back(res);
                yi.push_back(lnS);

                clipper::ftype weight = 1.0/eps;
                if (weight > 0.0) {
                    wi.push_back(1.0/weight);
                } else {
                    wi.push_back(0.0);
                }
            }
        } 
        clipper::ftype siga,sigb;
        
        a = b = 0.0f;
        straight_line_fit(xi,yi,wi,nobs,a,b,siga,sigb);
        return true;
    }

    
    template<class T> bool Iscale_llaniso<T>::operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& f_target, const clipper::HKL_data<clipper::datatypes::F_phi<T> >& f_ref )
	{
		typedef clipper::HKL_info::HKL_reference_index HRI;
        std::vector<clipper::ftype> param( 7, 0.0 );
        const clipper::HKL_info& hkls = f_target.hkl_info();
        // precondition
        {
            int nbins = 60; 
            std::vector<clipper::ftype> params_init1( nbins, 1.0 );
            clipper::BasisFn_binner basis_fo_beta1( f_ref, nbins, 1.0 );
            TargetFn_meanInth<clipper::datatypes::F_phi<T> > target_fo_beta1( f_ref, 1.0);
            clipper::ResolutionFn_nonlinear betaplot1( hkls, basis_fo_beta1, target_fo_beta1, params_init1 ); 
            
            std::vector<clipper::ftype> params_init2( nbins, 1.0 );
            clipper::BasisFn_binner basis_fo_beta2( f_target, nbins, 1.0 );
            TargetFn_meanInth<clipper::datatypes::F_sigF<T> > target_fo_beta2( f_target, 1.0);
            clipper::ResolutionFn_nonlinear betaplot2( hkls, basis_fo_beta2, target_fo_beta2, params_init2 ); 

            int nobs = f_target.num_obs();
            std::vector<clipper::ftype> xi, yi, wi; 
            xi.reserve(nobs) ; yi.reserve(nobs) ; wi.reserve(nobs);
            
            for ( HRI ih = f_target.first(); !ih.last(); ih.next() ) {
                if ( !f_target[ih].missing() && betaplot2.f(ih) > 0.0 && betaplot1.f(ih) > 0.00 ) {
                    clipper::ftype lnS = log(betaplot2.f(ih)) - log(betaplot1.f(ih));
                    clipper::ftype res = ih.invresolsq();
                    clipper::ftype eps = ih.hkl_class().epsilon();
                    xi.push_back(res);
                    yi.push_back(lnS);
                    
                    float weight = 1.0/eps;
                    if (weight > 0.0) {
                        wi.push_back(1.0/weight);
                    } else {
                        wi.push_back(0.0);
                    }
                }
            } 
            clipper::ftype a,b,siga,sigb;
            
            a = b = 0.0f;
            straight_line_fit(xi,yi,wi,nobs,a,b,siga,sigb);
            
            param[0] = b;
            param[1] = param[2] = param[3] = -a;
        }

        // expand to P1 in order to preserve symmetry
		clipper::Spacegroup spgrp1( clipper::Spacegroup::P1 );
		clipper::HKL_info hkl1( spgrp1, hkls.cell(), hkls.resolution(), true );
		clipper::HKL_data<clipper::datatypes::F_sigF<T> > fo1( hkl1 );
		clipper::HKL_data<clipper::datatypes::F_phi<T> >  fc1( hkl1 );
        clipper::ftype n=0.0;
		for ( HRI ih = hkl1.first(); !ih.last(); ih.next() ) {
			clipper::datatypes::F_sigF<T> f = f_target[ih.hkl()];
			if ( f.f() >= this->nsig_ * f.sigf() ) {
				fo1[ih] = f;
				fc1[ih] = f_ref[ih.hkl()];
                n += 1.0;
			}
		}
		// do the aniso scaling
		clipper::BasisFn_aniso_gaussian bfn;
        ctruncate::RestraintFn_sphericalU restU(10.0,n);
        ::TargetFn_scaleLogLikeF1F2<clipper::datatypes::F_phi<T>,clipper::datatypes::F_sigF<T> >
		tfn( fc1, fo1 );
		ctruncate::ResolutionFn_nonlinear_rest rfn( hkl1, bfn, tfn, param, restU );
        
        //invert the scaling
        param = rfn.params();
        for (std::vector<clipper::ftype>::iterator i=param.begin() ; i != param.end() ; ++i ) (*i) *= -1.0;
		this->u_i = bfn.u_aniso_orth( param );
		this->u_f = 0.5 * this->u_i;
		this->bscale_ = bfn.scale(param);
        clipper::datatypes::Compute_scale_u_aniso<clipper::datatypes::F_sigF<T> > compute_s(std::sqrt(std::fabs(this->bscale_)), -this->u_f ); //always take u_f
        f_target.compute( f_target, compute_s );

		return true;
	}
    
    
    template<class T> bool Iscale_llaniso<T>::operator() ( clipper::HKL_data<clipper::datatypes::F_phi<T> >& f_ref, const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& f_target)
	{
		typedef clipper::HKL_info::HKL_reference_index HRI;
        std::vector<clipper::ftype> param( 7, 0.0 );
        const clipper::HKL_info& hkls = f_target.hkl_info();
        
        // precondition
        {
            int nbins = 60; 
            std::vector<clipper::ftype> params_init1( nbins, 1.0 );
            clipper::BasisFn_binner basis_fo_beta1( f_ref, nbins, 1.0 );
            TargetFn_meanInth<clipper::datatypes::F_phi<T> > target_fo_beta1( f_ref, 1.0);
            clipper::ResolutionFn_nonlinear betaplot1( hkls, basis_fo_beta1, target_fo_beta1, params_init1 ); 
            
            std::vector<clipper::ftype> params_init2( nbins, 1.0 );
            clipper::BasisFn_binner basis_fo_beta2( f_target, nbins, 1.0 );
            TargetFn_meanInth<clipper::datatypes::F_sigF<T> > target_fo_beta2( f_target, 1.0);
            clipper::ResolutionFn_nonlinear betaplot2( hkls, basis_fo_beta2, target_fo_beta2, params_init2 ); 
            
            int nobs = f_target.num_obs();
            std::vector<clipper::ftype> xi, yi, wi; 
            xi.reserve(nobs) ; yi.reserve(nobs) ; wi.reserve(nobs);
            
            for ( HRI ih = f_target.first(); !ih.last(); ih.next() ) {
                if ( !f_target[ih].missing() && betaplot2.f(ih) > 0.0 && betaplot1.f(ih) > 0.00 ) {
                    clipper::ftype lnS = log(betaplot2.f(ih)) - log(betaplot1.f(ih));
                    clipper::ftype res = ih.invresolsq();
                    clipper::ftype eps = ih.hkl_class().epsilon();
                    xi.push_back(res);
                    yi.push_back(lnS);
                    
                    clipper::ftype weight = 1.0/eps;
                    if (weight > 0.0) {
                        wi.push_back(1.0/weight);
                    } else {
                        wi.push_back(0.0);
                    }
                }
            } 
            clipper::ftype a,b,siga,sigb;
            
            a = b = 0.0f;
            straight_line_fit(xi,yi,wi,nobs,a,b,siga,sigb);
            
            param[0] = b;
            param[1] = param[2] = param[3] = -a;
        }

        // expand to P1 in order to preserve symmetry
		clipper::Spacegroup spgrp1( clipper::Spacegroup::P1 );
		clipper::HKL_info hkl1( spgrp1, hkls.cell(), hkls.resolution(), true );
		clipper::HKL_data<clipper::datatypes::F_sigF<T> > fo1( hkl1 );
		clipper::HKL_data<clipper::datatypes::F_phi<T> >  fc1( hkl1 );
        clipper::ftype n = 0.0;
		for ( HRI ih = hkl1.first(); !ih.last(); ih.next() ) {
			clipper::datatypes::F_sigF<T> f = f_target[ih.hkl()];
			if ( f.f() >= this->nsig_ * f.sigf() ) {
				fo1[ih] = f;
				fc1[ih] = f_ref[ih.hkl()];
                n += 1.0;
			}
		}
		// do the aniso scaling
        ctruncate::RestraintFn_sphericalU restU(10.0, n );
		clipper::BasisFn_aniso_gaussian bfn;
        ::TargetFn_scaleLogLikeF1F2<clipper::datatypes::F_phi<T>,clipper::datatypes::F_sigF<T> >
		tfn( fc1, fo1 );
		ctruncate::ResolutionFn_nonlinear_rest rfn( hkl1, bfn, tfn, param, restU );
                
		this->u_i = bfn.u_aniso_orth( rfn.params() );
		this->u_f = 0.5 * this->u_i;
		this->bscale_ = bfn.scale(rfn.params() );
        double s = this->bscale_;
        clipper::datatypes::Compute_scale_u_aniso<clipper::datatypes::F_phi<T> > compute_s(std::sqrt(std::fabs(this->bscale_)), -this->u_f ); //always take u_f
        f_ref.compute( f_ref, compute_s );

		return true;
	}

	
	template<class T> bool Iscale_llaniso<T>::operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& I_ref, const clipper::HKL_data<clipper::datatypes::I_sigI<T> >& I_target )
	{
		typedef clipper::HKL_info::HKL_reference_index HRI;
        std::vector<clipper::ftype> param( 7, 0.0 );
        const clipper::HKL_info& hkls = I_target.hkl_info();
        
        // precondition
        {
            int nbins = 60; 
            std::vector<clipper::ftype> params_init1( nbins, 1.0 );
            clipper::BasisFn_binner basis_fo_beta1( I_ref, nbins, 1.0 );
            TargetFn_meanInth<clipper::datatypes::I_sigI<T> > target_fo_beta1( I_ref, 1.0);
            clipper::ResolutionFn_nonlinear betaplot1( hkls, basis_fo_beta1, target_fo_beta1, params_init1 ); 
            
            std::vector<clipper::ftype> params_init2( nbins, 1.0 );
            clipper::BasisFn_binner basis_fo_beta2( I_target, nbins, 1.0 );
            TargetFn_meanInth<clipper::datatypes::I_sigI<T> > target_fo_beta2( I_target, 1.0);
            clipper::ResolutionFn_nonlinear betaplot2( hkls, basis_fo_beta2, target_fo_beta2, params_init2 ); 
            
            int nobs = I_target.num_obs();
            std::vector<clipper::ftype> xi, yi, wi; 
            xi.reserve(nobs) ; yi.reserve(nobs) ; wi.reserve(nobs);
            
            for ( HRI ih = I_target.first(); !ih.last(); ih.next() ) {
                if ( !I_target[ih].missing() && betaplot2.f(ih) > 0.0 && betaplot1.f(ih) > 0.0 ) {
                    clipper::ftype lnS = log(betaplot2.f(ih)) - log(betaplot1.f(ih));
                    clipper::ftype res = ih.invresolsq();
                    clipper::ftype eps = ih.hkl_class().epsilon();
                    xi.push_back(res);
                    yi.push_back(lnS);
                    
                    clipper::ftype weight = 1.0/eps;
                    if (weight > 0.0) {
                        wi.push_back(1.0/weight);
                    } else {
                        wi.push_back(0.0);
                    }
                }
            } 
            clipper::ftype a,b,siga,sigb;
            
            a = b = 0.0f;
            straight_line_fit(xi,yi,wi,nobs,a,b,siga,sigb);
            
            param[0] = b;
            param[1] = param[2] = param[3] = -a;
        }
        
		// expand to P1 in order to preserve symmetry
		clipper::Spacegroup spgrp1( clipper::Spacegroup::P1 );
		clipper::HKL_info hkl1( spgrp1, hkls.cell(), hkls.resolution(), true );
		clipper::HKL_data<clipper::datatypes::I_sigI<T> > Io1( hkl1 );
		clipper::HKL_data<clipper::datatypes::I_sigI<T> > Ic1( hkl1 );
        clipper::ftype n = 0.0;
		for ( HRI ih = hkl1.first(); !ih.last(); ih.next() ) {
			clipper::datatypes::I_sigI<T> I = I_target[ih.hkl()];
			if ( I.I() >= this->nsig_ * I.sigI() ) {
				Io1[ih] = I;
				Ic1[ih] = I_ref[ih.hkl()];
                n+=1.0;
			}
		}
		// do the aniso scaling
        clipper::BasisFn_aniso_gaussian bfn;
        ctruncate::RestraintFn_sphericalU restU(10.0, n );
        ::TargetFn_scaleLogLikeI1I2<clipper::datatypes::I_sigI<T>,clipper::datatypes::I_sigI<T> >
        tfn( Ic1, Io1 );
		ctruncate::ResolutionFn_nonlinear_rest rfn( hkl1, bfn, tfn, param, restU, 0.0, false  );
        //clipper::ResolutionFn_nonlinear rfn( hkl1, bfn, tfn, param, 0.0, false  );
        
        this->u_i = bfn.u_aniso_orth( rfn.params() );
        this->u_f = 0.5*this->u_i;
        this->bscale_ = bfn.scale( rfn.params() );
		
        /* scale reference set */
        clipper::datatypes::Compute_scale_u_aniso<clipper::datatypes::I_sigI<T> > compute_s(std::sqrt(std::fabs(this->bscale_)), -this->u_f ); //always take u_f
		I_ref.compute( I_ref, compute_s );
        
        /*
        double maxres = hkls.resolution().invresolsq_limit();
        clipper::Spacegroup spg=hkls.spacegroup();
        clipper::Cell cell=hkls.cell();
        std::vector<double> v1(60,0.0),v2(60,0.0),v3(60,0.0),v1a(60,0.0),v2a(60,0.0),v3a(60,0.0);
        for ( HRI ih = I_ref.first(); !ih.last(); ih.next() ) {
            double reso = ih.invresolsq();
            if (!I_ref[ih].missing() && !I_target[ih].missing() ) {
                int bin = int( double(60) * reso /maxres  - 0.001);
                float epsiln = 1.0f/ih.hkl_class().epsilonc();
                for ( int jsym = 0; jsym != spg.num_primitive_symops() ; ++jsym ) {
                    for (int friedal = 0 ; friedal != 2 ; ++friedal) {
                        clipper::HKL ri = int(std::pow( -1.0f, float(friedal) ))*ih.hkl();
                        clipper::HKL rj = ri.transform( spg.primitive_symop( jsym ) );
                        if (clipper::HKL(rj.h(),0,0).invresolsq(cell)/reso < 0.2 ) {
                            v1[bin] += I_ref[ih].I();
                            v1a[bin] += I_target[ih].I();
                        }
                        if (clipper::HKL(0,rj.k(),0).invresolsq(cell)/reso < 0.2 ) {
                            v2[bin] += I_ref[ih].I();
                            v2a[bin] += I_target[ih].I();
                        }
                        if (clipper::HKL(0,0,rj.l()).invresolsq(cell)/reso < 0.2 ){
                            v3[bin] += I_ref[ih].I();
                            v3a[bin] += I_target[ih].I();
                        }                    
                    }
                }
            }
        }
        for (int i=0; i != 60 ; ++i) std::cout << "====>" << v1[i]/v1a[i] << " " << v2[i]/v2a[i] << " " << v3[i]/v3a[i] << std::endl; 
        */
		return true;
	}
    
	template<class T> bool Iscale_llaniso<T>::operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& f_ref, const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& f_target )
	{
		typedef clipper::HKL_info::HKL_reference_index HRI;
		std::vector<clipper::ftype> param( 7, 0.0 );
		const clipper::HKL_info& hkls = f_target.hkl_info();
        
        // precondition
        {
            int nbins = 60; 
            std::vector<clipper::ftype> params_init1( nbins, 1.0 );
            clipper::BasisFn_binner basis_fo_beta1( f_ref, nbins, 1.0 );
            TargetFn_meanInth<clipper::datatypes::F_sigF<T> > target_fo_beta1( f_ref, 1.0);
            clipper::ResolutionFn_nonlinear betaplot1( hkls, basis_fo_beta1, target_fo_beta1, params_init1 ); 
            
            std::vector<clipper::ftype> params_init2( nbins, 1.0 );
            clipper::BasisFn_binner basis_fo_beta2( f_target, nbins, 1.0 );
            TargetFn_meanInth<clipper::datatypes::F_sigF<T> > target_fo_beta2( f_target, 1.0);
            clipper::ResolutionFn_nonlinear betaplot2( hkls, basis_fo_beta2, target_fo_beta2, params_init2 ); 
            
            int nobs = f_target.num_obs();
            std::vector<clipper::ftype> xi, yi, wi; 
            xi.reserve(nobs) ; yi.reserve(nobs) ; wi.reserve(nobs);
            
            for ( HRI ih = f_target.first(); !ih.last(); ih.next() ) {
                if ( !f_target[ih].missing() && betaplot2.f(ih) > 0.0 && betaplot1.f(ih) > 0.00 ) {
                    clipper::ftype lnS = log(betaplot2.f(ih)) - log(betaplot1.f(ih));
                    clipper::ftype res = ih.invresolsq();
                    clipper::ftype eps = ih.hkl_class().epsilon();
                    xi.push_back(res);
                    yi.push_back(lnS);
                    
                    //float weight = ( res < 0.04 ) ? std::pow(0.04/res, 2.0)*eps : eps;  //poorer fit at resolution below 7.5A
                    clipper::ftype weight = 1.0/eps;
                    if (weight > 0.0) {
                        wi.push_back(1.0/weight);
                    } else {
                        wi.push_back(0.0);
                    }
                }
            } 
            clipper::ftype a,b,siga,sigb;
            
            a = b = 0.0f;
            straight_line_fit(xi,yi,wi,nobs,a,b,siga,sigb);
            
            param[0] = b;
            param[1] = param[2] = param[3] = -a;
            
        }

        // expand to P1 in order to preserve symmetry
		clipper::Spacegroup spgrp1( clipper::Spacegroup::P1 );
		clipper::HKL_info hkl1( spgrp1, hkls.cell(), hkls.resolution(), true );
		clipper::HKL_data<clipper::datatypes::F_sigF<T> > fo1( hkl1 );
		clipper::HKL_data<clipper::datatypes::F_sigF<T> >  fc1( hkl1 );
        clipper::ftype n=0.0f;
		for ( HRI ih = hkl1.first(); !ih.last(); ih.next() ) {
			clipper::datatypes::F_sigF<T> f = f_target[ih.hkl()];
			if ( f.f() >= this->nsig_ * f.sigf() ) {
				fo1[ih] = f;
				fc1[ih] = f_ref[ih.hkl()];
                n += 1.0f;
			}
		}
		// do the aniso scaling
        ctruncate::RestraintFn_sphericalU restU(10.0, n );
		clipper::BasisFn_aniso_gaussian bfn;
		::TargetFn_scaleLogLikeF1F2<clipper::datatypes::F_sigF<T>,clipper::datatypes::F_sigF<T> >
		tfn( fc1, fo1 );
		ctruncate::ResolutionFn_nonlinear_rest rfn( hkl1, bfn, tfn, param, restU);
        
        param = rfn.params();
 		this->u_i = bfn.u_aniso_orth( param );
		this->u_f = 0.5 * this->u_i;
		this->bscale_ = bfn.scale( param );
        clipper::datatypes::Compute_scale_u_aniso<clipper::datatypes::F_sigF<T> > compute_s(std::sqrt(std::fabs(this->bscale_)), -this->u_f ); //always take u_f
        f_ref.compute( f_ref, compute_s );

		return true;
	}
	
	template<class T> bool Iscale_llaniso<T>::operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo )
	{
		typedef clipper::datatypes::F_sigF<T> DATA;
		typedef clipper::TargetFn_scaleF1F2<DATA,DATA> TGT1;
		typedef ::TargetFn_scaleLogLikeF1F2<DATA,DATA> TGT2;
        typedef ctruncate::RestraintFn_sphericalU RST2;
		typedef clipper::BasisFn_spline SCALETYPE;
		return scale<DATA,TGT1,TGT2,RST2,SCALETYPE>( fo, -1.0, 12 );
	}
	
	template<class T> bool Iscale_llaniso<T>::operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io )
	{
		typedef clipper::datatypes::I_sigI<T> DATA;
		typedef clipper::TargetFn_scaleI1I2<DATA,DATA> TGT1;
		typedef ::TargetFn_scaleLogLikeI1I2<DATA,DATA> TGT2;
        typedef ctruncate::RestraintFn_sphericalU RST2;
		typedef clipper::BasisFn_spline SCALETYPE;		
		return scale<DATA,TGT1,TGT2,RST2,SCALETYPE>( Io, -1.0, 12 );
	}
	
	template<class T> bool Iscale_llaniso<T>::operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, clipper::ftype resfilter, const int npar_scl )
	{
		typedef clipper::datatypes::F_sigF<T> DATA;
		typedef clipper::TargetFn_scaleF1F2<DATA,DATA> TGT1;
		typedef ::TargetFn_scaleLogLikeF1F2<DATA,DATA> TGT2;
        typedef ctruncate::RestraintFn_sphericalU RST2;
		typedef clipper::BasisFn_spline SCALETYPE;
		return scale<DATA,TGT1,TGT2,RST2,SCALETYPE>( fo, resfilter, npar_scl );
	}
	
	
	template<class T> bool Iscale_llaniso<T>::operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io,clipper:: ftype resfilter, const int npar_scl )
	{
		typedef clipper::datatypes::I_sigI<T> DATA;
		typedef clipper::TargetFn_scaleI1I2<DATA,DATA> TGT1;
		typedef ::TargetFn_scaleLogLikeI1I2<DATA,DATA> TGT2;
        typedef ctruncate::RestraintFn_sphericalU RST2;
		typedef clipper::BasisFn_spline SCALETYPE;
		return scale<DATA,TGT1,TGT2,RST2,SCALETYPE>( Io, resfilter, npar_scl );
	}
	
	
	template<class T> template<class D, class T1, class T2, class R2, class S> bool Iscale_llaniso<T>::scale( clipper::HKL_data<D>& fo, const clipper::ftype resfilter, const int npar_scl )
	{
		typedef clipper::HKL_info::HKL_reference_index HRI;
		// expand to P1 in order to preserve symmetry
		const clipper::HKL_info& hkls = fo.hkl_info();
		clipper::Spacegroup spgrp1( clipper::Spacegroup::P1 );
		clipper::HKL_info hkl1( spgrp1, hkls.cell(), hkls.resolution(), true );
		clipper::HKL_data<D> fo1( hkl1 ), fs1( hkl1 ), fc1( hkl1 );
        clipper::ftype n=0.0;
		for ( HRI ih = hkl1.first(); !ih.last(); ih.next() ) {
			D f = fo[ih.hkl()];
			if ( obs(f) >= this->nsig_ * sigobs(f) ) {
                fo1[ih] = f;
                n += 1.0f;
            }
		}
		
		// calc target values
		fc1 = D( 1.0, 1.0 );
		for ( HRI ih = fc1.first(); !ih.last(); ih.next() )
			fc1[ih].scale( sqrt( ih.hkl_class().epsilon() ) );
		
		// do correction
		this->u_i = this->u_f = clipper::U_aniso_orth( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
		std::vector<clipper::ftype> param( 7, 0.0 );
		clipper::BasisFn_aniso_gaussian bfn;
		clipper::ftype dp;
        R2 restU(10.0, n );
		for ( int c = 0; c != 2; ++c ) {
			// remove current anistropy estimate
			clipper::datatypes::Compute_scale_u_aniso<D> compute_s(1.0,-(this->u_f) );
			fs1.compute( fo1, compute_s );
			
			// and calc E scale
			S basissc( fs1, npar_scl, 1.0 );
			std::vector<clipper::ftype> param_fo( basissc.num_params(), 1.0 );
			T1 tfn_fo( fs1, fc1 );
			clipper::ResolutionFn rfn_fo( hkl1, basissc, tfn_fo, param_fo );
			param_fo = rfn_fo.params();
			
			// prescale F to E-like scale
			for ( HRI ih = fs1.first(); !ih.last(); ih.next() ) {
				fs1[ih] = fo1[ih];
				fs1[ih].scale( sqrt( basissc.f_s( ih.invresolsq(), param_fo ) ) );
			}
						
			// do the aniso scaling
			T2 tfn( fc1, fs1 );
			param = std::vector<clipper::ftype>( 7, 0.0 );
			ctruncate::ResolutionFn_nonlinear_rest rfn( hkl1, bfn, tfn, param, restU );
			param = rfn.params();
			// set trace to zero (i.e. no isotropic correction)
			dp = (param[1]+param[2]+param[3])/3.0;
			param[1] -= dp; param[2] -= dp; param[3] -= dp;
			this->u_i = bfn.u_aniso_orth( param );
			this->u_f = 0.5 * this->u_i;
			//std::cout << c << " | " << param[1] << " " << param[2] << " " << param[3] << " " << param[4] << " " << param[5] << " " << param[6] << "\n"; std::cout << " DP " << dp << "\n";
		}
		//invert the scaling
        for (std::vector<clipper::ftype>::iterator i = param.begin() ; i != param.end() ; ++i ) (*i) *= -1.0;
		// sharpen or smooth as required
		clipper::Matrix<clipper::ftype> m(3,3);
		m(0,0)=       param[1]; m(1,1)=       param[2]; m(2,2)=       param[3];
		m(0,1)=m(1,0)=param[4]; m(0,2)=m(2,0)=param[5]; m(1,2)=m(2,1)=param[6];
		std::vector<clipper::ftype> ev = m.eigen();
		//std::cout << "EIGEN " << param[1] << " " << param[2] << " " << param[3] << " " << ev[0] << " " << ev[1] << " " << ev[2] << std::endl;
		dp = 0.0;
		if ( this->mode_ == Scaling::SHARPEN   ) dp = ev[2];
		if ( this->mode_ == Scaling::UNSHARPEN ) dp = ev[0];
		param[1] -= dp; param[2] -= dp; param[3] -= dp;
		this->u_i = bfn.u_aniso_orth( param );
		this->u_f = 0.5 * this->u_i;
		this->bscale_ = bfn.scale(param);
		// store the results
		clipper::datatypes::Compute_scale_u_aniso<D> compute_s(1.0,-(this->u_f) );
		fo.compute( fo, compute_s );
		
		return true;
    }
    
    //--------Restraints ------------------
    
    const RestraintFn_base::Aderiv& RestraintFn_sphericalU::aderiv(const clipper::Cell& cell, const std::vector<clipper::ftype>&  params) const
    {
        //currently orthogonal approximation
        const int& np = num_params();
        clipper::ftype scale = 1.0/(2.0*sigma_*sigma_);
        clipper::ftype Uiso = (params[1]+params[2]+params[3])/3.0;
        clipper::ftype Adet2 = params[1]*params[2]+params[1]*params[3]+params[2]*params[3];
        clipper::ftype norm2 = std::pow(clipper::Util::twopi(),3.0)*(sigma_*sigma_/Adet2);
        
        clipper::ftype d = std::fabs(params[1]-Uiso)+std::fabs(params[2]-Uiso)+std::fabs(params[3]-Uiso)+2.0*(params[4]+params[5]+params[6]);
        clipper::ftype d2 = std::fabs(params[1]-Uiso)*std::fabs(params[1]-Uiso)
            +std::fabs(params[2]-Uiso)*std::fabs(params[2]-Uiso)
            +std::fabs(params[3]-Uiso)*std::fabs(params[3]-Uiso)  
            +2.0*(params[4]*params[4]+params[5]*params[5]+params[6]*params[6]);
        
        result().f = (norm2 > 0.00001) ? 
            scale*d2 + 0.5*std::log(norm2) :
            scale*d2;
        
        //std::cout << " f= " << result().f << std::endl;
        for (int i = 1; i != 4 ; ++i) result().df[i] = 
            scale*2.0*(params[i]-Uiso);
        for (int i = 4; i != np ; ++i) result().df[i] =
            scale*2.0*params[i];
    
        //for (int i = 1; i != np ; ++i) std::cout << "  df= " << result().df[i] << std::endl;
            
        for (int i = 0; i != np ; ++ i )
            for (int j = 1; j != np ; ++ j )
                result().df2(i,j) = 0.0;
        for (int i = 1; i != 4 ; ++ i ) result().df2(i,i) =
            2.0*scale;
        for (int i = 4; i != np ; ++ i ) result().df2(i,i) =
            2.0*scale;
        /*
        for (int i = 1; i != 4 ; ++i) {
            for (int j = 1; j != 4 ; ++j ) result().df2(i,j) +=
                0.0;
            for (int j = 4; j != np ; ++j ) result().df2(i,j) +=
                0.0;
        }
        for (int i = 4; i != np ; ++i) {
            for (int j = 1; j != 4 ; ++j ) result().df2(i,j) +=
                0.0;
            for (int j = 4; j != np ; ++j ) result().df2(i,j) +=
                0.0;
        }
        */
        result().f *= nobs_;
        for (int i = 1; i != np ; ++i) result().df[i] *= nobs_;
        for (int i = 1; i != np ; ++i) 
            for (int j = 1; j != np ; ++j ) result().df2(i,j) *= nobs_;
        //for (int i = 1; i != np ; ++i)
        //    for (int j = 1; j != np ; ++j) std::cout << "  df2= " << result().df2(i,j) << std::endl;
        
        return result();
    }

    clipper::ftype RestraintFn_sphericalU::f_s( const std::vector<clipper::ftype>& params ) const
    {
        const int& np = num_params();
        clipper::ftype scale = 1.0/(2.0*sigma_*sigma_);
        clipper::ftype Uiso = (params[1]+params[2]+params[3])/3.0;
        clipper::ftype Adet2 = params[1]*params[2]+params[1]*params[3]+params[2]*params[3];
        clipper::ftype norm2 = std::pow(clipper::Util::twopi(),3.0)*(sigma_*sigma_/Adet2);
        
        clipper::ftype d2 = std::fabs(params[1]-Uiso)*std::fabs(params[1]-Uiso)
            +std::fabs(params[2]-Uiso)*std::fabs(params[2]-Uiso)
            +std::fabs(params[3]-Uiso)*std::fabs(params[3]-Uiso)  
            +2.0*(params[4]*params[4]+params[5]*params[5]+params[6]*params[6]);
        
        clipper::ftype f = (norm2 > 0.00001) ? scale*d2 + 0.5*std::log(norm2) : scale*d2;
        f *= nobs_;
        //std::cout << "f_s " << f << std::endl;
        
        return f;
    }
        
    //------------------------------------------------------
    /*! The constructor performs the full minimisation calculation.
     \param hkl_info HKL_info object which provides the reflection list.
     \param basisfn The basis function used to describe the desired property.
     \param targetfn The target function to be minimised.
     \param restraint A restraints function.
     \param damp_ If > 0.0, shifts are fdamped during early cycles to help
     convergence with difficult bases/target conbinations */
    ResolutionFn_nonlinear_rest::ResolutionFn_nonlinear_rest( const clipper::HKL_info& hkl_info, const clipper::BasisFn_base& basisfn, const clipper::TargetFn_base& targetfn, 
        const std::vector<clipper::ftype>& params, 
        const ctruncate::RestraintFn_base& restraint,const clipper::ftype damp, const bool debug )
    {
        clipper::ftype r0, r1, w, scale, g, s, dotprod;
        
        hkl_info_ = &hkl_info;
        basisfn_  = &basisfn;
        targetfn_ = &targetfn;
        restfn_ = &restraint;
        params_ = params;
        cell_ = hkl_info.cell();
        
        int nparams = basisfn_->num_params();
        clipper::Matrix<> dfdp2( nparams, nparams );
        clipper::Matrix<> drdp2( nparams, nparams );
        std::vector<clipper::ftype> drdp( nparams ), dfdp( nparams );
        std::vector<clipper::ftype> shiftn( nparams ), shiftg( nparams );
        std::vector<clipper::ftype> params_old( nparams );
        params_.resize( nparams );
        
        // loop for 20 cycles refining the params
        for ( int n = 0; n != 20; ++n ) {
            params_old = params_;
            
            // calc target fn and derivs
            calc_derivs( params_, r0, drdp, drdp2 );
            
            // include the effect of the restraints on the parameters
            const RestraintFn_base::Aderiv& aderiv = restfn_->aderiv(cell_,params_);
            
            r0 += aderiv.f;
            for (int i = 0; i != nparams ; ++i ) {
                drdp[i] += aderiv.df[i];
                for (int j = 0; j != nparams ; ++j ) {
                    drdp2(i,j) += aderiv.df2(i,j);
                }
            }
            
            // solve for Newton-Raphson shift
            shiftn = drdp2.solve( drdp );
            
            // get step sizes and check direction
            g = s = dotprod = 0.0;
            for ( int k = 0; k != nparams; ++k ) {
                g += drdp[k]*drdp[k];
                s += shiftn[k]*shiftn[k];
                dotprod += drdp[k]*shiftn[k];
            }
            g = sqrt(g);
            s = sqrt(s);
            
            // make gradient shift to match NR shift
            for ( int k = 0; k != nparams; ++k ) shiftg[k] = drdp[k] * s / g;
            
            // if NR shift opposes gradient, then ignore the NR shift
            if ( dotprod < 0.0 ) shiftn = shiftg;
            
            if ( debug ) {
                std::cout << "\nResolution function cycle: " << n << "\n";
                if ( dotprod < 0.0 ) std::cout << " Using scaled grad " << s / g << "\n";
                std::cout << " Gradient " << g << "\n";
                for ( int j = 0; j != nparams; ++j ) {
                    for ( int k = 0; k !=  nparams; ++k ) std::cout << " " << drdp2(j,k);
                    std::cout << "\n";
                }
                for ( int k = 0; k != nparams; ++k ) std::cout << " " << k << " " << params_[k] << " " << drdp[k] << " " << shiftn[k] << "\n";
            }
            
            // now try the step and if necessary reduce the shift
            scale = (1.0+clipper::ftype(n)) / (1.0+clipper::ftype(n)+damp);
            for ( int j = 0; j != 20; ++j ) {
                for ( int i = 0; i != nparams; ++i )
                    params_[i] = params_old[i] - scale*shiftn[i];
                r1 = 0.0;
                for ( clipper::HKL_info::HKL_reference_index ih = hkl_info.first(); !ih.last(); ih.next() ) {
                    clipper::HKL_class cls = ih.hkl_class();
                    w = 2.0 / cls.epsilonc();
                    r1 += w * targetfn.rderiv( ih, f( ih ) ).r;
                }
                r1 += restfn_->f(cell_,params_);
                //std::cout << " sub-cycle(" << j << "): " << r1 << std::endl;
                if ( clipper::Util::is_nan(r1) ) {
                    scale *= 0.5;
                } else {
                    scale *= 0.5;
                    if ( (r1-r0) <= 0.0 ) break;
                }
            }
            
            if ( debug ) std::cout << " Resolution function cycle " << n << " " << r0 << " " << r1 << " " << scale << "\n";
            
            clipper::ftype sig = (r0 > 1.0) ? r0 : 1.0;
            sig *= 10000.0*std::numeric_limits<clipper::ftype>::epsilon();
            if ( std::fabs( r1 - r0 ) < sig ) break;
        }
    }

    /*------------------aniso scaling using wilson approx----------------------*/
    
    template<class T> bool Iscale_wilson_aniso<T>::operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, const clipper::HKL_data<clipper::datatypes::F_phi<T> >& fc )
	{
		typedef clipper::HKL_info::HKL_reference_index HRI;
        typedef clipper::datatypes::F_sigF<T> DATA;
		typedef clipper::TargetFn_scaleF1F2<DATA,DATA > TGT1;
		typedef clipper::TargetFn_scaleF1F2<DATA,DATA > TGT2;
		typedef clipper::BasisFn_spline SCALETYPE;
		// expand to P1 in order to preserve symmetry
		const clipper::HKL_info& hkls = fo.hkl_info();
        
        // wilson element of scaling
        clipper::ftype beta, scaling;
        { 
            int nbins = 60; 
            std::vector<clipper::ftype> params_init1( nbins, 1.0 );
            clipper::BasisFn_binner basis_fo_beta1( fo, nbins, 1.0 );
            TargetFn_meanInth<clipper::datatypes::F_sigF<T> > target_fo_beta1( fo, 1.0);
            clipper::ResolutionFn_nonlinear betaplot1( hkls, basis_fo_beta1, target_fo_beta1, params_init1 ); 
            
            std::vector<clipper::ftype> params_init2( nbins, 1.0 );
            clipper::BasisFn_binner basis_fo_beta2( fc, nbins, 1.0 );
            TargetFn_meanInth<clipper::datatypes::F_phi<T> > target_fo_beta2( fc, 1.0);
            clipper::ResolutionFn_nonlinear betaplot2( hkls, basis_fo_beta2, target_fo_beta2, params_init2 ); 
            
            int nobs = fo.num_obs();
            std::vector<clipper::ftype> xi, yi, wi; 
            xi.reserve(nobs) ; yi.reserve(nobs) ; wi.reserve(nobs);
            
            for ( HRI ih = fc.first(); !ih.last(); ih.next() ) {
                if ( !fo[ih].missing() && betaplot2.f(ih) > 0.0 && betaplot1.f(ih) > 0.0 ) {
                    clipper::ftype lnS = 2.0*log(betaplot2.f(ih)) - 2.0*log(betaplot1.f(ih));                    clipper::ftype res = ih.invresolsq();
                    clipper::ftype eps = ih.hkl_class().epsilon();
                    xi.push_back(res);
                    yi.push_back(lnS);
                    
                    //float weight = ( res < 0.04 ) ? std::pow(0.04/res, 2.0)*eps : eps;  //poorer fit at resolution below 7.5A
                    clipper::ftype weight = 1.0/eps;
                    if (weight > 0.0) {
                        wi.push_back(1.0/weight);
                    } else {
                        wi.push_back(0.0);
                    }
                }
            } 
            clipper::ftype a,b,siga,sigb;
            
            a = b = 0.0f;
            straight_line_fit(xi,yi,wi,nobs,a,b,siga,sigb);
            
            scaling = b;
            beta = 4.0*a;
            
        }

		clipper::HKL_data<clipper::datatypes::F_sigF<T> > fo1( hkls );
		for ( HRI ih = hkls.first(); !ih.last(); ih.next() ) {
            fo1[ih] = fo[ih];
		}
		// do the aniso scaling
		scale<DATA ,TGT1 ,TGT2 ,SCALETYPE >( fo1, -1.0, 12 );
        
        this->bscale_ = std::exp(scaling);
        clipper::U_aniso_orth ua0 = this->u_i;
        clipper::U_aniso_orth ua1(clipper::Util::b2u(beta),clipper::Util::b2u(beta),clipper::Util::b2u(beta),0.0,0.0,0.0);
        this->u_i = ua1+(-1.0*ua0);
        this->u_f = 0.5 * this->u_i;
        
		// store the results
		clipper::datatypes::Compute_scale_u_aniso<DATA> compute_s(std::sqrt(this->kscale()),this->u_f );
		fo.compute( fo, compute_s );
    
		return true;
	}
	
	template<class T> bool Iscale_wilson_aniso<T>::operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io, const clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Ic )
	{
		typedef clipper::HKL_info::HKL_reference_index HRI;
		// expand to P1 in order to preserve symmetry
		const clipper::HKL_info& hkls = Io.hkl_info();
        
        // wilson element of scaling using Io as expectation value
        clipper::ftype beta, scaling;
        {
            int nbins = 60; 
            std::vector<clipper::ftype> params_init1( nbins, 1.0 );
            clipper::BasisFn_binner basis_fo_beta1( Io, nbins, 1.0 );
            TargetFn_meanInth<clipper::datatypes::I_sigI<T> > target_fo_beta1( Io, 1.0);
            clipper::ResolutionFn_nonlinear betaplot1( hkls, basis_fo_beta1, target_fo_beta1, params_init1 ); 
            
            std::vector<clipper::ftype> params_init2( nbins, 1.0 );
            clipper::BasisFn_binner basis_fo_beta2( Ic, nbins, 1.0 );
            TargetFn_meanInth<clipper::datatypes::I_sigI<T> > target_fo_beta2( Ic, 1.0);
            clipper::ResolutionFn_nonlinear betaplot2( hkls, basis_fo_beta2, target_fo_beta2, params_init2 ); 
            
            int nobs = Ic.num_obs();
            std::vector<clipper::ftype> xi, yi, wi; 
            xi.reserve(nobs) ; yi.reserve(nobs) ; wi.reserve(nobs);
            
            for ( HRI ih = Ic.first(); !ih.last(); ih.next() ) {
                if ( !Ic[ih].missing() && betaplot2.f(ih) > 0.0 && betaplot1.f(ih) > 0.0 ) {
                    clipper::ftype lnS = log(betaplot2.f(ih)) - log(betaplot1.f(ih));
                    clipper::ftype res = ih.invresolsq();
                    clipper::ftype eps = ih.hkl_class().epsilon();
                    xi.push_back(res);
                    yi.push_back(lnS);
                    
                    //float weight = ( res < 0.04 ) ? std::pow(0.04/res, 2.0)*eps : eps;  //poorer fit at resolution below 7.5A
                    clipper::ftype weight = 1.0/eps;
                    if (weight > 0.0) {
                        wi.push_back(1.0/weight);
                    } else {
                        wi.push_back(0.0);
                    }
                }
            } 
            clipper::ftype a,b,siga,sigb;
            
            a = b = 0.0f;
            straight_line_fit(xi,yi,wi,nobs,a,b,siga,sigb);
            
            scaling = b;
            beta = 4.0*a;
            
        }
                
		clipper::HKL_data<clipper::datatypes::I_sigI<T> > Ic1( hkls );
		for ( HRI ih = hkls.first(); !ih.last(); ih.next() ) {
            Ic1[ih] = Ic[ih];
		}
		// do the aniso scaling
		typedef clipper::datatypes::I_sigI<T> DATA;
		typedef clipper::TargetFn_scaleI1I2<DATA,DATA> TGT1;
		typedef clipper::TargetFn_scaleI1I2<DATA,DATA> TGT2;
		typedef clipper::BasisFn_spline SCALETYPE;		
		scale<DATA,TGT1,TGT2,SCALETYPE>( Ic1, -1.0, 12 );
        
        this->bscale_ = std::exp(scaling);
        clipper::U_aniso_orth ua0 = this->u_i;
        clipper::U_aniso_orth ua1(clipper::Util::b2u(beta),clipper::Util::b2u(beta),clipper::Util::b2u(beta),0.0,0.0,0.0);
        
        this->u_i = ua1+ua0;
        this->u_f = 0.5 * this->u_i;
        
		// store the results
		clipper::datatypes::Compute_scale_u_aniso<DATA> compute_s(std::sqrt(this->kscale()),this->u_f );
		Io.compute( Io, compute_s );
        
        /*
        double maxres = hkls.resolution().invresolsq_limit();
        clipper::Spacegroup spg=hkls.spacegroup();
        clipper::Cell cell=hkls.cell();
        std::vector<double> v1(60,0.0),v2(60,0.0),v3(60,0.0),v1a(60,0.0),v2a(60,0.0),v3a(60,0.0);
        for ( HRI ih = Io.first(); !ih.last(); ih.next() ) {
            double reso = ih.invresolsq();
            if (!Io[ih].missing() && !Ic[ih].missing() ) {
                int bin = int( double(60) * reso /maxres  - 0.001);
                float epsiln = 1.0f/ih.hkl_class().epsilonc();
                for ( int jsym = 0; jsym != spg.num_primitive_symops() ; ++jsym ) {
                    for (int friedal = 0 ; friedal != 2 ; ++friedal) {
                        clipper::HKL ri = int(std::pow( -1.0f, float(friedal) ))*ih.hkl();
                        clipper::HKL rj = ri.transform( spg.primitive_symop( jsym ) );
                        if (clipper::HKL(rj.h(),0,0).invresolsq(cell)/reso < 0.2 ) {
                            v1[bin] += Io[ih].I();
                            v1a[bin] += Ic[ih].I();
                        }
                        if (clipper::HKL(0,rj.k(),0).invresolsq(cell)/reso < 0.2 ) {
                            v2[bin] += Io[ih].I();
                            v2a[bin] += Ic[ih].I();
                        }
                        if (clipper::HKL(0,0,rj.l()).invresolsq(cell)/reso < 0.2 ){
                            v3[bin] += Io[ih].I();
                            v3a[bin] += Ic[ih].I();
                        }                    }
                }
            }
        }
        for (int i=0; i != 60 ; ++i) std::cout << "====>" << v1[i]/v1a[i] << " " << v2[i]/v2a[i] << " " << v3[i]/v3a[i] << std::endl; 
         */
		return true;
	}
    
	template<class T> bool Iscale_wilson_aniso<T>::operator() ( clipper::HKL_data<clipper::datatypes::F_phi<T> >& fc, const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo )
	{
		typedef clipper::HKL_info::HKL_reference_index HRI;
		// expand to P1 in order to preserve symmetry
		const clipper::HKL_info& hkls = fo.hkl_info();        
        
        // wilson element of scaling
        clipper::ftype beta, scaling;
        {
            int nbins = 60; 
            std::vector<clipper::ftype> params_init1( nbins, 1.0 );
            clipper::BasisFn_binner basis_fo_beta1( fc, nbins, 1.0 );
            TargetFn_meanInth<clipper::datatypes::F_phi<T> > target_fo_beta1( fc, 1.0);
            clipper::ResolutionFn_nonlinear betaplot1( hkls, basis_fo_beta1, target_fo_beta1, params_init1 ); 
            
            std::vector<clipper::ftype> params_init2( nbins, 1.0 );
            clipper::BasisFn_binner basis_fo_beta2( fo, nbins, 1.0 );
            TargetFn_meanInth<clipper::datatypes::F_sigF<T> > target_fo_beta2( fo, 1.0);
            clipper::ResolutionFn_nonlinear betaplot2( hkls, basis_fo_beta2, target_fo_beta2, params_init2 ); 

            int nobs = fo.num_obs();
            std::vector<clipper::ftype> xi, yi, wi; 
            xi.reserve(nobs) ; yi.reserve(nobs) ; wi.reserve(nobs);
            
            for ( HRI ih = fc.first(); !ih.last(); ih.next() ) {
                if ( !fc[ih].missing() && betaplot2.f(ih) > 0.0 && betaplot1.f(ih) > 0.0 ) {
                    clipper::ftype lnS = 2.0*log(betaplot2.f(ih)) - 2.0*log(betaplot1.f(ih));
                    clipper::ftype res = ih.invresolsq();
                    clipper::ftype eps = ih.hkl_class().epsilon();
                    xi.push_back(res);
                    yi.push_back(lnS);
                    
                    //float weight = ( res < 0.04 ) ? std::pow(0.04/res, 2.0)*eps : eps;  //poorer fit at resolution below 7.5A
                    clipper::ftype weight = 1.0/eps;
                    if (weight > 0.0) {
                        wi.push_back(1.0/weight);
                    } else {
                        wi.push_back(0.0);
                    }
                }
            } 
            clipper::ftype a,b,siga,sigb;
            
            a = b = 0.0f;
            straight_line_fit(xi,yi,wi,nobs,a,b,siga,sigb);
            
            scaling = b;
            beta = 4.0*a;
            
        }
        
		clipper::HKL_data<clipper::datatypes::F_sigF<T> > fo1( hkls );
		for ( HRI ih = hkls.first(); !ih.last(); ih.next() ) {
            fo1[ih] = fo[ih];
		}
		// do the aniso scaling
 		typedef clipper::datatypes::F_sigF<T> DATA;
		typedef clipper::TargetFn_scaleF1F2<DATA,DATA> TGT1;
		typedef clipper::TargetFn_scaleF1F2<DATA,DATA> TGT2;
		typedef clipper::BasisFn_spline SCALETYPE;		
		scale<DATA,TGT1,TGT2,SCALETYPE>( fo1, -1.0, 12 );
        
        this->bscale_ = std::exp(scaling); // reverse sign
        clipper::U_aniso_orth ua0 = this->u_i;
        clipper::U_aniso_orth ua1(clipper::Util::b2u(beta),clipper::Util::b2u(beta),clipper::Util::b2u(beta),0.0,0.0,0.0);
        this->u_i = ua1+ua0; //reverse sign in this case
        this->u_f = 0.5 * this->u_i;
        
		// store the results
        clipper::datatypes::Compute_scale_u_aniso<clipper::datatypes::F_phi<T> > compute_s(std::sqrt(this->kscale()),this->u_f );
		fc.compute( fc, compute_s );
        
		return true;
	}
	
  	template<class T> bool Iscale_wilson_aniso<T>::operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, const clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fc )
	{
		typedef clipper::HKL_info::HKL_reference_index HRI;
		// expand to P1 in order to preserve symmetry
		const clipper::HKL_info& hkls = fo.hkl_info();
        
        // wilson element of scaling
        clipper::ftype beta, scaling;
        {
            int nbins = 60; 
            std::vector<clipper::ftype> params_init1( nbins, 1.0 );
            clipper::BasisFn_binner basis_fo_beta1( fo, nbins, 1.0 );
            TargetFn_meanInth<clipper::datatypes::F_sigF<T> > target_fo_beta1( fo, 1.0);
            clipper::ResolutionFn_nonlinear betaplot1( hkls, basis_fo_beta1, target_fo_beta1, params_init1 ); 
            
            std::vector<clipper::ftype> params_init2( nbins, 1.0 );
            clipper::BasisFn_binner basis_fo_beta2( fc, nbins, 1.0 );
            TargetFn_meanInth<clipper::datatypes::F_sigF<T> > target_fo_beta2( fc, 1.0);
            clipper::ResolutionFn_nonlinear betaplot2( hkls, basis_fo_beta2, target_fo_beta2, params_init2 ); 
            
            int nobs = fo.num_obs();
            std::vector<clipper::ftype> xi, yi, wi; 
            xi.reserve(nobs) ; yi.reserve(nobs) ; wi.reserve(nobs);
            
            for ( HRI ih = fc.first(); !ih.last(); ih.next() ) {
                if ( !fo[ih].missing() && betaplot2.f(ih) > 0.0 && betaplot1.f(ih) > 0.0 ) {
                    clipper::ftype lnS = 2.0*log(betaplot2.f(ih)) - 2.0*log(betaplot1.f(ih));
                    clipper::ftype res = ih.invresolsq();
                    clipper::ftype eps = ih.hkl_class().epsilon();
                    xi.push_back(res);
                    yi.push_back(lnS);
                    
                    //float weight = ( res < 0.04 ) ? std::pow(0.04/res, 2.0)*eps : eps;  //poorer fit at resolution below 7.5A
                    clipper::ftype weight = 1.0/eps;
                    if (weight > 0.0) {
                        wi.push_back(1.0/weight);
                    } else {
                        wi.push_back(0.0);
                    }
                }
            } 
            clipper::ftype a,b,siga,sigb;
            
            a = b = 0.0f;
            straight_line_fit(xi,yi,wi,nobs,a,b,siga,sigb);
            
            scaling = b;
            beta = 4.0*a;
            
        }
        
		clipper::HKL_data<clipper::datatypes::F_sigF<T> > fc1( hkls );
		for ( HRI ih = hkls.first(); !ih.last(); ih.next() ) {
            fc1[ih] = fc[ih];
		}
		// do the aniso scaling
		typedef clipper::datatypes::F_sigF<T> DATA;
		typedef clipper::TargetFn_scaleF1F2<DATA,DATA> TGT1;
		typedef clipper::TargetFn_scaleF1F2<DATA,DATA> TGT2;
		typedef clipper::BasisFn_spline SCALETYPE;		
		scale<DATA,TGT1,TGT2,SCALETYPE>( fc1, -1.0, 12 );
        
        this->bscale_ = std::exp(scaling);
        clipper::U_aniso_orth ua0 = this->u_i;
        clipper::U_aniso_orth ua1(clipper::Util::b2u(beta),clipper::Util::b2u(beta),clipper::Util::b2u(beta),0.0,0.0,0.0);
        this->u_i = ua1+ua0;
        this->u_f = 0.5 * this->u_i;
        
		// store the results
		clipper::datatypes::Compute_scale_u_aniso<DATA> compute_s(std::sqrt(this->kscale()),this->u_f );
		fo.compute( fo, compute_s );
        
		return true;
	}


	template<class T> bool Iscale_wilson_aniso<T>::operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo )
	{
		typedef clipper::datatypes::F_sigF<T> DATA;
		typedef clipper::TargetFn_scaleF1F2<DATA,DATA> TGT1;
		typedef clipper::TargetFn_scaleF1F2<DATA,DATA> TGT2;
		typedef clipper::BasisFn_spline SCALETYPE;
		return scale<DATA,TGT1,TGT2,SCALETYPE>( fo, -1.0, 12 );
	}
	
	template<class T> bool Iscale_wilson_aniso<T>::operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io )
	{
		typedef clipper::datatypes::I_sigI<T> DATA;
		typedef clipper::TargetFn_scaleI1I2<DATA,DATA> TGT1;
		typedef clipper::TargetFn_scaleI1I2<DATA,DATA> TGT2;
		typedef clipper::BasisFn_spline SCALETYPE;		
		return scale<DATA,TGT1,TGT2,SCALETYPE>( Io, -1.0, 12 );
	}
	
	template<class T> bool Iscale_wilson_aniso<T>::operator() ( clipper::HKL_data<clipper::datatypes::F_sigF<T> >& fo, clipper::ftype resfilter, const int npar_scl )
	{
		typedef clipper::datatypes::F_sigF<T> DATA;
		typedef clipper::TargetFn_scaleF1F2<DATA,DATA> TGT1;
		typedef clipper::TargetFn_scaleF1F2<DATA,DATA> TGT2;
		typedef clipper::BasisFn_spline SCALETYPE;
		return scale<DATA,TGT1,TGT2,SCALETYPE>( fo, resfilter, npar_scl );
	}
	
	
	template<class T> bool Iscale_wilson_aniso<T>::operator() ( clipper::HKL_data<clipper::datatypes::I_sigI<T> >& Io,clipper:: ftype resfilter, const int npar_scl )
	{
		typedef clipper::datatypes::I_sigI<T> DATA;
		typedef clipper::TargetFn_scaleI1I2<DATA,DATA> TGT1;
		typedef clipper::TargetFn_scaleI1I2<DATA,DATA> TGT2;
		typedef clipper::BasisFn_spline SCALETYPE;
		return scale<DATA,TGT1,TGT2,SCALETYPE>( Io, resfilter, npar_scl );
	}
	
	
	template<class T> template<class D, class T1, class T2, class S> bool Iscale_wilson_aniso<T>::scale( clipper::HKL_data<D>& fo, const clipper::ftype resfilter, const int npar_scl )
	{
		typedef clipper::HKL_info::HKL_reference_index HRI;
		// expand to P1 in order to preserve symmetry
		const clipper::HKL_info& hkls = fo.hkl_info();
		clipper::Spacegroup spgrp1( clipper::Spacegroup::P1 );
		clipper::HKL_info hkl1( spgrp1, hkls.cell(), hkls.resolution(), true );
		clipper::HKL_data<D> fo1( hkl1 ), fs1( hkl1 ), fc1( hkl1 );
		for ( HRI ih = hkl1.first(); !ih.last(); ih.next() ) {
			D f = fo[ih.hkl()];
			if ( obs(f) >= this->nsig_ * sigobs(f) ) fo1[ih] = f;
		}
		
		// calc target values
		fc1 = D( 1.0, 1.0 );
		for ( HRI ih = fc1.first(); !ih.last(); ih.next() )
			fc1[ih].scale( sqrt( ih.hkl_class().epsilon() ) );
		
		// do correction
		this->u_i = this->u_f = clipper::U_aniso_orth( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
		std::vector<clipper::ftype> param( 7, 0.0 );
		clipper::BasisFn_aniso_gaussian bfn;
		clipper::ftype dp;
		for ( int c = 0; c < 2; c++ ) {
			// remove current anistropy estimate
			clipper::datatypes::Compute_scale_u_aniso<D> compute_s(1.0,-(this->u_f) );
			fs1.compute( fo1, compute_s );
            
			// and calc E scale
			S basissc( fs1, npar_scl, 1.0 );
			std::vector<clipper::ftype> param_fo( basissc.num_params(), 1.0 );
			T1 tfn_fo( fs1, fc1 );
			clipper::ResolutionFn rfn_fo( hkl1, basissc, tfn_fo, param_fo );
			param_fo = rfn_fo.params();
			
			// prescale F to E-like scale
			for ( HRI ih = fs1.first(); !ih.last(); ih.next() ) {
				fs1[ih] = fo1[ih];
				fs1[ih].scale( sqrt( basissc.f_s( ih.invresolsq(), param_fo ) ) );
			}
            
			// do the aniso scaling
			T2 tfn( fs1, fc1 );
			param = std::vector<clipper::ftype>( 7, 0.0 );
			clipper::ResolutionFn_nonlinear rfn( hkl1, bfn, tfn, param );
			param = rfn.params();
			// set trace to zero (i.e. no isotropic correction)
			dp = (param[1]+param[2]+param[3])/3.0;
			param[1] -= dp; param[2] -= dp; param[3] -= dp;
			this->u_i = bfn.u_aniso_orth( param );
			this->u_f = 0.5 * this->u_i;
		}
		
		// sharpen or smooth as required
		clipper::Matrix<clipper::ftype> m(3,3);
		m(0,0)=       param[1]; m(1,1)=       param[2]; m(2,2)=       param[3];
		m(0,1)=m(1,0)=param[4]; m(0,2)=m(2,0)=param[5]; m(1,2)=m(2,1)=param[6];
		std::vector<clipper::ftype> ev = m.eigen();
		dp = 0.0;
		if ( this->mode_ == Scaling::SHARPEN   ) dp = ev[2];
		if ( this->mode_ == Scaling::UNSHARPEN ) dp = ev[0];
		param[1] -= dp; param[2] -= dp; param[3] -= dp;
		this->u_i = bfn.u_aniso_orth( param );
		this->u_f = 0.5 * this->u_i;
		this->bscale_ = bfn.scale(param);
		// store the results
		clipper::datatypes::Compute_scale_u_aniso<D> compute_s(1.0,-(this->u_f) );
		fo.compute( fo, compute_s );
		
		return true;
    }

    //------------------------------------------------------
    
	// compile templates
	
    	template class Iscale_aniso_base<clipper::ftype32>;
    	template class Iscale_aniso_base<clipper::ftype64>;
    
	template class Iscale_aniso<clipper::ftype32>;
	template bool Iscale_aniso<clipper::ftype32>::scale<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::TargetFn_scaleF1F2<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::datatypes::F_sigF<clipper::ftype32> >,clipper::TargetFn_scaleLogF1F2<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::datatypes::F_sigF<clipper::ftype32> >,clipper::BasisFn_binner>( clipper::HKL_data<clipper::datatypes::F_sigF<clipper::ftype32> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_aniso<clipper::ftype32>::scale<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::TargetFn_scaleF1F2<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::datatypes::F_sigF<clipper::ftype32> >,clipper::TargetFn_scaleLogF1F2<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::datatypes::F_sigF<clipper::ftype32> >,clipper::BasisFn_linear>( clipper::HKL_data<clipper::datatypes::F_sigF<clipper::ftype32> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_aniso<clipper::ftype32>::scale<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::TargetFn_scaleF1F2<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::datatypes::F_sigF<clipper::ftype32> >,clipper::TargetFn_scaleLogF1F2<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::datatypes::F_sigF<clipper::ftype32> >,clipper::BasisFn_spline>( clipper::HKL_data<clipper::datatypes::F_sigF<clipper::ftype32> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_aniso<clipper::ftype32>::scale<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::TargetFn_scaleI1I2<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::datatypes::I_sigI<clipper::ftype32> >,clipper::TargetFn_scaleLogI1I2<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::datatypes::I_sigI<clipper::ftype32> >,clipper::BasisFn_binner>( clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype32> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_aniso<clipper::ftype32>::scale<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::TargetFn_scaleI1I2<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::datatypes::I_sigI<clipper::ftype32> >,clipper::TargetFn_scaleLogI1I2<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::datatypes::I_sigI<clipper::ftype32> >,clipper::BasisFn_linear>( clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype32> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_aniso<clipper::ftype32>::scale<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::TargetFn_scaleI1I2<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::datatypes::I_sigI<clipper::ftype32> >,clipper::TargetFn_scaleLogI1I2<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::datatypes::I_sigI<clipper::ftype32> >,clipper::BasisFn_spline>( clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype32> >& fo, const clipper::ftype resfilter, const int npar_scl );
	
	template class Iscale_aniso<clipper::ftype64>;
	template bool Iscale_aniso<clipper::ftype64>::scale<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::TargetFn_scaleF1F2<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::datatypes::F_sigF<clipper::ftype64> >,clipper::TargetFn_scaleLogF1F2<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::datatypes::F_sigF<clipper::ftype64> >,clipper::BasisFn_binner>( clipper::HKL_data<clipper::datatypes::F_sigF<clipper::ftype64> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_aniso<clipper::ftype64>::scale<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::TargetFn_scaleF1F2<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::datatypes::F_sigF<clipper::ftype64> >,clipper::TargetFn_scaleLogF1F2<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::datatypes::F_sigF<clipper::ftype64> >,clipper::BasisFn_linear>( clipper::HKL_data<clipper::datatypes::F_sigF<clipper::ftype64> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_aniso<clipper::ftype64>::scale<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::TargetFn_scaleF1F2<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::datatypes::F_sigF<clipper::ftype64> >,clipper::TargetFn_scaleLogF1F2<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::datatypes::F_sigF<clipper::ftype64> >,clipper::BasisFn_spline>( clipper::HKL_data<clipper::datatypes::F_sigF<clipper::ftype64> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_aniso<clipper::ftype64>::scale<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::TargetFn_scaleI1I2<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::datatypes::I_sigI<clipper::ftype64> >,clipper::TargetFn_scaleLogI1I2<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::datatypes::I_sigI<clipper::ftype64> >,clipper::BasisFn_binner>( clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype64> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_aniso<clipper::ftype64>::scale<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::TargetFn_scaleI1I2<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::datatypes::I_sigI<clipper::ftype64> >,clipper::TargetFn_scaleLogI1I2<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::datatypes::I_sigI<clipper::ftype64> >,clipper::BasisFn_linear>( clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype64> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_aniso<clipper::ftype64>::scale<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::TargetFn_scaleI1I2<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::datatypes::I_sigI<clipper::ftype64> >,clipper::TargetFn_scaleLogI1I2<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::datatypes::I_sigI<clipper::ftype64> >,clipper::BasisFn_spline>( clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype64> >& fo, const clipper::ftype resfilter, const int npar_scl );
	
    template class Iscale_llaniso<clipper::ftype32>;
	template bool Iscale_llaniso<clipper::ftype32>::scale<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::TargetFn_scaleF1F2<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::datatypes::F_sigF<clipper::ftype32> >,clipper::TargetFn_scaleLogF1F2<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::datatypes::F_sigF<clipper::ftype32> >, ctruncate::RestraintFn_sphericalU, clipper::BasisFn_binner>( clipper::HKL_data<clipper::datatypes::F_sigF<clipper::ftype32> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_llaniso<clipper::ftype32>::scale<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::TargetFn_scaleF1F2<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::datatypes::F_sigF<clipper::ftype32> >,clipper::TargetFn_scaleLogF1F2<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::datatypes::F_sigF<clipper::ftype32> >, ctruncate::RestraintFn_sphericalU, clipper::BasisFn_linear>( clipper::HKL_data<clipper::datatypes::F_sigF<clipper::ftype32> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_llaniso<clipper::ftype32>::scale<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::TargetFn_scaleF1F2<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::datatypes::F_sigF<clipper::ftype32> >,clipper::TargetFn_scaleLogF1F2<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::datatypes::F_sigF<clipper::ftype32> >, ctruncate::RestraintFn_sphericalU, clipper::BasisFn_spline>( clipper::HKL_data<clipper::datatypes::F_sigF<clipper::ftype32> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_llaniso<clipper::ftype32>::scale<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::TargetFn_scaleI1I2<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::datatypes::I_sigI<clipper::ftype32> >,clipper::TargetFn_scaleLogI1I2<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::datatypes::I_sigI<clipper::ftype32> >, ctruncate::RestraintFn_sphericalU, clipper::BasisFn_binner>( clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype32> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_llaniso<clipper::ftype32>::scale<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::TargetFn_scaleI1I2<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::datatypes::I_sigI<clipper::ftype32> >,clipper::TargetFn_scaleLogI1I2<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::datatypes::I_sigI<clipper::ftype32> >, ctruncate::RestraintFn_sphericalU, clipper::BasisFn_linear>( clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype32> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_llaniso<clipper::ftype32>::scale<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::TargetFn_scaleI1I2<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::datatypes::I_sigI<clipper::ftype32> >,clipper::TargetFn_scaleLogI1I2<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::datatypes::I_sigI<clipper::ftype32> >, ctruncate::RestraintFn_sphericalU, clipper::BasisFn_spline>( clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype32> >& fo, const clipper::ftype resfilter, const int npar_scl );
	
	template class Iscale_llaniso<clipper::ftype64>;
	template bool Iscale_llaniso<clipper::ftype64>::scale<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::TargetFn_scaleF1F2<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::datatypes::F_sigF<clipper::ftype64> >,clipper::TargetFn_scaleLogF1F2<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::datatypes::F_sigF<clipper::ftype64> >, ctruncate::RestraintFn_sphericalU, clipper::BasisFn_binner>( clipper::HKL_data<clipper::datatypes::F_sigF<clipper::ftype64> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_llaniso<clipper::ftype64>::scale<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::TargetFn_scaleF1F2<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::datatypes::F_sigF<clipper::ftype64> >,clipper::TargetFn_scaleLogF1F2<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::datatypes::F_sigF<clipper::ftype64> >, ctruncate::RestraintFn_sphericalU, clipper::BasisFn_linear>( clipper::HKL_data<clipper::datatypes::F_sigF<clipper::ftype64> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_llaniso<clipper::ftype64>::scale<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::TargetFn_scaleF1F2<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::datatypes::F_sigF<clipper::ftype64> >,clipper::TargetFn_scaleLogF1F2<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::datatypes::F_sigF<clipper::ftype64> >,ctruncate::RestraintFn_sphericalU, clipper::BasisFn_spline>( clipper::HKL_data<clipper::datatypes::F_sigF<clipper::ftype64> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_llaniso<clipper::ftype64>::scale<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::TargetFn_scaleI1I2<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::datatypes::I_sigI<clipper::ftype64> >,clipper::TargetFn_scaleLogI1I2<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::datatypes::I_sigI<clipper::ftype64> >, ctruncate::RestraintFn_sphericalU, clipper::BasisFn_binner>( clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype64> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_llaniso<clipper::ftype64>::scale<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::TargetFn_scaleI1I2<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::datatypes::I_sigI<clipper::ftype64> >,clipper::TargetFn_scaleLogI1I2<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::datatypes::I_sigI<clipper::ftype64> >, ctruncate::RestraintFn_sphericalU, clipper::BasisFn_linear>( clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype64> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_llaniso<clipper::ftype64>::scale<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::TargetFn_scaleI1I2<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::datatypes::I_sigI<clipper::ftype64> >,clipper::TargetFn_scaleLogI1I2<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::datatypes::I_sigI<clipper::ftype64> >, ctruncate::RestraintFn_sphericalU, clipper::BasisFn_spline>( clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype64> >& fo, const clipper::ftype resfilter, const int npar_scl );

    template class Iscale_loganiso<clipper::ftype32>;
	template bool Iscale_loganiso<clipper::ftype32>::scale<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::TargetFn_scaleF1F2<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::datatypes::F_sigF<clipper::ftype32> >,clipper::TargetFn_scaleLogF1F2<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::datatypes::F_sigF<clipper::ftype32> >,clipper::BasisFn_binner>( clipper::HKL_data<clipper::datatypes::F_sigF<clipper::ftype32> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_loganiso<clipper::ftype32>::scale<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::TargetFn_scaleF1F2<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::datatypes::F_sigF<clipper::ftype32> >,clipper::TargetFn_scaleLogF1F2<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::datatypes::F_sigF<clipper::ftype32> >,clipper::BasisFn_linear>( clipper::HKL_data<clipper::datatypes::F_sigF<clipper::ftype32> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_loganiso<clipper::ftype32>::scale<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::TargetFn_scaleF1F2<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::datatypes::F_sigF<clipper::ftype32> >,clipper::TargetFn_scaleLogF1F2<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::datatypes::F_sigF<clipper::ftype32> >,clipper::BasisFn_spline>( clipper::HKL_data<clipper::datatypes::F_sigF<clipper::ftype32> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_loganiso<clipper::ftype32>::scale<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::TargetFn_scaleI1I2<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::datatypes::I_sigI<clipper::ftype32> >,clipper::TargetFn_scaleLogI1I2<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::datatypes::I_sigI<clipper::ftype32> >,clipper::BasisFn_binner>( clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype32> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_loganiso<clipper::ftype32>::scale<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::TargetFn_scaleI1I2<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::datatypes::I_sigI<clipper::ftype32> >,clipper::TargetFn_scaleLogI1I2<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::datatypes::I_sigI<clipper::ftype32> >,clipper::BasisFn_linear>( clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype32> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_loganiso<clipper::ftype32>::scale<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::TargetFn_scaleI1I2<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::datatypes::I_sigI<clipper::ftype32> >,clipper::TargetFn_scaleLogI1I2<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::datatypes::I_sigI<clipper::ftype32> >,clipper::BasisFn_spline>( clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype32> >& fo, const clipper::ftype resfilter, const int npar_scl );
	
	template class Iscale_loganiso<clipper::ftype64>;
	template bool Iscale_loganiso<clipper::ftype64>::scale<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::TargetFn_scaleF1F2<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::datatypes::F_sigF<clipper::ftype64> >,clipper::TargetFn_scaleLogF1F2<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::datatypes::F_sigF<clipper::ftype64> >,clipper::BasisFn_binner>( clipper::HKL_data<clipper::datatypes::F_sigF<clipper::ftype64> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_loganiso<clipper::ftype64>::scale<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::TargetFn_scaleF1F2<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::datatypes::F_sigF<clipper::ftype64> >,clipper::TargetFn_scaleLogF1F2<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::datatypes::F_sigF<clipper::ftype64> >,clipper::BasisFn_linear>( clipper::HKL_data<clipper::datatypes::F_sigF<clipper::ftype64> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_loganiso<clipper::ftype64>::scale<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::TargetFn_scaleF1F2<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::datatypes::F_sigF<clipper::ftype64> >,clipper::TargetFn_scaleLogF1F2<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::datatypes::F_sigF<clipper::ftype64> >,clipper::BasisFn_spline>( clipper::HKL_data<clipper::datatypes::F_sigF<clipper::ftype64> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_loganiso<clipper::ftype64>::scale<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::TargetFn_scaleI1I2<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::datatypes::I_sigI<clipper::ftype64> >,clipper::TargetFn_scaleLogI1I2<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::datatypes::I_sigI<clipper::ftype64> >,clipper::BasisFn_binner>( clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype64> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_loganiso<clipper::ftype64>::scale<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::TargetFn_scaleI1I2<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::datatypes::I_sigI<clipper::ftype64> >,clipper::TargetFn_scaleLogI1I2<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::datatypes::I_sigI<clipper::ftype64> >,clipper::BasisFn_linear>( clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype64> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_loganiso<clipper::ftype64>::scale<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::TargetFn_scaleI1I2<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::datatypes::I_sigI<clipper::ftype64> >,clipper::TargetFn_scaleLogI1I2<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::datatypes::I_sigI<clipper::ftype64> >,clipper::BasisFn_spline>( clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype64> >& fo, const clipper::ftype resfilter, const int npar_scl );

    template class Iscale_wilson_aniso<clipper::ftype32>;
	template bool Iscale_wilson_aniso<clipper::ftype32>::scale<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::TargetFn_scaleF1F2<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::datatypes::F_sigF<clipper::ftype32> >,clipper::TargetFn_scaleLogF1F2<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::datatypes::F_sigF<clipper::ftype32> >,clipper::BasisFn_binner>( clipper::HKL_data<clipper::datatypes::F_sigF<clipper::ftype32> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_wilson_aniso<clipper::ftype32>::scale<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::TargetFn_scaleF1F2<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::datatypes::F_sigF<clipper::ftype32> >,clipper::TargetFn_scaleLogF1F2<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::datatypes::F_sigF<clipper::ftype32> >,clipper::BasisFn_linear>( clipper::HKL_data<clipper::datatypes::F_sigF<clipper::ftype32> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_wilson_aniso<clipper::ftype32>::scale<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::TargetFn_scaleF1F2<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::datatypes::F_sigF<clipper::ftype32> >,clipper::TargetFn_scaleLogF1F2<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::datatypes::F_sigF<clipper::ftype32> >,clipper::BasisFn_spline>( clipper::HKL_data<clipper::datatypes::F_sigF<clipper::ftype32> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_wilson_aniso<clipper::ftype32>::scale<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::TargetFn_scaleI1I2<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::datatypes::I_sigI<clipper::ftype32> >,clipper::TargetFn_scaleLogI1I2<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::datatypes::I_sigI<clipper::ftype32> >,clipper::BasisFn_binner>( clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype32> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_wilson_aniso<clipper::ftype32>::scale<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::TargetFn_scaleI1I2<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::datatypes::I_sigI<clipper::ftype32> >,clipper::TargetFn_scaleLogI1I2<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::datatypes::I_sigI<clipper::ftype32> >,clipper::BasisFn_linear>( clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype32> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_wilson_aniso<clipper::ftype32>::scale<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::TargetFn_scaleI1I2<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::datatypes::I_sigI<clipper::ftype32> >,clipper::TargetFn_scaleLogI1I2<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::datatypes::I_sigI<clipper::ftype32> >,clipper::BasisFn_spline>( clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype32> >& fo, const clipper::ftype resfilter, const int npar_scl );
	
	template class Iscale_wilson_aniso<clipper::ftype64>;
	template bool Iscale_wilson_aniso<clipper::ftype64>::scale<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::TargetFn_scaleF1F2<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::datatypes::F_sigF<clipper::ftype64> >,clipper::TargetFn_scaleLogF1F2<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::datatypes::F_sigF<clipper::ftype64> >,clipper::BasisFn_binner>( clipper::HKL_data<clipper::datatypes::F_sigF<clipper::ftype64> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_wilson_aniso<clipper::ftype64>::scale<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::TargetFn_scaleF1F2<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::datatypes::F_sigF<clipper::ftype64> >,clipper::TargetFn_scaleLogF1F2<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::datatypes::F_sigF<clipper::ftype64> >,clipper::BasisFn_linear>( clipper::HKL_data<clipper::datatypes::F_sigF<clipper::ftype64> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_wilson_aniso<clipper::ftype64>::scale<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::TargetFn_scaleF1F2<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::datatypes::F_sigF<clipper::ftype64> >,clipper::TargetFn_scaleLogF1F2<clipper::datatypes::F_sigF<clipper::ftype64>,clipper::datatypes::F_sigF<clipper::ftype64> >,clipper::BasisFn_spline>( clipper::HKL_data<clipper::datatypes::F_sigF<clipper::ftype64> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_wilson_aniso<clipper::ftype64>::scale<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::TargetFn_scaleI1I2<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::datatypes::I_sigI<clipper::ftype64> >,clipper::TargetFn_scaleLogI1I2<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::datatypes::I_sigI<clipper::ftype64> >,clipper::BasisFn_binner>( clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype64> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_wilson_aniso<clipper::ftype64>::scale<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::TargetFn_scaleI1I2<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::datatypes::I_sigI<clipper::ftype64> >,clipper::TargetFn_scaleLogI1I2<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::datatypes::I_sigI<clipper::ftype64> >,clipper::BasisFn_linear>( clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype64> >& fo, const clipper::ftype resfilter, const int npar_scl );
	template bool Iscale_wilson_aniso<clipper::ftype64>::scale<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::TargetFn_scaleI1I2<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::datatypes::I_sigI<clipper::ftype64> >,clipper::TargetFn_scaleLogI1I2<clipper::datatypes::I_sigI<clipper::ftype64>,clipper::datatypes::I_sigI<clipper::ftype64> >,clipper::BasisFn_spline>( clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype64> >& fo, const clipper::ftype resfilter, const int npar_scl );
} // namespace clipper
