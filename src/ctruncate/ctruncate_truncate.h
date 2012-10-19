//
//     CTRUNCATE
//     Copyright (C) 2006-2011 Norman Stein, Charles Ballard
//
//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Application.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//


#ifndef __CTRUNCATE_TRUNCATE_H
#define __CTRUNCATE_TRUNCATE_H

#include "alt_hkl_datatypes.h"
#include "clipper/clipper.h"
#include "clipper/clipper-ccp4.h"
#include "ccp4/csymlib.h"

namespace ctruncate {
	
	int truncate(clipper::HKL_data<clipper::data32::I_sigI>& isig, clipper::HKL_data<clipper::data32::I_sigI>& jsig,
				 clipper::HKL_data<clipper::data32::F_sigF>& fsig, clipper::ResolutionFn& Sigma, float scalef, 
				 CSym::CCP4SPG *spg1, clipper::Resolution& reso, int& nrej, bool debug);
	int truncate(clipper::HKL_data<clipper::data32::J_sigJ_ano>& isig, clipper::HKL_data<clipper::data32::J_sigJ_ano>& jsig, 
				 clipper::HKL_data<clipper::data32::G_sigG_ano>& fsig, clipper::ResolutionFn& Sigma, float scalef, 
				 CSym::CCP4SPG *spg1, clipper::Resolution& reso, int& nrej, bool debug);
	int truncate(clipper::HKL_data<clipper::data32::I_sigI>& isig, clipper::HKL_data<clipper::data32::I_sigI>& jsig,
				 clipper::HKL_data<clipper::data32::F_sigF>& fsig, clipper::HKL_data<clipper::data32::I_sigI>& Sigma, float scalef, 
				 CSym::CCP4SPG *spg1, clipper::Resolution& reso, int& nrej, bool debug);
	int truncate(clipper::HKL_data<clipper::data32::J_sigJ_ano>& isig, clipper::HKL_data<clipper::data32::J_sigJ_ano>& jsig, 
				 clipper::HKL_data<clipper::data32::G_sigG_ano>& fsig, clipper::HKL_data<clipper::data32::I_sigI>& Sigma, float scalef, 
				 CSym::CCP4SPG *spg1, clipper::Resolution& reso, int& nrej, bool debug);

	int truncate(clipper::HKL_data<clipper::data32::I_sigI>& isig, clipper::HKL_data<clipper::data32::I_sigI>& jsig,
				 clipper::HKL_data<clipper::data32::F_sigF>& fsig, float scalef, 
				 CSym::CCP4SPG *spg1, clipper::Resolution& reso, int& nrej, bool debug);
	int truncate(clipper::HKL_data<clipper::data32::J_sigJ_ano>& isig, clipper::HKL_data<clipper::data32::J_sigJ_ano>& jsig, 
				 clipper::HKL_data<clipper::data32::G_sigG_ano>& fsig, float scalef, 
				 CSym::CCP4SPG *spg1, clipper::Resolution& reso, int& nrej, bool debug);
}

#endif
