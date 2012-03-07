//
//     CTRUNCATE
//     Copyright (C) 2006-2011 Norman Stein, Charles Ballard
//
//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Application.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//

#ifndef __CTRUNCATE_MOMENTS_H
#define __CTRUNCATE_MOMENTS_H

#include "clipper/clipper.h"
#include "clipper/clipper-ccp4.h"

namespace ctruncate {

	void moments_Z(clipper::HKL_data<clipper::data32::I_sigI>& isig, float maxres, int nbins);
	
	void moments_E(clipper::HKL_data<clipper::data32::I_sigI>& isig, float maxres, int nbins, CCP4Program& prog);
}

#endif

