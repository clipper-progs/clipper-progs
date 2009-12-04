//
//     ctruncate_utils.cpp
//     Copyright (C) 2006-2008 Norman Stein
//
//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Application.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//

#include "clipper/ccp4/ccp4_mtz_io.h"
#include <algorithm>
#include <iostream>
#include <vector>

extern "C" {
#include <math.h>
#include <stdio.h>
}

int bisect(double (*f)(double), double x1, double x2, double &xmid)
{
	double epsilon = 1.0e-5;
	double f1 = (*f)(x1);
	double f2 = (*f)(x2);
	double fmid;
	if ( f1*f2 > 0.0 ) {
		printf("Bisect: root not bracketed\n");
		return(0);
	}

	for (int i=0; i<50; i++) {
		xmid = 0.5*(x1+x2);
		fmid = (*f)(xmid);
		if ( fabs(fmid) < epsilon ) return(1); 
		if ( f1*fmid < 0.0 ) {
			x2 = xmid;
			f2 = fmid;
		}
		else {
			x1 = xmid;
			f1 = fmid;
		}
	}
	printf("Bisect: too many iterations\n");
	return(0);
}



void straight_line_fit(std::vector<float> x, std::vector<float> y, std::vector<float> w, int n, float &a, float &b, float &siga, float &sigb)
{
  // fits a straight line through a set of points (yi,xi) using least squares
	float d;
	float sx,sy,sxx,sxy,sw;
	int i;
	sx = 0;
	sy = 0;
	sw = 0;
	sxx = 0;
	sxy = 0;
	for (i=0;i<n;i++) {
		sxx += w[i]*x[i]*x[i];
        sx += w[i]*x[i];
		sy += w[i]*y[i];
		sw += w[i];
		sxy += w[i]*x[i]*y[i];
	}
	d = sxx*sw - sx*sx;
	//printf("%e %e %e %e %e %e\n", sxx,sx,sy,sw,sxy,d);
	if ( fabs(d) < 1.0e-3 * fabs(sxx*sw) ) {
		clipper::Message::message( clipper::Message_fatal( "least squares fit: zero denominator" ) );
		return;
	}
	a = (sxy*sw - sx*sy)/d;
	b = (sy*sxx - sx*sxy)/d;
	siga = sqrt(sw/d);
	sigb = sqrt(sxx/d);
	return;
}



