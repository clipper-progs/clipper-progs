//
//     ctruncate_utils.cpp
//     Copyright (C) 2006-2008 Norman Stein
//
//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Application.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//

#include "ctruncate_utils.h"
#include <cstdlib>

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



void straight_line_fit(std::vector<clipper::ftype>& x, std::vector<clipper::ftype>& y, std::vector<clipper::ftype>& w, int n, clipper::ftype &a, clipper::ftype &b, clipper::ftype &siga, clipper::ftype &sigb)
{
  // fits a straight line through a set of points (yi,xi) using least squares
    clipper::ftype d;
    clipper::ftype sx,sy,sxx,sxy,sw;
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

void tricart(clipper::Cell& cell, clipper::Mat33<clipper::ftype>& transf)
{
	/* Calculates the matrix that transforms coordinates relative to the
	 triclinic axes a1, a2 and a3 to a Cartesian set of axes. a2(cart)
	 runs along a2, a1(cart) lies in the plane of a1 and a2 and a3(cart)
	 runs along a1 x a2.
	 
	 I.e. X || b* x c,  Y || b*,  Z || a* x b* || c
	 This does not agree with any of the standard orthogonalisation
	 codes, and so we cannot use library functions to get transf.
	 
	 Alpha, beta, gamma must be given in radians.
	 Lit.: M.G.Rossman & D.M.Blow, Acta Cryst.(1962) Vol.15,24
	 formula (9).
	 */
    float c1,c2,c3,s1,s3,cw,sw;
	
    c1 = cos( cell.alpha_star() );
    c2 = cos( cell.beta_star() );
    c3 = cos( cell.gamma_star() );
    s1 = sin( cell.alpha_star() );
    s3 = sin( cell.gamma_star() );
	
    cw = (c2-c1*c3)/(s1*s3);
    //sw = sin(acos(cw));   // ?? use simpler formula
	sw = 0.0;
	if (fabs(cw) < 1.0) sw = sqrt(1.0-cw*cw);
	
    transf(0,0) = cell.a_star()*s3*sw;
    transf(0,1) = 0.0;
    transf(0,2) = 0.0;
    transf(1,0) = cell.a_star()*c3;
    transf(1,1) = cell.b_star();
    transf(1,2) = cell.c_star()*c1;
    transf(2,0) = cell.a_star()*s3*cw;
    transf(2,1) = 0.0;
    transf(2,2) = cell.c_star()*s1;
	return;
}



// convert twinning operator from a matrix to a string
// code modified from Symop::format

void MatrixToString( clipper::Mat33<int>& op, clipper::String &s )
{
	clipper::String t, hkl="hkl";
	for ( int i = 0; i < 3; i++ ) {
		t = "";
		for ( int j = 0; j < 3; j++ ) {
			if ( op(i,j) != 0 ) {
				t += ( op(i,j) > 0 ) ? "+" : "-";
				if ( std::abs( op(i,j) ) != 12 )
					t += clipper::String::rational( fabs( float( op(i,j) )/12.0 ), 24 );
				t += hkl[j];
			}
		}
		s += t.substr( ( t[0] == '+' ) ? 1 : 0 );
		if ( i < 2 ) s+= ", ";
	}
}



#include <assert.h>
#define ASSERT assert
namespace ctruncate
{
	
	
	//--------------------------------------------------------------
	void Rings::DefaultIceRings()
	{
		// Set up default ice rings: 3.90, 3.67, 3.44A
		// clear any existing ones first
		nrings = 0;
		rings.clear();
		// Revised figures from Garman and Schneider
		// use a constant width in reciprocal space
		// rings should get wider at higher resolution, but they
		// probably get weaker as well 
		const double RWIDTH = 0.005;
		// resolution in A, full width in d* 1/A
		AddRing(3.8996, RWIDTH);
		AddRing(3.6697, RWIDTH);
		AddRing(3.4398, RWIDTH);
		AddRing(2.6699, RWIDTH);
		AddRing(2.2499, RWIDTH);
		AddRing(2.0800, RWIDTH);
		AddRing(1.9499, RWIDTH);
		AddRing(1.9200, RWIDTH);
		AddRing(1.8900, RWIDTH);
		AddRing(1.7250, RWIDTH);
	}
	//--------------------------------------------------------------
	void Rings::CheckRing(const int& Iring) const
	{
		ASSERT (Iring < nrings && Iring >= 0);
	}
	//--------------------------------------------------------------
	// Resolution in A, width in 1/d^2 units
	void Rings::AddRing(const double& Resolution, const double& width)
	{
		rings.push_back(IceRing(Resolution, width));
		nrings++;
	}
	//--------------------------------------------------------------
	// Copy rejected rings only
	void Rings::CopyRejRings(const Rings& other)
	{
		rings.clear();
		for (int i=0;i<other.nrings;i++)
		{
			if (other.rings[i].Reject())
			{rings.push_back(other.rings[i]);}
		}
		nrings = rings.size();
	}
	//--------------------------------------------------------------
	// Clear list
	void Rings::Clear()
	{
		nrings = 0;
		rings.clear();
	}
	//--------------------------------------------------------------
	// If in ring, returns ring number (0,n-1), else = -1 
	int Rings::InRing(const double& invresolsq) const
	{
		for (int i=0;i<rings.size();i++)
		{
			if (rings[i].InRing(invresolsq))
			{return i;}
		}
		return -1;
	}
	//--------------------------------------------------------------
	void Rings::ClearSums()
	{
		for (int i=0;i<rings.size();i++)
		{
			rings[i].ClearSums();
		}
	}
	//--------------------------------------------------------------
	void Rings::AddObs(const int& Iring, const clipper::datatypes::I_sigI<float>& I_sigI,
					   const double& invresolsq)
	{
		CheckRing(Iring);
		rings[Iring].AddObs(I_sigI, invresolsq);
	}
	//--------------------------------------------------------------
	void Rings::SetReject(const int& Iring)
	{
		CheckRing(Iring);
		rings[Iring].SetReject();
	}
	//--------------------------------------------------------------
	void Rings::SetReject(const int& Iring, const bool& Rej)
	{
		CheckRing(Iring);
		rings[Iring].SetReject(Rej);
	}
	//--------------------------------------------------------------
	bool Rings::Reject(const int& Iring) const
	{
		CheckRing(Iring);
		return rings[Iring].Reject();
	}
	//--------------------------------------------------------------
	double Rings::MeanI(const int& Iring) const
	{
		CheckRing(Iring);
		return rings[Iring].MeanI();
	}
	//--------------------------------------------------------------
	double Rings::MeanSigI(const int& Iring) const
	{
		CheckRing(Iring);
		return rings[Iring].MeanSigI();
	}
	//--------------------------------------------------------------
	double Rings::MeanSSqr(const int& Iring) const
	{
		CheckRing(Iring);
		return rings[Iring].MeanSSqr();
	}
	//--------------------------------------------------------------
	int Rings::N(const int& Iring) const
	{
		CheckRing(Iring);
		return rings[Iring].N();
	}
	//--------------------------------------------------------------
	//--------------------------------------------------------------
	IceRing::IceRing(const double& Resolution, const double& width)
	{
		ring_invressqr = 1./(Resolution*Resolution);
		halfwidth_invressqr = 0.5*width;
		ClearSums();
	}
	//--------------------------------------------------------------
	bool IceRing::InRing(const double& invresolsq) const
	{
		if (Close<double,double>(invresolsq,
								 ring_invressqr,halfwidth_invressqr))
		{return true;}
		return false;
	}
	//--------------------------------------------------------------
	void IceRing::ClearSums()
	{
		sum_I = 0.0;
		sum_sigI = 0.0;
		sum_sSqr = 0.0;
		nI = 0;
		reject = false;
	}
	//--------------------------------------------------------------
	void IceRing::AddObs(const clipper::datatypes::I_sigI<float>& I_sigI, const double& invresolsq)
	{
		sum_I += I_sigI.I();
		sum_sigI += I_sigI.sigI();
		sum_sSqr += invresolsq;
		nI++;
	}
	//--------------------------------------------------------------
	double IceRing::MeanI() const
	{
		if (nI > 0)
		{return sum_I/nI;}
		else
		{return 0.0;}
	}
	//--------------------------------------------------------------
	double IceRing::MeanSigI() const
	{
		if (nI > 0)
		{return sum_sigI/nI;}
		else
		{return 0.0;}
	}
	//--------------------------------------------------------------
	double IceRing::MeanSSqr() const
	{
		if (nI > 0)
		{return sum_sSqr/nI;}
		else
		{return 0.0;}
	}

}

