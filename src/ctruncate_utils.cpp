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



void straight_line_fit(std::vector<float>& x, std::vector<float>& y, std::vector<float>& w, int n, float &a, float &b, float &siga, float &sigb)
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

void tricart(clipper::Cell& cell, clipper::Mat33<float>& transf)
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
				if ( abs( op(i,j) ) != 12 )
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
	//--------------------------------------------------------------
	
	
	
	PattPeak::PattPeak(float maxinvres, int nbins, float temp ) : _maxinvres(maxinvres), _nbins(nbins), _patterson(nbins,0.0f)
	{
		float coef = 1.5f;
		float dmax = (sqrt(1.0f/maxinvres)/3.0f)*2.0f*1.5f;
		float btemp = ( temp > 0.0f ) ? temp : 0.0f ;
		float dmax1 =  std::sqrt(btemp/(clipper::Util::twopi()*clipper::Util::twopi()*2.0f) )*2.0f*1.5f ;
		_width = ( dmax1 > dmax ) ? dmax1 : dmax ;
	}
	
	float PattPeak::operator()(clipper::BasisFn_spline& basis_fo, clipper::ResolutionFn& resol_fn)
	{
		calcOriginPeak(basis_fo, resol_fn);
		fitOriginPeak(basis_fo, resol_fn);
		
		return optRes();
	}
	
	float PattPeak::optRes()
	{
		float width_res = 0.715*1.0f/_maxinvres;
		float width_patt = 2.0f*_sigma;
		
		return std::sqrt((width_patt*width_patt+width_res*width_res)/2.0f);
	}
	
	// calculate the patterson origin peak in 1-d using fitted data
	
	void PattPeak::calcOriginPeak(clipper::BasisFn_spline& basis_fo, clipper::ResolutionFn& resol_fn) 
	/* -----------------------------------------------------------
	 
	 <I(s)> = average intensity , dS = DETRH
	 
	 P(r) = 4pi * Int <I> * (sin(2pisr)/(2pisr)) * s^2 * ds
	 
	 -----------------------------------------------------------*/
	{
		float widthd = _width/float(_nbins);
		float widthr = _maxinvres/float(_nbins);

		for (int id = 0 ; id != _nbins ; ++id ) {
			float d = widthd*(float(id)+0.5f);
			
			for ( int ir=0; ir!=_nbins; ++ir ) {
				float res = widthr*(float(ir)+0.5f); 
				float rsq = res*res;
				float intensity = basis_fo.f_s( rsq, resol_fn.params() );
				float sr = clipper::Util::twopi()*res*d;
				
				_patterson[id] += 2.0f*intensity * std::sin(sr)*rsq*widthr/(res*d);			}
		}
		return;
	}
	
	void PattPeak::fitOriginPeak(clipper::BasisFn_spline& basis_fo, clipper::ResolutionFn& resol_fn) 
	/* fit gaussain to OriginPeak
	 */
	{
		std::vector<float> weights(_nbins,1.0f);
		std::vector<float> x(_nbins);
		std::vector<float> y(_nbins);
		
		for( int i = 0 ; i != _nbins ; ++i) {
			float dist = (float(i)+0.5f)*_width/float(_nbins);
			x[i] = 0.25*dist*dist;
			y[i] = std::log(1.0f/_patterson[i]);
		}
		
		float a, b, siga, sigb;
		
		straight_line_fit(x,y,weights,_nbins,a,b,siga,sigb);
		
		_P000 = std::exp(-a);
		
		float btemp = 0.25*b;
		
		_sigma = std::sqrt(1.0f/std::abs(2.0f*btemp) );
		
		return;	
	}
	
	//                       A   R   N   D   C   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
	const char Matthews::_ResidueName1[21] = {'A','R','N','D','C','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'};
	const int Matthews::_Catoms[21] = { 3,  6,  4,  4,  3,  3,  5,  5,  2,  6,  6,  6,  6,  5,  9,  5,  3,  4, 11,  9,  5 };
	const int Matthews::_Hatoms[21] = { 7, 14,  8,  7,  7,  7, 10,  9,  5,  9, 13, 13, 14, 11, 11,  9,  7,  9, 12, 11, 11 };
	const int Matthews::_Natoms[21] = { 1,  4,  2,  1,  1,  1,  2,  1,  1,  3,  1,  1,  2,  1,  1,  1,  1,  1,  2,  1,  1 };
	const int Matthews::_Oatoms[21] = { 2,  2,  3,  4,  2,  2,  3,  4,  2,  2,  2,  2,  2,  2,  2,  2,  3,  3,  2,  3,  2 };
	const int Matthews::_Satoms[21] = { 0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0 };
	const double Matthews::_rbin[15] = {1.199,1.501,1.650,1.801,1.901,2.001,2.201,2.401, 
		2.601,2.801,3.101,3.501,5.001,5.001,5.001};
	const double Matthews::_p0[15] = {0.085,0.312,0.400,0.503,0.597,0.729,1.052,1.781,
		2.852,3.386,3.841,4.281,4.592,1.503,0.257};
	const double Matthews::_vmbar[15] = {2.052,2.102,2.122,2.132,2.140,2.155,2.171,2.182,
		2.191,2.192,2.205,2.211,2.210,2.256,2.324};
	const double Matthews::_wcoeff[15] = {0.213,0.214,0.220,0.225,0.226,0.231,0.236,0.241,
		0.242,0.244,0.244,0.244,0.245,0.446,0.327};
	const double Matthews::_acoeff[15] = {28.38,102.7,187.5,339.3,434.1,540.5,686.2,767.9,
		835.9,856.9,854.0,849.6,846.7,136.6,47.10};
	const double Matthews::_scoeff[15] = {0.953,0.807,0.775,0.702,0.648,0.640,0.635,0.589,
		0.584,0.542,0.500,0.485,0.480,1.180,0.466};
	
	const float Matthews::_densProt = (1.0f/0.74f);
	const float Matthews::_densDNA = (1.0f/0.5f);
	
	Matthews::Matthews(bool protein, bool rna) : _protein(protein), _rna(rna)
	{
		if (_protein ) {
			if (_rna ) {
				_density = 0.25*_densDNA + 0.75*_densProt; //assume 25%/75% DNA/protein, as in Kantardjieff
			} else {
				_density = _densProt;
			}
		} else {
			_density = _densDNA;
		}
	}
	
	int Matthews::operator() (clipper::Cell& cell, clipper::Spacegroup& spacegroup, clipper::SEQfile& file, float reso)
	{
		clipper::MMoleculeSequence seq;
		file.import_molecule_sequence( seq );
		clipper::MPolymerSequence poly = seq[0];
		clipper::String sequence = poly.sequence();
		
		std::vector<int> numatoms(5,0);
		for (int i=0; i<sequence.length(); i++) {
			for (int j=0; j<21; j++) {
				if (sequence[i] == _ResidueName1[j]) {
					numatoms[0] += _Catoms[j];
					numatoms[1] += _Natoms[j];
					numatoms[2] += _Oatoms[j];
					numatoms[3] += _Hatoms[j];
					numatoms[4] += _Satoms[j];
					break;
				}
			}
		}
		_weight = 12.0*numatoms[0] + 14.0*numatoms[1] + 16.0*numatoms[2] +1.0*numatoms[3] +32.0*numatoms[4];
		
		int nmol = 0;
		
		if (_protein ) {
			if (_rna ) {
				nmol = calcComp(cell, spacegroup, reso);
			} else {
				nmol = calcProt(cell, spacegroup, reso);
			}
		} else {
			nmol = calcDNA(cell, spacegroup, reso);
		}
		return nmol;
		
	}
	
	int Matthews::operator() (clipper::Cell& cell, clipper::Spacegroup& spacegroup, int nresidues, float reso)
	{
		_reso = reso;
		// use average weight of C G A T residue (NOT base pair)
		float dna_weight = 325.96*nresidues;
		
		float prot_weight = (12.0*5 + 14.0*1.35 + 16.0*1.5 +1.0*8 +32.0*0.05)*nresidues;
		
		int nmol = 0;
		
		if (_protein ) {
			if (_rna ) {
				_weight = 0.25*dna_weight + 0.75*prot_weight; //assume 25%/75% DNA/protein, as in Kantardjieff
				nmol = calcComp(cell, spacegroup, reso);
			} else {
				_weight = prot_weight;
				nmol = calcProt(cell, spacegroup, reso);
			}
		} else {
			_weight = dna_weight;
			nmol = calcDNA(cell, spacegroup, reso);
		}
		return nmol;
	}
	
	int Matthews::calcProt(clipper::Cell& cell, clipper::Spacegroup& spacegroup, float reso)
	{
		float volume = cell.volume();
		float nsym = spacegroup.num_symops();
		int maxmols1 = std::floor(0.602f*volume*(_density)/(_weight*nsym) )+1;
		int nmol = 0;
		_reso = reso;
		
		float maxprob = 0.0f;
		
		_cmath.resize(maxmols1);
		_solvent.resize(maxmols1);
		_prob.resize(maxmols1);
		_probt.resize(maxmols1);
		
		for (int i = 0 ; i != maxmols1 ; ++i ) {
			_prob[i] = 0.0f;
			_probt[i] = 0.0f;
		}
		
		int imols = 1;
		
		for ( ; imols != maxmols1 ; ++imols ) {
			_cmath[imols] = volume/(_weight*float(imols)*nsym);
			float solvent = (1.0-1.0f/(0.602f*_cmath[imols]*_density))*100.0;
			if (_solvent[imols] < 0.0f ) break;
			_solvent[imols] = solvent;
			
			for (int j = 0; j != 13 ; ++j ) {
				_probt[imols] += vmProb(_cmath[imols],_p0[j],_vmbar[j],_wcoeff[j],_acoeff[j],_scoeff[j]);
			}
			_probt[0] += _probt[imols];
			if ( _probt[imols] >= maxprob) {
				nmol = imols;
				maxprob = _probt[imols];
			}
		}
		
		if ( reso < 99.0f ) {
			maxprob = 0.0f;
			int ireso = 0;
			for ( ; ireso != 13 ; ++ireso ) {
				if ( reso <= _rbin[ireso] ) break;
			}
			for (int i = 1 ; i != imols ; ++i ) {
				_prob[i] = vmProb(_cmath[i],_p0[ireso],_vmbar[ireso],_wcoeff[ireso],_acoeff[ireso],_scoeff[ireso]);
				_prob[0] += _prob[i];
				if ( _prob[i] >= maxprob) {
					nmol = i;
					maxprob = _prob[i];
				}				
			}
			for (int i = 1 ; i != imols ; ++i ) _prob[i]/=_prob[0];
		}
		for (int i = 1 ; i != imols ; ++i ) _probt[i]/=_probt[0];	
		
		return nmol;
	}
	
	
	int Matthews::calcComp(clipper::Cell& cell, clipper::Spacegroup& spacegroup, float reso)
	{
		float volume = cell.volume();
		float nsym = spacegroup.num_symops();
		int maxmols1 = std::floor(0.602f*volume*(_density)/(0.602f*_weight*nsym) )+1;
		int nmol = 0;
		int COMP = 13;
		
		float maxprob = 0.0f;
		
		_cmath.resize(maxmols1);
		_solvent.resize(maxmols1);
		_prob.resize(maxmols1);
		_probt.resize(maxmols1);
		
		for (int i = 1 ; i != maxmols1 ; ++i ) {
			_cmath[i] = volume/(_weight*i*nsym);
			_solvent[i] = (1.0-1.0f/(0.602f*_cmath[i]*_density))*100.0;

			_prob[i] = vmProb(_cmath[i],_p0[COMP],_vmbar[COMP],_wcoeff[COMP],_acoeff[COMP],_scoeff[COMP]);
			_prob[0] += _prob[i];
			if ( _prob[i] >= maxprob) {
				nmol = i;
				maxprob = _probt[i];
			}				
		}
		for (int i = 1 ; i != maxmols1 ; ++i ) _prob[i]/=_prob[0];
		
		return nmol;
	}
	
	
	
	int Matthews::calcDNA(clipper::Cell& cell, clipper::Spacegroup& spacegroup, float reso)
	{
		float volume = cell.volume();
		float nsym = spacegroup.num_symops();
		int maxmols1 = std::floor(0.602f*volume*(_density)/(0.602f*_weight*nsym) )+1;
		int nmol = 0;
		int DNA = 14;
		
		float maxprob = 0.0f;
		
		_cmath.resize(maxmols1);
		_solvent.resize(maxmols1);
		_prob.resize(maxmols1);
		_probt.resize(maxmols1);
		
		for (int i = 1 ; i != maxmols1 ; ++i ) {
			_cmath[i] = volume/(_weight*i*nsym);
			_solvent[i] = (1.0-1.0f/(0.602f*_cmath[i]*_density))*100.0;
			_prob[i] = vmProb(_cmath[i],_p0[DNA],_vmbar[DNA],_wcoeff[DNA],_acoeff[DNA],_scoeff[DNA]);
			_prob[0] += _prob[i];
			if ( _prob[i] >= maxprob) {
				nmol = i;
				maxprob = _probt[i];
			}				
		}
		for (int i = 1 ; i != maxmols1 ; ++i ) _prob[i]/=_prob[0];
		
		return nmol;
	}
	

	
	float Matthews::matthews(clipper::Cell& cell, clipper::Spacegroup& spacegroup, int nmol)
	{
		/*  Assume 64% of cell volume solvent in dna crystal, with density 1.0/0.5
			Assume 60% of cell volume solvent in dna/protein crystal..
			Assume 47% of cell volume solvent in protein crytal, with density 1.0/0.74
			*/
		float volume = cell.volume();
		float nsym = spacegroup.num_symops();
			
		float dna_weight = (1.0-0.64)*_densDNA*(0.602*volume)/(nmol*nsym);
		float prot_weight = (1.0-0.47)*_densProt*(0.602*volume)/(nmol*nsym);
		if (_protein ) {
			if (_rna ) {
				_weight = 0.25*dna_weight + 0.75*prot_weight; //assume 25%/75% DNA/protein, as in Kantardjieff
			} else {
				_weight = prot_weight;
			}
		} else {
			_weight = dna_weight;
		}
		return volume/(_weight*nmol*nsym);
	}
		
	void Matthews::summary(){
		int maxmol1 = _cmath.size();
		printf("$TABLE: Matthews coefficients:\n");
		printf("$GRAPHS");
		if (_protein ) {
			if (_rna ) {
				printf(": 25%% DNA/75%% protein, as in Kantardjieff");
			} else {
				printf(": Protein crystal");
			}
		} else {
			printf(": DNA crystal");
		}
		if ( _reso < 99.0f ) {
			printf(" computed at resolution of %5.3f :A:1,2,3,4\n", _reso,maxmol1);
			printf("$$ Nmol/asym Matthews_Coeff sovlent_frac P(%5.3f) $$\n$$\n",_reso);
			
			for (int i=1; i != maxmol1 ; ++i) {
				printf("  %3d       %6.2f          %6.2f       %6.2f\n", i, _cmath[i], _solvent[i]/100.0f, _prob[i]);
			}
		} else {
			printf(" :A:1,2,3,4 \n", maxmol1);
			printf("$$ Nmol/asym Matthews Coeff solvent_frac P() $$\n$$\n");
			
			for (int i=1; i != maxmol1 ; ++i) {
				printf("  %3d       %6.2f          %6.2f       %6.2f\n", i, _cmath[i], _solvent[i]/100.0f, _probt[i]);
			}
		}
		printf("$$\n\n");
		
	}
}

