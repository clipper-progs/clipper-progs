#ifndef CTRUNCATE_UTILS_H
#define CTRUNCATE_UTILS_H

#include "clipper/clipper.h"
#include "clipper/clipper-contrib.h"
#include "clipper/clipper-ccp4.h"
#include "clipper/clipper-minimol.h"

#include <vector>
#include <cmath>

void tricart(clipper::Cell& cell, clipper::Mat33<float>& transf);
void MatrixToString( clipper::Mat33<int>& twinoper, clipper::String &s );
void straight_line_fit(std::vector<float>& x, std::vector<float>& y, std::vector<float>& w, int n, float &a, float &b, float &siga, float &sigb);


namespace ctruncate {
	
	// Close(a,b[,tol])  true if a == b within tolerance
	template<class T1, class T2> inline static bool
	Close(const T1& a, const T2& b,
		  const T1& tol=1.0e-6)
	{ return std::abs(a-b)<=tol; }
	
	class IceRing
	{
	public:
		// Resolution in A, width in 1/d^2 units
		IceRing(const double& Resolution, const double& width);
		
		// Resolution of ring: return centre of ring as d* = 1/d
		double Dstar() const {return std::sqrt(ring_invressqr);}
		
		// If in ring, returns true
		bool InRing(const double& invresolsq) const;
		// Clear intensity sums
		void ClearSums();
		// Add in IsigI, use clipper version for now rather than pointless version
		void AddObs(const clipper::datatypes::I_sigI<float>& I_sigI, const double& invresolsq);
		
		void SetReject() {reject=true;}
		void SetReject(const bool& Rej) {reject=Rej;}
		
		// Results
		double MeanI() const;
		double MeanSigI() const;
		double MeanSSqr() const;
		int N() const {return nI;}
		
		// Reject flag, true means reject reflection in this range
		bool Reject() const {return reject;}
		
	private:
		double ring_invressqr;  // centre of ring in 1/d^2
		double halfwidth_invressqr;  // halfwidth of ring in 1/d^2
		double sum_I;
		double sum_sigI;
		double sum_sSqr;
		int nI;
		bool reject;
	};
	//--------------------------------------------------------------
	class Rings
	{
	public:
		Rings() : nrings(0) {}
		
		// Resolution in A, width in 1/d^2 units
		void AddRing(const double& Resolution, const double& width);
		
		// Set up default ice rings: 3.90, 3.67, 3.44A
		void DefaultIceRings();
		
		// Clear list
		void Clear();
		
		int Nrings() const {return nrings;}
		
		// Resolution of iring'th ring as d*
		double Dstar(const int& iring) const
		{return rings.at(iring).Dstar();}
		
		// Copy rejected rings only
		void CopyRejRings(const Rings& other);
		
		// If in ring, returns ring number (0,n-1), else = -1 
		int InRing(const double& invresolsq) const;
		// Clear intensity sums
		void ClearSums();
		// Add in IsigI
		void AddObs(const int& Iring, const clipper::datatypes::I_sigI<float>& I_sigI, const double& invresolsq);
		
		void SetReject(const int& Iring);
		void SetReject(const int& Iring, const bool& Rej);
		
		// Results
		double MeanI(const int& Iring) const;
		double MeanSigI(const int& Iring) const;
		double MeanSSqr(const int& Iring) const;
		// Reject flag, true means reject reflection in this range
		bool Reject(const int& Iring) const; 
		int N(const int& Iring) const;
		
	private:
		int nrings;
		std::vector<IceRing> rings;
		
		void CheckRing(const int& Iring) const;
	};
	
	//--------------------------------------------------------------
	
	class PattPeak 
	{
	public:
		PattPeak( float maxinvres, int nbins = 20, float temp = -1.0f );
		
		float operator() (clipper::BasisFn_spline& basis_fo, clipper::ResolutionFn& resol_fn);
		void setBins(int nbins) { _nbins = nbins; _patterson.resize(nbins); }
		void setMaxInvRes(float maxinvres) { _maxinvres = maxinvres; }
		void setWidth(float width) { _width = width; }
		
		float p000() { return _P000; }
		float sigma() { return _sigma; }
		float optRes();
				   
	private:
		void calcOriginPeak(clipper::BasisFn_spline& basis_fo, clipper::ResolutionFn& resol_fn);
		void fitOriginPeak(clipper::BasisFn_spline& basis_fo, clipper::ResolutionFn& resol_fn);
		
		float _P000;
		float _sigma;
		
		int _nbins; // number of bins
		float _maxinvres;  // maximum 1/res to use
		float _width;
		std::vector<float> _patterson;
	};
	
	class Matthews
	{
	public:
		Matthews(bool protein = true, bool rna = false);
		
		int operator() (clipper::Cell& cell, clipper::Spacegroup& spacegroup, clipper::SEQfile& file, float reso = 99.0f );
		int operator() (clipper::Cell& cell, clipper::Spacegroup& spacegroup, int nresidues, float reso = 99.0f );
		
		float matthews(clipper::Cell& cell, clipper::Spacegroup& spacegroup, int nmol);
		float matthews(int i) { return _cmath[i]; }
		
		float solvent(int i) { return _solvent[i]; }
		
		float prob(int i) { return _prob[i]; }
		
		void summary();
		
	private:
		// prob calc
		float vmProb(double x, double y0, double xc, double wt, double a, double s) 
		{
			double z = (x-xc)/wt;
			return y0+a*(std::exp(-std::exp(-z)-z*s+1.0));
		}
		int calcProt(clipper::Cell& cell, clipper::Spacegroup& spacegroup, float reso);
		int calcComp(clipper::Cell& cell, clipper::Spacegroup& spacegroup, float reso);
		int calcDNA(clipper::Cell& cell, clipper::Spacegroup& spacegroup, float reso);
		
		/* array of contents based on mmdb single letter codes */
		static const char _ResidueName1[21];
		static const int _Catoms[21];
		static const int _Hatoms[21];
		static const int _Natoms[21];
		static const int _Oatoms[21];
		static const int _Satoms[21];
		
		/* Arrays   :       1: rbin, 2: p0, 3: vmbar, 4: w, 5: a, 6: s
		   columns :        1-13 protein data inclusive to corresponding binmax value
		                    14 DNA (all resolutions)
		                    15 DNA/Protein complexes(25%/75%) */
		static const double _rbin[15];
		static const double _p0[15];
		static const double _vmbar[15];
		static const double _wcoeff[15];
		static const double _acoeff[15];
		static const double _scoeff[15];
		
		static const float _densProt;
		static const float _densDNA;

		bool _protein;
		bool _rna;
		float _weight;
		float _density;
		float _reso;
		std::vector<float> _cmath;
		std::vector<float> _solvent;
		std::vector<float> _prob;
		std::vector<float> _probt;
	};
}
#endif
