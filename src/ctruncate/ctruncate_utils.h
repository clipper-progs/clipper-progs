#ifndef CTRUNCATE_UTILS_H
#define CTRUNCATE_UTILS_H

#include "clipper/clipper.h"
#include "clipper/clipper-contrib.h"
#include "clipper/clipper-ccp4.h"
#include "clipper/clipper-minimol.h"

#include <vector>
#include <cmath>

void tricart(clipper::Cell& cell, clipper::Mat33<clipper::ftype>& transf);
void MatrixToString( clipper::Mat33<int>& twinoper, clipper::String &s );
void straight_line_fit(std::vector<clipper::ftype>& x, std::vector<clipper::ftype>& y, std::vector<clipper::ftype>& w, int n, clipper::ftype &a, clipper::ftype &b, clipper::ftype &siga, clipper::ftype &sigb);


namespace ctruncate {
	
	//! class for program start and end.  Extend for some html
	class CCP4Program
	{
	public:
		CCP4Program( const char* name, const char* vers, const char* rcsdate );
		~CCP4Program();
		void summary_beg() const;
		void summary_end() const;
		void set_termination_message( std::string msg ) { msg_ = msg; }
		std::stringstream& xml_start(std::stringstream&) const;
		std::stringstream& xml_end(std::stringstream&) const;
	private:
		bool html, summ;
		std::string name_, msg_, vers_, rundate_;
	};
	
	//warning, cfile_ may go out of scope
	class ReflectionFile
	{
	public:
		ReflectionFile() : cfile_(NULL) {}
		ReflectionFile(clipper::CCP4MTZfile& cfile, std::string& cfilename) : cfile_(&cfile), cfilename_(cfilename)
		{ cell_ = cfile_->cell(); spgr_ = cfile_->spacegroup(); cdname_ = cfile_->assigned_paths()[0].notail(); }
		~ReflectionFile() {}
		
		void operator()(clipper::CCP4MTZfile& cfile, std::string& cfilename) 
		{ cfile_ = &cfile; cfilename_ = cfilename; cell_ = cfile_->cell(); spgr_ = cfile_->spacegroup(); cdname_ = cfile_->assigned_paths()[0].notail(); return; }
		
		void output() const;
		std::stringstream& xml_output(std::stringstream&) const;
		
	private:
		clipper::CCP4MTZfile* cfile_;
		std::string cfilename_;
		std::string cdname_;
		clipper::Cell cell_;
		clipper::Spacegroup spgr_;
	};
	
	//warning, cdata_ may go out of scope
	class ReflectionData
	{
	public:
		ReflectionData() : cdata_(NULL) {}
		ReflectionData(clipper::HKL_data_base& cdata) : cdata_(&cdata) {}
		~ReflectionData() {}
		
		void operator()(clipper::HKL_data_base& cdata) { cdata_ = &cdata; return; }
		//wrap HKL_data_base functions
		clipper::String type() const { return cdata_->type(); }
		int num_obs() const;
		int num_reflections() const { return cdata_->num_obs(); }
		int num_centric() const;
		int num_acentric() const;
		clipper::Range<clipper::ftype> invresolsq_range() const { return cdata_->invresolsq_range(); }
        clipper::HKL max_hkl() const;
        clipper::HKL max_sym_hkl() const;
	
		void output() const;
		std::stringstream& xml_output(std::stringstream&) const;
		
	private:
		clipper::HKL_data_base* cdata_;
	};
	
	class Utils
	{
		
	public:
		//! Convert +/-/cov to mean , with NaN checks and weights
		template<class T> inline static T mean( const T& pl, const T& mi, const T& spl, const T& smi )
		{
			if ( clipper::Util::is_nan(pl) ) return mi;
			else if (clipper::Util::is_nan(mi) ) return pl;
			else {
				T wpl = T(1.0)/(spl*spl);
				T wmi = T(1.0)/(smi*smi);
				return (wpl*pl+wmi*mi)/(wpl+wmi);
			}
			return clipper::Util::nan();
		}
		//! Convert sig+/sig-/cov to sig , with NaN checks and weights
		template<class T> inline static T sig_mean( const T& pl, const T& mi, const T& spl, const T& smi, const T& cov )
		{
			if ( clipper::Util::is_nan(pl) ) return smi;
			else if (clipper::Util::is_nan(mi) ) return spl;
			else if (clipper::Util::is_nan(cov) ) {
				T wpl = T(1.0)/(spl*spl);
				T wmi = T(1.0)/(smi*smi);
				T mean = (wpl*pl+wmi*mi)/(wpl+wmi);
				T sig1 = std::sqrt(1.0/(wpl+wmi));
                                return sig1;
			} else {
				T wpl = T(1.0)/(spl*spl);
				T wmi = T(1.0)/(smi*smi);
				T mean = (wpl*pl+wmi*mi)/(wpl+wmi);
				T sig1 = std::sqrt(1.0/(wpl+wmi+2.0*cov));
                                return sig1;
			}
			return clipper::Util::nan();
		}

		static double pbdv(double v, double x);
	};
	
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
		template <class T> void AddObs(const clipper::datatypes::I_sigI<T>& I_sigI, const double& invresolsq,const double& multiplicity=1.0) {
			if (!I_sigI.missing() ) {
				sum_I += multiplicity*I_sigI.I();
				sum_sigI += multiplicity*I_sigI.sigI();
				sum_sSqr += multiplicity*invresolsq;
				nI+=multiplicity;
			}
			nO += multiplicity;
		}		
		template <class T> void AddObs(const clipper::datatypes::F_sigF<T>& I_sigI, const double& invresolsq,const double& multiplicity=1.0) {
			if (!I_sigI.missing() ) {
				sum_I += multiplicity*I_sigI.f();
				sum_sigI += multiplicity*I_sigI.sigf();
				sum_sSqr += multiplicity*invresolsq;
				nI+=multiplicity;
			}
			nO += multiplicity;
		}
		void SetReject() {reject=true;}
		void SetReject(const bool& Rej) {reject=Rej;}
		
		// Results
		double MeanI() const;
		double MeanSigI() const;
		double MeanSSqr() const;
		double N() const {return nI;}
		double Comp() const;
		
		// Reject flag, true means reject reflection in this range
		bool Reject() const {return reject;}
		
	private:
		double ring_invressqr;  // centre of ring in 1/d^2
		double halfwidth_invressqr;  // halfwidth of ring in 1/d^2
		double sum_I;
		double sum_sigI;
		double sum_sSqr;
		double nI;
		double nO;
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
	        
                // Copy rings
                void Copy(const Rings& other);
	
		// Copy rejected rings only
		void CopyRejRings(const Rings& other);
		
		// If in ring, returns ring number (0,n-1), else = -1 
		int InRing(const double& invresolsq) const;
		// Clear intensity sums
		void ClearSums();
		// Add in IsigI
		template <class T> void AddObs(const int& Iring, const clipper::datatypes::I_sigI<T>& I_sigI, const double& invresolsq, const double& multiplicity=1.0) {
			CheckRing(Iring);
			rings[Iring].AddObs(I_sigI, invresolsq, multiplicity );
		}
		template <class T> void AddObs(const int& Iring, const clipper::datatypes::F_sigF<T>& I_sigI, const double& invresolsq, const double& multiplicity=1.0) {
			CheckRing(Iring);
			rings[Iring].AddObs(I_sigI, invresolsq, multiplicity );
		}
		
		void SetReject(const int& Iring);
		void SetReject(const int& Iring, const bool& Rej);
		
		// Results
		double MeanI(const int& Iring) const;
		double MeanSigI(const int& Iring) const;
		double MeanSSqr(const int& Iring) const;
		// Reject flag, true means reject reflection in this range
		bool Reject(const int& Iring) const; 
		double N(const int& Iring) const;
		double Comp(const int& Iring) const;
		
	private:
		int nrings;
		std::vector<IceRing> rings;
		
		void CheckRing(const int& Iring) const;
	};
	
	
	//! Resolution ordinal gernerator
	/*! This class is a helper class for functions which need to divide
	 reflections up by resolution whilst guaranteeing a certain
	 distribution of number of reflections per range. It takes a list
	 of reflections, one at a time, and calculates a function to get
	 the approximate ordinal number of a reflection in a list sorted by
	 resolution.
	 */
	class ResolutionRange_ordinal : public clipper::Generic_ordinal
	{
	public:
		//! initialiser: takes an HKL_info and uses all reflections.
		template<class T> void init( const clipper::HKL_info& hklinfo, const clipper::Range<T>&, const clipper::ftype& power );
		//! initialiser: takes an HKL_data & uses non-missing reflections.
		template<class T> void init( const clipper::HKL_data_base& hkldata, const clipper::Range<T>&, const clipper::ftype& power );
		//! initialiser: takes an HKL_data + Cell & uses non-missing reflections.
		template<class T> void init( const clipper::HKL_data_base& hkldata, const clipper::Cell& cell, const clipper::Range<T>&, const clipper::ftype& power );
	};
	
	// resolution ordinal
	
	template <class T> void ResolutionRange_ordinal::init( const clipper::HKL_info& hklinfo, const clipper::Range<T>& range, const clipper::ftype& power )
	{
		clipper::Generic_ordinal::init( range, 1000 );
		for (clipper::HKL_info::HKL_reference_index ih = hklinfo.first(); !ih.last(); ih.next() ) 
			if (range.contains(ih.invresolsq() ) ) accumulate( ih.invresolsq() );
		prep_ordinal();
		
		for ( int i = 0; i < hist.size(); i++ )
			hist[i] = pow( hist[i], 1.0/power );
	}
	
	template<class T> void ResolutionRange_ordinal::init( const clipper::HKL_data_base& hkldata, const clipper::Range<T>& range, const clipper::ftype& power )
	{
		Generic_ordinal::init( range, 1000 );
		for (clipper::HKL_info::HKL_reference_index ih = hkldata.first_data(); !ih.last(); hkldata.next_data(ih) )
			if (range.contains(ih.invresolsq() ) ) accumulate( ih.invresolsq() );
		prep_ordinal();
		
		for ( int i = 0; i < hist.size(); i++ )
			hist[i] = pow( hist[i], 1.0/power );
	}
	
	template<class T> void ResolutionRange_ordinal::init( const clipper::HKL_data_base& hkldata, const clipper::Cell& cell, const clipper::Range<T>& range, const clipper::ftype& power )
	{
		Generic_ordinal::init( range, 1000 );
		for (clipper::HKL_info::HKL_reference_index ih = hkldata.first_data(); !ih.last(); hkldata.next_data(ih) )
			if (range.contains(ih.invresolsq() ) ) accumulate( ih.hkl().invresolsq( cell ) );
		prep_ordinal();
		
		for ( int i = 0; i < hist.size(); i++ )
			hist[i] = pow( hist[i], 1.0/power );
	}
	
	
}
#endif
