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
#include <ccp4/ccp4_general.h>
#include <ccp4/ccp4_program.h>
#include <cstdlib>
#include <iomanip>

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
	
	CCP4Program::CCP4Program( const char* name, const char* vers, const char* rcsdate )
	{
		name_ = name;
		vers_ = vers;
		rundate_ = rcsdate;
		html = ( getenv( "CCP_SUPPRESS_HTML" ) == NULL );
		summ = ( getenv( "CCP_SUPPRESS_SUMMARY" ) == NULL );
		CCP4::ccp4ProgramName( (char*)name );
		CCP4::ccp4_prog_vers( (char*)vers );
		CCP4::ccp4RCSDate( (char*)rcsdate );
		summary_beg();
		if ( html ) std::cout << "<html> <!-- CCP4 HTML LOGFILE -->" << std::endl
			<< "<hr>" << std::endl << "<pre>" << std::endl;
		CCP4::ccp4_banner();
		summary_end();
		CCP4::ccp4ProgramTime(1);
	}
	
	
	CCP4Program::~CCP4Program()
	{
		std::cout << std::endl;
		summary_beg();
		std::cout << name_ << ": " << msg_ << std::endl;
		CCP4::ccp4ProgramTime(0);
		if ( html ) std::cout << "</pre>" << std::endl << "</html>" << std::endl;
		summary_end();
	}
	
	
	void CCP4Program::summary_beg() const
	{
		if ( summ ) {
			if ( html )
				std::cout << "<B><FONT COLOR='#FF0000'><!--SUMMARY_BEGIN-->" << std::endl;
			else
				std::cout << "<!--SUMMARY_BEGIN-->" << std::endl;
		}
	}
	
	
	void CCP4Program::summary_end() const
	{
		if ( summ ) {
			if ( html )
				std::cout << "<!--SUMMARY_END--></FONT></B>" << std::endl;
			else
				std::cout << "<!--SUMMARY_END-->" << std::endl;
		}
	}
	
	std::stringstream& CCP4Program::xml_start(std::stringstream& ss) const 
	{
                std::string name = name_;
                std::use_facet<std::ctype<char> >(std::locale()).toupper(&name[0], &name[0] + name.size());
		ss << "<" << name << " version=\"" << vers_ <<  "\" RunTime=\"" << rundate_ << "\" >" << std::endl;
		return ss;
	}
	
	std::stringstream& CCP4Program::xml_end(std::stringstream& ss) const 
	{
                std::string name = name_;
                std::use_facet<std::ctype<char> >(std::locale()).toupper(&name[0], &name[0] + name.size());
		ss << "</" << name << ">" << std::endl;
		return ss;
	}
	
	//--------------------------------------------------------------
	
	void ReflectionFile::output() const
	{
		printf("\nReflection File INFO:\n\n");
		std::cout << "Reflection file name: " << cfilename_ << std::endl;
		std::cout << "Crystal/dataset names: " << cdname_ << "\n"; 
		printf("Spacegroup: %s (number %4d)\n", spgr_.symbol_hm().c_str(), spgr_.spacegroup_number() );
		printf("Cell parameters: %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n", cell_.a(), cell_.b(), cell_.c(), 
			   clipper::Util::rad2d( cell_.alpha() ), clipper::Util::rad2d( cell_.beta() ), clipper::Util::rad2d( cell_.gamma() ) );
	}
	
	std::stringstream& ReflectionFile::xml_output(std::stringstream& ss) const
	/*
	<ReflectionFile stream="HKLIN" name="/tmp/ccb/test-results/aucn.mtz">
	<cell>
	<a>  88.91</a>
	<b>  88.91</b>
	<c>  229.2</c>
	<alpha>     90</alpha>
	<beta>     90</beta>
	<gamma>     90</gamma>
	</cell>
	<MergedData>False</MergedData>
	<SpacegroupName> P 41 2 2</SpacegroupName>
	</ReflectionFile>
	*/
	{
		ss << "<ReflectionFile name=\"" << cfilename_ << " \">" << std::endl;
		ss << "<CrystalDatasetId>" << cdname_ << "</CrystalDatasetId>" << std::endl;
		ss << "<cell>" << std::endl;
		ss << "<a>" << std::fixed << std::setprecision(2) << cell_.a() << "</a>" << std::endl;
		ss << "<b>" << std::fixed << std::setprecision(2) << cell_.b() << "</b>" << std::endl;
		ss << "<c>" << std::fixed << std::setprecision(2) << cell_.c() << "</c>" << std::endl;
		ss << "<alpha>" << std::fixed << std::setprecision(2) << clipper::Util::rad2d(cell_.alpha() ) << "</alpha>" << std::endl;
		ss << "<beta>" << std::fixed << std::setprecision(2) << clipper::Util::rad2d(cell_.beta() ) << "</beta>" << std::endl;
		ss << "<gamma>" << std::fixed << std::setprecision(2) << clipper::Util::rad2d(cell_.gamma() ) << "</gamma>" << std::endl;
		ss << "</cell>" << std::endl;
		ss << "<SpacegroupName>" << spgr_.symbol_hm() << "</SpacegroupName>" << std::endl;
		ss << "</ReflectionFile>" << std::endl;
		return ss;
	}
	
	
	//--------------------------------------------------------------
	
	
	/*! \return The number of centric data in the object. */
	int ReflectionData::num_centric() const
	{
		int num = 0;
		for ( clipper::HKL_info::HKL_reference_index ih = cdata_->first_data(); !ih.last(); cdata_->next_data(ih) ) if (ih.hkl_class().centric() ) ++num;
		return num;
	}
	
	/*! \return The number of acentric data in the object. */
	int ReflectionData::num_acentric() const
	{
		int num = 0;
		for ( clipper::HKL_info::HKL_reference_index ih = cdata_->first_data(); !ih.last(); cdata_->next_data(ih) ) if (!ih.hkl_class().centric() ) ++num;
		return num;
	}

	/*! \return The number of reflections in the object. */
	int ReflectionData::num_obs() const
	{
		int num = 0;
		clipper::ftype a,s;
		clipper::xtype working[cdata_->data_size()];
		for ( clipper::HKL_info::HKL_reference_index ih = cdata_->first_data(); !ih.last(); cdata_->next_data(ih) ) {
			if (cdata_->data_size() == 2) {
				a = working[0];
				s = working[1];
				if (!clipper::Util::is_nan(a) && ! s >= 0.0 ) ++num;
			} else {
				a = working[0];
				s = working[1];
				if (!clipper::Util::is_nan(a) && ! s >= 0.0 ) ++num;
				a = working[2];
				s = working[3];
				if (!clipper::Util::is_nan(a) && ! s >= 0.0 ) ++num;	}
		}
		return num;
	}
	
    /*! \return limit of HKL values */
    clipper::HKL ReflectionData::max_hkl() const
    {
        int h(0), k(0), l(0);
        for ( clipper::HKL_data_base::HKL_reference_index ih = cdata_->first_data(); !ih.last(); cdata_->next_data(ih) ) {
            clipper::HKL t = ih.hkl();
            h = std::max(std::abs(t.h()),h);
            k = std::max(std::abs(t.k()),k);
            l = std::max(std::abs(t.l()),l);
        }
        return clipper::HKL(h,k,l);
    }
    
    /*! \return limit of HKL values */
    clipper::HKL ReflectionData::max_sym_hkl() const
    {
        int h(0), k(0), l(0);
        clipper::Spacegroup spgr(cdata_->hkl_info().spacegroup() );
        for ( clipper::HKL_data_base::HKL_reference_index ih = cdata_->first_data(); !ih.last(); cdata_->next_data(ih) ) {
            clipper::HKL t(ih.hkl() );
            h = std::max(std::abs(t.h()),h);
            k = std::max(std::abs(t.k()),k);
            l = std::max(std::abs(t.l()),l);
            for (int i=1 ; i != spgr.num_primops() ; ++i ) {
                clipper::HKL t1(t.transform(spgr.primitive_symop(i)) );
                h = std::max(std::abs(t1.h()),h);
                k = std::max(std::abs(t1.k()),k);
                l = std::max(std::abs(t1.l()),l);
            }
        } 
        return clipper::HKL(h,k,l);
    }
    
	void ReflectionData::output() const
	{
        clipper::HKL hkl(max_hkl() );
        clipper::HKL shkl(max_sym_hkl() );
        clipper::Cell cell(cdata_->base_cell() );
		printf("\nReflection Data INFO:\n\n");
		std::cout << "Reflection data type: " << cdata_->type() << std::endl;
		std::cout << "Number of observations (including Freidal mates): " << num_obs() << std::endl;
		std::cout << "Number of unique reflections (excluding Freidal): " << cdata_->num_obs() << " (Acentric: " << num_acentric() << ", Centric: " << num_centric() << ")" << std::endl;
		std::cout << "Resolution range of data: " << std::fixed << std::setprecision(3) << 1.0/std::sqrt(invresolsq_range().min() ) << " - " << 1.0/std::sqrt(invresolsq_range().max() ) << " A" << std::endl;
        std::cout << "Maximum index h (a*): " << std::fixed << std::setprecision(3) << hkl.h() << " (" << 1.0/std::sqrt( clipper::HKL(hkl.h(),0,0).invresolsq(cell) ) << " A)   - by symm. -  " << shkl.h() << " (" << 1.0/std::sqrt( clipper::HKL(shkl.h(),0,0).invresolsq(cell) ) << " A)" << std::endl;
        std::cout << "Maximum index k (b*): " << std::fixed << std::setprecision(3) << hkl.k() << " (" << 1.0/std::sqrt( clipper::HKL(0,hkl.k(),0).invresolsq(cell) ) << " A)   - by symm. -  " << shkl.k() << " (" << 1.0/std::sqrt( clipper::HKL(0,shkl.k(),0).invresolsq(cell) ) << " A)" << std::endl;
        std::cout << "Maximum index l (c*): " << std::fixed << std::setprecision(3) << hkl.l() << " (" << 1.0/std::sqrt( clipper::HKL(0,0,hkl.l()).invresolsq(cell) ) << " A)   - by symm. -  " << shkl.l() << " (" << 1.0/std::sqrt( clipper::HKL(0,0,shkl.l()).invresolsq(cell) ) << " A)" << std::endl;
        
		std::cout << std::endl;
	}
	
	std::stringstream& ReflectionData::xml_output(std::stringstream& ss) const
	/*<ReflectionData>
	 <NumberReflections>6701</NumberReflections>
	 <NumberObservations>9286</NumberObservations>
	 <NumberParts>9478</NumberParts>
	 <ResolutionHigh>    3.00</ResolutionHigh>
	 <NumberLattices>0</NumberLattices>
	 <NumberBatches>6</NumberBatches>
	 <NumberDatasets>1</NumberDatasets>
	 <Dataset  name="DMSO/DMSO/red_aucn">
	 <Run>
	 <number>1</number>
	 <Datasetname>DMSO/DMSO/red_aucn</Datasetname>
	 <BatchRange>5 10</BatchRange>
	 <BatchOffset>0</BatchOffset>
	 <PhiRange>343 349</PhiRange>
	 <FileStream>HKLIN</FileStream>
	 <Used>true</Used>
	 </Run>
	 </Dataset>
	 </ReflectionData>	*/
	{
        clipper::HKL hkl(max_hkl() );
        clipper::HKL shkl(max_sym_hkl() );
        clipper::Cell cell(cdata_->base_cell() );
		ss << "<ReflectionData>" << std::endl;
		ss << "<ReflectionDataType>" << cdata_->type() << "</ReflectionDataType>" << std::endl;
		ss << "<NumberObservations>" << num_obs() << "</NumberObservations>" << std::endl;
		ss << "<NumberReflections>" << cdata_->num_obs() << "</NumberReflections>" << std::endl;
		ss << "<NumberAcentric>" << num_acentric() << "</NumberAcentric>" << std::endl;
		ss << "<NumberCentric>" << num_centric() << "</NumberCentric>" << std::endl;
		ss << "<ResolutionLow>" << std::fixed << std::setprecision(3) << 1.0/std::sqrt(invresolsq_range().min() ) << "</ResolutionLow>" << std::endl;
		ss << "<ResolutionHigh>" << std::fixed << std::setprecision(3) << 1.0/std::sqrt(invresolsq_range().max() ) << "</ResolutionHigh>" << std::endl;
        ss << "<FileReflectionIndexMax id=\"h\">" << std::endl;
        ss << "<Index>" << hkl.h() << "</Index>" << std::endl;
        ss << "<Resolution unit=\"Angstrom\">" << std::fixed << std::setprecision(3) << 1.0/std::sqrt( clipper::HKL(hkl.h(),0,0).invresolsq(cell) ) << "</Resolution>" << std::endl;
        ss << "</FileReflectionIndexMax>" << std::endl;
        ss << "<FileReflectionIndexMax id=\"k\">" << std::endl;
        ss << "<Index>" << hkl.k() << "</Index>" << std::endl;
        ss << "<Resolution unit=\"Angstrom\">" << std::fixed << std::setprecision(3) << 1.0/std::sqrt( clipper::HKL(0,hkl.k(),0).invresolsq(cell) ) << "</Resolution>" << std::endl;
        ss << "</FileReflectionIndexMax>" << std::endl;
        ss << "<FileReflectionIndexMax id=\"l\">" << std::endl;
        ss << "<Index>" << hkl.l() << "</Index>" << std::endl;
        ss << "<Resolution unit=\"Angstrom\">" << std::fixed << std::setprecision(3) << 1.0/std::sqrt( clipper::HKL(0,0,hkl.l()).invresolsq(cell) ) << "</Resolution>" << std::endl;
        ss << "</FileReflectionIndexMax>" << std::endl;
        ss << "<ReflectionIndexMax id=\"h\">" << std::endl;
        ss << "<Index>" << shkl.h() << "</Index>" << std::endl;
        ss << "<Resolution unit=\"Angstrom\">" << std::fixed << std::setprecision(3) << 1.0/std::sqrt( clipper::HKL(shkl.h(),0,0).invresolsq(cell) ) << "</Resolution>" << std::endl;
        ss << "</ReflectionIndexMax>" << std::endl;
        ss << "<ReflectionIndexMax id=\"k\">" << std::endl;
        ss << "<Index>" << shkl.k() << "</Index>" << std::endl;
        ss << "<Resolution unit=\"Angstrom\">" << std::fixed << std::setprecision(3) << 1.0/std::sqrt( clipper::HKL(0,shkl.k(),0).invresolsq(cell) ) << "</Resolution>" << std::endl;
        ss << "</ReflectionIndexMax>" << std::endl;
        ss << "<ReflectionIndexMax id=\"l\">" << std::endl;
        ss << "<Index>" << shkl.l() << "</Index>" << std::endl;
        ss << "<Resolution unit=\"Angstrom\">" << std::fixed << std::setprecision(3) << 1.0/std::sqrt( clipper::HKL(0,0,shkl.l()).invresolsq(cell) ) << "</Resolution>" << std::endl;
        ss << "</ReflectionIndexMax>" << std::endl;

		ss << "</ReflectionData>" << std::endl;
		return ss;
	}
	
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
        // Copy rings
        void Rings::Copy(const Rings& other)
        {
                rings.clear();
                for (int i=0;i<other.nrings;i++)
                {
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
	/*template <class T> void Rings::AddObs(const int& Iring, const clipper::datatypes::I_sigI<T>& I_sigI,
					   const double& invresolsq, const double& multiplicity)
	{
		CheckRing(Iring);
		rings[Iring].AddObs(I_sigI, invresolsq, multiplicity );
	} */
	//--------------------------------------------------------------
	/*template <class T> void Rings::AddObs(const int& Iring, const clipper::datatypes::F_sigF<T>& I_sigI,
					   const double& invresolsq, const double& multiplicity)
	{
		CheckRing(Iring);
		rings[Iring].AddObs(I_sigI, invresolsq, multiplicity );
	}*/
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
	double Rings::N(const int& Iring) const
	{
		CheckRing(Iring);
		return rings[Iring].N();
	}
	//--------------------------------------------------------------
	double Rings::Comp(const int& Iring) const
	{
		CheckRing(Iring);
		return rings[Iring].Comp();
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
		nO = 0;
		reject = false;
	}
	//--------------------------------------------------------------
	/*template <class T> void IceRing::AddObs(const clipper::datatypes::I_sigI<T>& I_sigI, const double& invresolsq, const double& multiplicity)
	{
		if (!I_sigI.missing() ) {
			sum_I += multiplicity*I_sigI.I();
			sum_sigI += multiplicity*I_sigI.sigI();
			sum_sSqr += multiplicity*invresolsq;
			nI+=multiplicity;
		}
		nO += multiplicity;
	} */
	//--------------------------------------------------------------
	/*template <class T> void IceRing::AddObs(const clipper::datatypes::F_sigF<T>& I_sigI, const double& invresolsq, const double& multiplicity)
	{
		if (!I_sigI.missing() ) {
			sum_I += multiplicity*I_sigI.f();
			sum_sigI += multiplicity*I_sigI.sigf();
			sum_sSqr += multiplicity*invresolsq;
			nI+=multiplicity;
		}
		nO += multiplicity;
	}	 */
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
	double IceRing::Comp() const
	{
		if ( nO > 0)
		{return nI/nO;}
		else
		{return 0.0;}
	}
	//--------------------------------------------------------------
	
}

#include "ctruncate_pcf.cpp"

