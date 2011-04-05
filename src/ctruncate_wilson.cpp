//
//     CTRUNCATE
//     Copyright (C) 2006-2011 Norman Stein, Charles Ballard
//
//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Application.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//

#include "ctruncate_wilson.h"
#include "ctruncate_utils.h"
#include "intensity_scale.h"
#include "intensity_target.h"

namespace ctruncate {

	std::string name[5] = { "C", "N", "O", "H", "S" };
	
	std::vector<float> wilson_calc(clipper::HKL_data<clipper::data32::I_sigI>& isig, std::vector<int>& numatoms, float maxres, 
								   int nprm, CCP4Program& prog)
	// common part
	{
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		// Wilson plot
		
		const clipper::HKL_info& hklinf = isig.hkl_info();
		int nsym = hklinf.spacegroup().num_symops();
		
		ctruncate::Rings icerings;
		icerings.DefaultIceRings();
		
		std::vector<double> params_init( nprm, 1.0 );
		clipper::BasisFn_linear basis_fo_wilson( isig, nprm, 2.0 );
		TargetFn_meanInth<clipper::data32::I_sigI> target_fo_wilson( isig, 1);
		clipper::ResolutionFn wilsonplot( hklinf, basis_fo_wilson, target_fo_wilson, params_init );
		
		std::vector<float> xi, yi, wi, xxi, yyi, yy2i; 
		float totalscat; 
		const float minres_scaling = 0.0625;   // 4 Angstroms
		const float maxres_scaling = 0.0816;    // 3.5 Angstroms
		
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			if ( !isig[ih].missing() && wilsonplot.f(ih) > 0.0) {
				float lnS = -log(wilsonplot.f(ih));
				float res = ih.invresolsq();
				
				totalscat = 0;
				for (int i=0;i!=5;++i) {
					clipper::Atom atom;
					atom.set_occupancy(1.0);
					atom.set_element(name[i]);
					atom.set_u_iso(0.0);
					atom.set_u_aniso_orth( clipper::U_aniso_orth( clipper::U_aniso_orth::null() ) ); // need this o/w next line hangs
					clipper::AtomShapeFn sf(atom);
					float scat = sf.f(res);
					totalscat +=  float( nsym * numatoms[i] ) * scat * scat;
				}
				lnS += log(totalscat);
				
				if (res > minres_scaling && ( icerings.InRing(ih.hkl().invresolsq(isig.base_cell() ) ) == -1 ) ) {  
					xi.push_back(res);
					yi.push_back(lnS);
					//float weight = pow(isig[ih].sigI(),2);
					float weight = isig[ih].sigI();
					//if (res > 0.1) printf("%f\n",weight);
					if (weight > 0.0) {
						wi.push_back(1.0/weight);
					}
					else {
						wi.push_back(0.0);
					}
				}
			}
		}
		
		int nobs = xi.size();
		//printf("%d %d %d\n", xi.size(), yi.size(), wi.size());
		float a,b,siga,sigb,a1,b1;
		b = 0.0; a = 0.0f;
		bool line(false);
		if ( wi.size() > 200 && maxres > maxres_scaling) {               // 3.5 Angstroms
			straight_line_fit(xi,yi,wi,nobs,a,b,siga,sigb);
			prog.summary_beg();
			printf("\nResults Wilson plot:\n");
			a *= 2.0;
			printf ("B = %6.3f intercept = %6.3f siga = %6.3f sigb = %6.3f\n",a,b,siga,sigb);
			printf("scale factor on intensity = %10.4f\n\n",(exp(b)));
			prog.summary_end();
			line = true;
		} else {
			printf("Too few high resolution points to determine B factor and Wilson scale factor\n");
		}
		
		// Sigma or Normalisation curve
		// calculate Sigma (mean intensity in resolution shell) 
		// use intensities uncorrected for anisotropy
		
		int nprm2 = 12;
		
		clipper::HKL_data<clipper::data32::I_sigI> xsig(hklinf);  // knock out ice rings and centric
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			double reso = ih.hkl().invresolsq(hklinf.cell());
			xsig[ih] = clipper::data32::I_sigI( isig[ih].I(), isig[ih].sigI() );
			if ( ih.hkl_class().centric() ) xsig[ih].I() = xsig[ih].sigI() = clipper::Util::nan(); // loose centrics
			if ( icerings.InRing(reso) != -1 )
				xsig[ih].I() = xsig[ih].sigI() = clipper::Util::nan(); // loose ice rings
		}		
		std::vector<double> params_ice( nprm2, 1.0 );
		clipper::BasisFn_spline basis_fo( xsig, nprm2, 2.0 );
		TargetFn_meanInth<clipper::data32::I_sigI> target_fo( xsig, 1 );
		clipper::ResolutionFn Sigma( hklinf, basis_fo, target_fo, params_ice );
		
		// end of Norm calc
		
		
		// wilson plot plus Norm curve
		printf("$TABLE: Wilson plot:\n");
		printf("$GRAPHS");
		//printf(": Wilson plot:0|0.1111x-7|-5:1,2:\n$$");  // limits hardwired
		if (line) {
			printf(": Wilson plot - estimated B factor = %5.1f :A:1,2,3,4:\n$$", a);  
			printf(" 1/resol^2 ln(I/I_th) Sigma Overall-B $$\n$$\n");
		} else {
			printf(": Wilson plot :A:1,2,3:\n$$");  
			printf(" 1/resol^2 ln(I/I_th) Sigma $$\n$$\n");
		}
		
		int nbins = 60;
		for ( int i=0; i!=nbins; ++i ) {
			float res = maxres*(float(i)+0.5)/float(nbins); 
			float totalscat = 0;
			for (int i=0;i!=5;++i) {
				clipper::Atom atom;
				atom.set_occupancy(1.0);
				atom.set_element(name[i]);
				atom.set_u_iso(0.0);
				atom.set_u_aniso_orth( clipper::U_aniso_orth( clipper::U_aniso_orth::null() ) ); // need this o/w next line hangs
				clipper::AtomShapeFn sf(atom);
				float scat = sf.f(res);
				totalscat +=  float( nsym * numatoms[i] ) * scat * scat;
			}
			if (line) printf("%10.5f %10.5f %10.5f %10.5f \n", res,log(basis_fo_wilson.f_s( res, wilsonplot.params() ))-log(totalscat),
				   log(basis_fo.f_s( res, Sigma.params() ))-log(totalscat),-0.5*a*res-b);
			else printf("%10.5f %10.5f %10.5f \n", res,log(basis_fo_wilson.f_s( res, wilsonplot.params() ))-log(totalscat),
						log(basis_fo.f_s( res, Sigma.params() ))-log(totalscat));
		}
		
		printf("$$\n\n");
		
		std::vector<float> params(2,0.0f); params[0]=b ; params[1] = -a;
		return params;
	}
	
	std::vector<float> wilson_plot(clipper::HKL_data<clipper::data32::I_sigI>& isig, float maxres, int nbins, CCP4Program& prog)
	{
		const clipper::HKL_info& hklinf = isig.hkl_info();
		int nsym = hklinf.spacegroup().num_symops();
		
		int nresidues = int(0.5*hklinf.cell().volume()/(nsym*157));
		prog.summary_beg();
		printf("\n\nWILSON SCALING:\n");
		printf("Estimated number of residues = %d\n",nresidues);
		prog.summary_end();
		
		std::vector<int> numatoms(5,0); 
		numatoms[0] = 5*nresidues;
		numatoms[1] = int(1.35*nresidues);
		numatoms[2] = int(1.5*nresidues);
		numatoms[3] = 8*nresidues;
		numatoms[4] = int(0.05*nresidues);
		
		wilson_calc(isig, numatoms, maxres, nbins, prog);
	}
	
	std::vector<float> wilson_plot(clipper::HKL_data<clipper::data32::I_sigI>& isig, clipper::MPolymerSequence& poly, float maxres, 
					 int nbins, CCP4Program& prog)
	{
		std::vector<int> numatoms(5,0); 
		
		// Use single letter residue names from mmdb - note that C appears twice
		
		//                   A   R   N   D   C   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
		char ResidueName1[21] = {'A','R','N','D','C','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'};
		int Catoms[21] = { 3,  6,  4,  4,  3,  3,  5,  5,  2,  6,  6,  6,  6,  5,  9,  5,  3,  4, 11,  9,  5 };
		int Hatoms[21] = { 7, 14,  8,  7,  7,  7, 10,  9,  5,  9, 13, 13, 14, 11, 11,  9,  7,  9, 12, 11, 11 };
		int Natoms[21] = { 1,  4,  2,  1,  1,  1,  2,  1,  1,  3,  1,  1,  2,  1,  1,  1,  1,  1,  2,  1,  1 };
		int Oatoms[21] = { 2,  2,  3,  4,  2,  2,  3,  4,  2,  2,  2,  2,  2,  2,  2,  2,  3,  3,  2,  3,  2 };
		int Satoms[21] = { 0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0 };
		
		
		clipper::String sequence = poly.sequence();
		for (int i=0; i != sequence.length(); ++i) {
			for (int j=0; j<21; j++) {
				if (sequence[i] == ResidueName1[j]) {
					numatoms[0] += Catoms[j];
					numatoms[1] += Natoms[j];
					numatoms[2] += Oatoms[j];
					numatoms[3] += Hatoms[j];
					numatoms[4] += Satoms[j];
					break;
				}
			}
		}
		prog.summary_beg();
		printf("\n\nWILSON SCALING:\n");
		printf("User supplied sequence contains %d C, %d N, %d O, %d H, %d S atoms\n", 
			   numatoms[0],numatoms[1],numatoms[2],numatoms[3],numatoms[4]);
		prog.summary_end();
		
		wilson_calc(isig, numatoms, maxres, nbins, prog);
		
	}
	
	std::vector<float> wilson_plot(clipper::HKL_data<clipper::data32::I_sigI>& isig, int nresidues, float maxres, int nbins, 
					 CCP4Program& prog)
	{
		prog.summary_beg();
		printf("\n\nWILSON SCALING:\n");
		printf("User supplied number of residues = %d\n",nresidues);
		prog.summary_end();
		
		std::vector<int> numatoms(5,0); 
		numatoms[0] = 5*nresidues;
		numatoms[1] = int(1.35*nresidues);
		numatoms[2] = int(1.5*nresidues);
		numatoms[3] = 8*nresidues;
		numatoms[4] = int(0.05*nresidues);
		
		
		wilson_calc(isig, numatoms, maxres, nbins, prog);
	}	
	
	// Truncate style Wilson plot
	/*(
	 std::vector<int> N_all(nbins,0);
	 std::vector<int> N_obs(nbins,0);
	 std::vector<float> I_obs(nbins,0.0);
	 
	 std::vector<float> xtr, ytr, wtr; 
	 
	 xxi.clear();
	 yyi.clear();
	 
	 
	 for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
	 int bin = int( float(nbins) * ih.invresolsq() / maxres - 0.5 );
	 if (bin >= nbins || bin < 0) printf("Warning: (Wilson 2) illegal bin number %d\n", bin);
	 N_all[bin]++;
	 if ( !isig[ih].missing() ) {
	 I_obs[bin] += isig[ih].I();
	 N_obs[bin]++;
	 }
	 }
	 
	 for ( int j=0; j!=nbins; ++j ) {
	 float res = maxres*(float(j)+0.5)/float(nbins);
	 totalscat = 0;
	 for (int i=0;i!=5;++i) {
	 Atom atom;
	 atom.set_occupancy(1.0);
	 atom.set_element(name[i]);
	 atom.set_u_iso(0.0);
	 atom.set_u_aniso_orth( U_aniso_orth( U_aniso_orth::null() ) ); // need this o/w next line hangs
	 AtomShapeFn sf(atom);
	 float scat = sf.f(res);
	 totalscat +=  float( nsym * numatoms[i] ) * scat * scat;
	 }
	 
	 if (I_obs[j] > 0.0) {
	 float w1 = log( I_obs[j] / (float(N_obs[j]) * totalscat) );
	 float w2 = log( I_obs[j] / (float(N_all[j]) * totalscat) );
	 
	 xxi.push_back(res);
	 yyi.push_back(w1);
	 yy2i.push_back(w2);
	 
	 if (res > minres_scaling) {  
	 xtr.push_back(res);
	 ytr.push_back(w1);
	 wtr.push_back(1.0);
	 }
	 }
	 }
	 
	 nobs = xtr.size();
	 if ( wi.size() > 200 && maxres > 0.0816) {
	 straight_line_fit(xtr,ytr,wtr,nobs,a1,b1,siga,sigb);
	 prog.summary_beg();
	 printf("\nresults from fitting Truncate style Wilson plot\n");
	 printf ("B = %6.3f intercept = %6.3f siga = %6.3f sigb = %6.3f\n",-2.0*a1,-b1,siga,sigb);
	 printf("scale factor on intensity = %10.4f\n\n", exp(-b1));
	 prog.summary_end();
	 }
	 
	 printf("$TABLE: Truncate style Wilson plot:\n");
	 printf("$GRAPHS");
	 printf(": Wilson plot - estimated B factor = %5.1f :A:1,2,3:\n$$", -2.0*a1);  
	 //printf(": Wilson plot:0|0.1111x-8|-5.5:1,2,3:\n$$");  // limits hardwired
	 printf(" 1/resol^2 obs all $$\n$$\n");
	 for ( int i=0; i<xxi.size(); i++ ) {
	 printf( "%10.6f %10.6f %10.6f\n", xxi[i], yyi[i], yy2i[i]);
	 }
	 printf("$$\n\n");
*/	 

}