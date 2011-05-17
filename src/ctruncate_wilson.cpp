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
#include "best.h"

namespace ctruncate {

	std::string name[5] = { "C", "N", "O", "H", "S" };
	
	std::vector<float> wilson_calc(clipper::HKL_data<clipper::data32::I_sigI>& isig, std::vector<int>& numatoms, float maxres, 
								   int nprm, CCP4Program& prog, clipper::HKL_data<clipper::data32::I_sigI>& ref)
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
				
		std::vector<double> params_ref( nprm, 1.0 );
		clipper::BasisFn_spline basis_ref( ref, nprm, 2.0 );
		TargetFn_meanInth<clipper::data32::I_sigI> target_ref( ref, 1 );
		clipper::ResolutionFn Sigma( hklinf, basis_ref, target_ref, params_ref );
		
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
				   log(basis_ref.f_s( res, Sigma.params() ))-log(totalscat),-0.5*a*res-b);
			else printf("%10.5f %10.5f %10.5f \n", res,log(basis_fo_wilson.f_s( res, wilsonplot.params() ))-log(totalscat),
						log(basis_ref.f_s( res, Sigma.params() ))-log(totalscat));
		}
		
		printf("$$\n\n");
		
		std::vector<float> params(2,0.0f); params[0]=b ; params[1] = -a;
		return params;
	}
	
	std::vector<float> wilson_plot(clipper::HKL_data<clipper::data32::I_sigI>& isig, float maxres, int nbins, CCP4Program& prog,
								   clipper::HKL_data<clipper::data32::I_sigI>& ref)
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
		
		wilson_calc(isig, numatoms, maxres, nbins, prog, ref);
	}
	
	std::vector<float> wilson_plot(clipper::HKL_data<clipper::data32::I_sigI>& isig, clipper::MPolymerSequence& poly, float maxres, 
								   int nbins, CCP4Program& prog, clipper::HKL_data<clipper::data32::I_sigI>& ref)
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
		
		wilson_calc(isig, numatoms, maxres, nbins, prog, ref);
		
	}
	
	std::vector<float> wilson_plot(clipper::HKL_data<clipper::data32::I_sigI>& isig, int nresidues, float maxres, int nbins, 
					 CCP4Program& prog, clipper::HKL_data<clipper::data32::I_sigI>& ref)
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
		
		
		wilson_calc(isig, numatoms, maxres, nbins, prog, ref);
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
	
	const std::string WilsonB::AtomNames[5] = { "C", "N", "O", "H", "S" };
	const char WilsonB::ResidueNames[21] = {'A','R','N','D','C','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'};
	const int WilsonB::Catoms[21] = { 3,  6,  4,  4,  3,  3,  5,  5,  2,  6,  6,  6,  6,  5,  9,  5,  3,  4, 11,  9,  5 };
	const int WilsonB::Hatoms[21] = { 7, 14,  8,  7,  7,  7, 10,  9,  5,  9, 13, 13, 14, 11, 11,  9,  7,  9, 12, 11, 11 };
	const int WilsonB::Natoms[21] = { 1,  4,  2,  1,  1,  1,  2,  1,  1,  3,  1,  1,  2,  1,  1,  1,  1,  1,  2,  1,  1 };
	const int WilsonB::Oatoms[21] = { 2,  2,  3,  4,  2,  2,  3,  4,  2,  2,  2,  2,  2,  2,  2,  2,  3,  3,  2,  3,  2 };
	const int WilsonB::Satoms[21] = { 0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0 };
	
	float WilsonB::operator()(clipper::HKL_data<clipper::data32::I_sigI>& isig, ctruncate::Rings *ice)
	{
		const clipper::HKL_info& hklinf = isig.hkl_info();
		int nsym = hklinf.spacegroup().num_symops();
		maxres = hklinf.resolution().invresolsq_limit();
		intensity = &isig;
		std::cout << "resolution: " << maxres << std::endl;
		
		nresidues = int(0.5*hklinf.cell().volume()/(nsym*157)); 
		numatoms[0] = 5*nresidues;
		numatoms[1] = int(1.35*nresidues);
		numatoms[2] = int(1.5*nresidues);
		numatoms[3] = 8*nresidues;
		numatoms[4] = int(0.05*nresidues);
		
		ctruncate::Rings icerings;
		icerings.DefaultIceRings();
		
		if (ice == NULL ) ice = &icerings;
			
		if (mode == BEST ) {
			wilson_best(isig, *ice);
		} else {
			wilson_straight(isig, *ice);
		}
		
		return this->B();
	}
	
	float WilsonB::operator()(clipper::HKL_data<clipper::data32::I_sigI>& isig, int nresidues, ctruncate::Rings *ice)
	{
		nresidue_supplied = true;
		const clipper::HKL_info& hklinf = isig.hkl_info();
		int nsym = hklinf.spacegroup().num_symops();
		maxres = hklinf.resolution().invresolsq_limit();
		intensity = &isig;
		
		numatoms[0] = 5*nresidues;
		numatoms[1] = int(1.35*nresidues);
		numatoms[2] = int(1.5*nresidues);
		numatoms[3] = 8*nresidues;
		numatoms[4] = int(0.05*nresidues);
		
		ctruncate::Rings icerings;
		icerings.DefaultIceRings();
		
		if (ice == NULL ) ice = &icerings;
			
			if (mode == BEST ) {
				wilson_best(isig, *ice);
			} else {
				wilson_straight(isig, *ice);
			}
		
		return this->B();
	}
	
	float WilsonB::operator()(clipper::HKL_data<clipper::data32::I_sigI>& isig, clipper::MPolymerSequence& poly, ctruncate::Rings *ice)
	{
		nresidue_supplied = true;
		const clipper::HKL_info& hklinf = isig.hkl_info();
		int nsym = hklinf.spacegroup().num_symops();
		maxres = hklinf.resolution().invresolsq_limit();
		intensity = &isig;
		
		clipper::String sequence = poly.sequence();
		for (int i=0; i != sequence.length(); ++i) {
			for (int j=0; j<21; j++) {
				if (sequence[i] == ResidueNames[j]) {
					numatoms[0] += Catoms[j];
					numatoms[1] += Natoms[j];
					numatoms[2] += Oatoms[j];
					numatoms[3] += Hatoms[j];
					numatoms[4] += Satoms[j];
					break;
				}
			}
		}
		
		ctruncate::Rings icerings;
		icerings.DefaultIceRings();
		
		if (ice == NULL ) ice = &icerings;
			
			if (mode == BEST ) {
				wilson_best(isig, *ice);
			} else {
				wilson_straight(isig, *ice);
			}
				
		return this->B();
	}
	
	void WilsonB::summary()
	{
		printf("\n\nWILSON SCALING:\n");
		if ( nresidue_supplied) {
			printf("User supplied number of residues = %d\n",nresidues);
		} else if ( sequence_supplied) {
			printf("User supplied sequence contains %d C, %d N, %d O, %d H, %d S atoms in %d residues\n", 
				   numatoms[0],numatoms[1],numatoms[2],numatoms[3],numatoms[4], nresidues);
			printf("User supplied number of residues = %d\n",nresidues);
		} else {
			printf("Estimated number of residues = %d\n",nresidues);
		}
		printf("\nResults Wilson plot:\n");
		if ( _a != 0.0 ) {
			if ( mode == BEST ) printf("Computed using Popov & Bourenkov, Acta D (2003) D59, 1145\n");
			printf ("B = %6.3f intercept = %6.3f siga = %6.3f sigb = %6.3f\n",this->B(),_b,_siga,_sigb);
			printf("scale factor on intensity = %10.4f\n\n",this->intercept());
		} else {
			printf("Too few high resolution points to determine B factor and Wilson scale factor\n");
		}
	}
	
	void WilsonB::plot(clipper::HKL_data<clipper::data32::I_sigI>& ref, clipper::String& name)
	{
		int nbins = 60;
		const clipper::HKL_info& hklinf = intensity->hkl_info();
		const clipper::HKL_info& hklinf_ref = ref.hkl_info();
		int nsym = hklinf.spacegroup().num_symops();
		
		float scale = 1.0f;
		float off = 0.0f;
		if (mode == WilsonB::STRAIGHT) {
			std::vector<float> xi(200), yi(200), wi(200,1.0), xxi(200), yyi(200), yy2i(200);
			for (int i = 0 ; i != 200 ; ++i ) {
				float res = 0.13223 + i*0.004091973;
				float totalscat = 0;
				for (int j=0;j!=5;++j) {
					clipper::Atom atom;
					atom.set_occupancy(1.0);
					atom.set_element(WilsonB::AtomNames[j]);
					atom.set_u_iso(0.0);
					atom.set_u_aniso_orth( clipper::U_aniso_orth( clipper::U_aniso_orth::null() ) ); // need this o/w next line hangs
					clipper::AtomShapeFn sf(atom);
					float scat = sf.f(res);
					totalscat +=  float( nsym * numatoms[j] ) * scat * scat;
				}
				yi[i] = totalscat;
				xi[i] = ctruncate::BEST(res);
			}
			int nobs = 200;
			float siga,sigb;
			straight_line_fit(xi,yi,wi,nobs,scale,off,siga,sigb);
		}
		
		//calc wilson plot (double handling unfortunately)
		std::vector<double> params_init( nbins, 1.0 );
		clipper::BasisFn_linear basis_fo_wilson( *intensity, nbins, 2.0 );
		TargetFn_meanInth<clipper::data32::I_sigI> target_fo_wilson( *intensity, 1);
		clipper::ResolutionFn wilsonplot( hklinf, basis_fo_wilson, target_fo_wilson, params_init );
		
		// Sigma or Normalisation curve
		// calculate Sigma (mean intensity in resolution shell) 
		// use intensities uncorrected for anisotropy
		std::vector<double> params_ref( nbins, 1.0 );
		clipper::BasisFn_spline basis_ref( ref, nbins, 2.0 );
		TargetFn_meanInth<clipper::data32::I_sigI> target_ref( ref, 1 );
		clipper::ResolutionFn Sigma( hklinf_ref, basis_ref, target_ref, params_ref );
		
		// end of Norm calc
		
		// wilson plot plus Norm curve
		printf("$TABLE: Wilson plot:\n");
		printf("$GRAPHS");
		//printf(": Wilson plot:0|0.1111x-7|-5:1,2:\n$$");  // limits hardwired
		if ( _a != -1.0 ) {
			printf(": Wilson plot - estimated B factor = %5.1f :A:1,2,3,4,5:\n$$", this->B());  
			printf(" 1/resol^2 ln(I/I_th) %s Linear Best $$\n$$\n", name.c_str());
		} else {
			printf(": Wilson plot :A:1,2,3:\n$$");  
			printf(" 1/resol^2 ln(I/I_th) %s $$\n$$\n", name.c_str());
		}
		
		for ( int i=0; i!=nbins; ++i ) {
			float res = maxres*(float(i)+0.5)/float(nbins); 
			float totalscat = 0;
			for (int j=0;j!=5;++j) {
				clipper::Atom atom;
				atom.set_occupancy(1.0);
				atom.set_element(WilsonB::AtomNames[j]);
				atom.set_u_iso(0.0);
				atom.set_u_aniso_orth( clipper::U_aniso_orth( clipper::U_aniso_orth::null() ) ); // need this o/w next line hangs
				clipper::AtomShapeFn sf(atom);
				float scat = sf.f(res);
				totalscat +=  float( nsym * numatoms[j] ) * scat * scat;
			}
			if ( _a != -1.0) 
				printf("%10.5f %10.5f %10.5f %10.5f %10.5f\n",
					   res,
					   log(basis_fo_wilson.f_s( res, wilsonplot.params() ))-log(totalscat),
					   log(basis_ref.f_s( res, Sigma.params() ))-log(totalscat),
					   -_a*res-_b, 
					   std::log(exp(-_a*res-_b)*(_bscale*scale*ctruncate::BEST(res)+_boff+off))-log(totalscat) );
			else printf("%10.5f %10.5f %10.5f \n", 
						res,
						log(basis_fo_wilson.f_s( res, wilsonplot.params() ))-log(totalscat),
						log(basis_ref.f_s( res, Sigma.params() ))-log(totalscat));
		}
		
		printf("$$\n\n");

	}
	
	void WilsonB::plot()
	{
		int nbins = 60;
		const clipper::HKL_info& hklinf = intensity->hkl_info();
		int nsym = hklinf.spacegroup().num_symops();
		
		float scale = 1.0f;
		float off = 0.0f;
		if (mode == WilsonB::STRAIGHT) {
			std::vector<float> xi(200), yi(200), wi(200,1.0), xxi(200), yyi(200), yy2i(200);
			for (int i = 0 ; i != 200 ; ++i ) {
				float res = 0.13223 + i*0.004091973;
				float totalscat = 0;
				for (int j=0;j!=5;++j) {
					clipper::Atom atom;
					atom.set_occupancy(1.0);
					atom.set_element(WilsonB::AtomNames[j]);
					atom.set_u_iso(0.0);
					atom.set_u_aniso_orth( clipper::U_aniso_orth( clipper::U_aniso_orth::null() ) ); // need this o/w next line hangs
					clipper::AtomShapeFn sf(atom);
					float scat = sf.f(res);
					totalscat +=  float( nsym * numatoms[j] ) * scat * scat;
				}
				yi[i] = totalscat;
				xi[i] = ctruncate::BEST(res);
			}
			int nobs = 200;
			float siga,sigb;
			straight_line_fit(xi,yi,wi,nobs,scale,off,siga,sigb);
		}
		
		//calc wilson plot (double handling unfortunately)
		std::vector<double> params_init( nbins, 1.0 );
		clipper::BasisFn_linear basis_fo_wilson( *intensity, nbins, 2.0 );
		TargetFn_meanInth<clipper::data32::I_sigI> target_fo_wilson( *intensity, 1);
		clipper::ResolutionFn wilsonplot( hklinf, basis_fo_wilson, target_fo_wilson, params_init );
		
		printf("$TABLE: Wilson plot:\n");
		printf("$GRAPHS");
		//printf(": Wilson plot:0|0.1111x-7|-5:1,2:\n$$");  // limits hardwired
		if ( _a != -1.0 ) {
			printf(": Wilson plot - estimated B factor = %5.1f :A:1,2,3,4:\n$$", this->B());  
			printf(" 1/resol^2 ln(I/I_th) Linear Best $$\n$$\n");
		} else {
			printf(": Wilson plot :A:1,2,3:\n$$");  
			printf(" 1/resol^2 ln(I/I_th) $$\n$$\n");
		}
					
		for ( int i=0; i!=nbins; ++i ) {
			float res = maxres*(float(i)+0.5)/float(nbins); 
			float totalscat = 0;
			for (int j=0;j!=5;++j) {
				clipper::Atom atom;
				atom.set_occupancy(1.0);
				atom.set_element(WilsonB::AtomNames[j]);
				atom.set_u_iso(0.0);
				atom.set_u_aniso_orth( clipper::U_aniso_orth( clipper::U_aniso_orth::null() ) ); // need this o/w next line hangs
				clipper::AtomShapeFn sf(atom);
				float scat = sf.f(res);
				totalscat +=  float( nsym * numatoms[j] ) * scat * scat;
			}
			if ( _a != -1.0) //printf("%10.5f %10.5f %10.5f %10.5f \n", res,log(basis_fo_wilson.f_s( res, wilsonplot.params() ))-log(totalscat),
							//		-_a*res-_b, -_a*res-_b+std::log(ctruncate::BEST(res))-std::log(totalscat));
					printf("%10.5f %10.5f %10.5f %10.5f \n", 
						   res,
						   basis_fo_wilson.f_s( res, wilsonplot.params() ),
						   exp(-_a*res-_b)*totalscat, 
						   exp(-_a*res-_b)*(scale*_bscale*ctruncate::BEST(res)+off+_boff));
			else printf("%10.5f %10.5f \n", 
						res,
						log(basis_fo_wilson.f_s( res, wilsonplot.params() ))-log(totalscat));
		}
		
		printf("$$\n\n");
		
	}
	
	void WilsonB::wilson_straight(clipper::HKL_data<clipper::data32::I_sigI>& isig, ctruncate::Rings& icerings)
	// common part
	{
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		// Wilson plot
		int nbins = 60;
		const clipper::HKL_info& hklinf = isig.hkl_info();
		int nsym = hklinf.spacegroup().num_symops();
		
		clipper::HKL_data<clipper::data32::I_sigI> xsig(hklinf);  // knock out ice rings and centric
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			double reso = ih.hkl().invresolsq(hklinf.cell());
			xsig[ih] = clipper::data32::I_sigI( (isig[ih].I()), isig[ih].sigI() );
			if ( ih.hkl_class().centric() ) xsig[ih].I() = xsig[ih].sigI() = clipper::Util::nan(); // loose centrics
			if ( icerings.InRing(reso) != -1 ) 
				if ( icerings.Reject( icerings.InRing(reso) ) ) xsig[ih].I() = xsig[ih].sigI() = clipper::Util::nan(); // loose ice rings
		}		
		
		std::vector<double> params_init( nbins, 1.0 );
		clipper::BasisFn_linear basis_fo_wilson( xsig, nbins, 2.0 );
		TargetFn_meanInth<clipper::data32::I_sigI> target_fo_wilson( xsig, 1);
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
					atom.set_element(WilsonB::AtomNames[i]);
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

		if ( wi.size() > 200 && maxres > maxres_scaling) {               // 3.5 Angstroms
			_a = _b = 0.0f;
			straight_line_fit(xi,yi,wi,nobs,_a,_b,_siga,_sigb);
		}
	}
	
	
	void WilsonB::wilson_best(clipper::HKL_data<clipper::data32::I_sigI>& isig, ctruncate::Rings& icerings)
	// common part
	{
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		// Wilson plot
		int nbins = 60;
		const clipper::HKL_info& hklinf = isig.hkl_info();
		int nsym = hklinf.spacegroup().num_symops();
		
		{
			std::vector<float> xi(200), yi(200), wi(200,1.0), xxi(200), yyi(200), yy2i(200);
			for (int i = 0 ; i != 200 ; ++i ) {
				float res = 0.13223 + i*0.004091973;
				float totalscat = 0;
				for (int j=0;j!=5;++j) {
					clipper::Atom atom;
					atom.set_occupancy(1.0);
					atom.set_element(WilsonB::AtomNames[j]);
					atom.set_u_iso(0.0);
					atom.set_u_aniso_orth( clipper::U_aniso_orth( clipper::U_aniso_orth::null() ) ); // need this o/w next line hangs
					clipper::AtomShapeFn sf(atom);
					float scat = sf.f(res);
					totalscat +=  float( nsym * numatoms[j] ) * scat * scat;
				}
				yi[i] = totalscat;
				xi[i] = ctruncate::BEST(res);
			}
			int nobs = 200;
			float siga,sigb;
			straight_line_fit(xi,yi,wi,nobs,_bscale,_boff,siga,sigb);
			std::cout << "scale and offset " << _bscale << " " << _boff << std::endl;
			
		}
		
		clipper::HKL_data<clipper::data32::I_sigI> xsig(hklinf);  // knock out ice rings and centric
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			double reso = ih.hkl().invresolsq(hklinf.cell());
			xsig[ih] = clipper::data32::I_sigI( (isig[ih].I()), isig[ih].sigI() );
			if ( ih.hkl_class().centric() ) xsig[ih].I() = xsig[ih].sigI() = clipper::Util::nan(); // loose centrics
			if ( icerings.InRing(reso) != -1 ) 
				if ( icerings.Reject( icerings.InRing(reso) ) ) xsig[ih].I() = xsig[ih].sigI() = clipper::Util::nan(); // loose ice rings
		}		
		
		std::vector<double> params_init( nbins, 1.0 );
		clipper::BasisFn_linear basis_fo_wilson( xsig, nbins, 2.0 );
		TargetFn_meanInth<clipper::data32::I_sigI> target_fo_wilson( xsig, 1);
		clipper::ResolutionFn wilsonplot( hklinf, basis_fo_wilson, target_fo_wilson, params_init );
		
		std::vector<float> xi, yi, wi, xxi, yyi, yy2i; 
		float totalscat; 
		
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			if ( !isig[ih].missing() && wilsonplot.f(ih) > 0.0) {
				float lnS = -log(wilsonplot.f(ih));
				float res = ih.invresolsq();
				
				totalscat = _bscale*ctruncate::BEST(res)+_boff;
				lnS += log(totalscat);
				
				if ( icerings.InRing(ih.hkl().invresolsq(isig.base_cell()  ) == -1 ) ) {  
					xi.push_back(res);
					yi.push_back(lnS);
					//float weight = pow(isig[ih].sigI(),2);
					float weight = ( res < 0.04 ) ? std::pow(0.04/res, 2.0)*isig[ih].sigI() : isig[ih].sigI();
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
		
		if ( wi.size() > 200 ) {
			_a = _b = 0.0f;
			straight_line_fit(xi,yi,wi,nobs,_a,_b,_siga,_sigb);
		}
	}
	
			
}