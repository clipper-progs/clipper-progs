//
//     CTRUNCATE
//     Copyright (C) 2006-2011 Norman Stein, Charles Ballard
//
//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Application.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//

#include "ctruncate_moments.h"

namespace ctruncate {
	
	void moments_Z(clipper::HKL_data<clipper::data32::I_sigI>& isig, float maxres, int nbins)
	{
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		
		int ncentric = 0;
		
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			if ( !isig[ih].missing() )
				if (ih.hkl_class().centric()) ncentric += 1;
		}
		
		int ncbins = std::min(nbins, ncentric/10);
		
		std::vector<int> Na(nbins,0),  Nc(ncbins,0);
		std::vector<double> E1a(nbins,0.0), E1c(ncbins,0.0);
		std::vector<double> E3a(nbins,0.0), E3c(ncbins,0.0);
		std::vector<double> I1a(nbins,0.0), I1c(ncbins,0.0);
		std::vector<double> I2a(nbins,0.0), I2c(ncbins,0.0);
		std::vector<double> I3a(nbins,0.0), I3c(ncbins,0.0);
		std::vector<double> I4a(nbins,0.0), I4c(ncbins,0.0);
		
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			if ( !isig[ih].missing() ) {
				//if( ih.hkl_class().epsilon() > 1.01) printf("%d %d %d epsilon = %f\n",ih.hkl().h(), ih.hkl().k(), ih.hkl().l(), ih.hkl_class().epsilon());
				double I = isig[ih].I() / ih.hkl_class().epsilon();
				int bin = int( nbins * pow( ih.invresolsq() / double(maxres), 1.0 ) - 0.5  );
				int cbin = int( ncbins * pow( ih.invresolsq() / double(maxres), 1.0 ) - 0.5  );
				if (bin >= nbins || bin < 0) printf("Warning: (moments) illegal bin number %d\n", bin);
				//printf("%3d %11.4f %11.6f\n",bin,I,ih.invresolsq());
				if (!ih.hkl_class().centric()) {
					Na[bin]++;
					if (I > 0.0) {
						E1a[bin] += sqrt(I);
						I1a[bin] += I;
						E3a[bin] += pow(I,1.5);
						I2a[bin] += I*I;
						I3a[bin] += I*I*I;
						I4a[bin] += I*I*I*I;
					}
				}
				else if (ncentric != 0) {
					Nc[cbin]++;
					if (cbin >= ncbins || cbin < 0) printf("Warning: (moments) illegal cbin number %d\n", cbin);
					if (I > 0.0) {
						E1c[cbin] += sqrt(I);
						I1c[cbin] += I;
						E3c[cbin] += pow(I,1.5);
						I2c[cbin] += I*I;
						I3c[cbin] += I*I*I;
						I4c[cbin] += I*I*I*I;
					}
				}
			}
		}
		
		printf("$TABLE: Acentric moments of E:\n");
		printf("$GRAPHS");
		printf(": 4th moment of E (Expected value = 2, Perfect Twin = 1.5):0|%5.3fx0|5:1,4:\n", maxres);
		printf(": 1st & 3rd moments of E (Expected values = 0.886, 1.329, Perfect twin = 0.94, 1.175):0|%5.3fx0|2:1,2,3:\n", maxres);
		//printf(": 6th & 8th moments of E (Expected value = 6, 24, Perfect Twin 3, 7.5):0|%5.3fx0|48:1,5,6:\n", maxres);
		printf("$$ 1/resol^2 <E> <E**3> <E**4> <E**6> <E**8> $$\n$$\n");
		
		for (int i=0; i<nbins; i++) {
			double res = maxres * pow((double(i) + 0.5)/double(nbins), 1.00000);
			double n = double(Na[i]);
			if (n > 0 && I1a[i] > 0.0) {
				E1a[i] /= sqrt(I1a[i]*n);
				E3a[i] *= sqrt(n) / pow(I1a[i],1.5);
				I2a[i] *= n /( I1a[i]*I1a[i] );
				I3a[i] *= n*n / pow(I1a[i],3);
				I4a[i] *= pow(n,3) / pow(I1a[i],4);
			}
			printf("%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", res, E1a[i], E3a[i], I2a[i], I3a[i], I4a[i]);
		}
		printf("$$\n\n");
		
		if (ncentric != 0) {
			printf("$TABLE: Centric moments of E:\n");
			printf("$GRAPHS");
			printf(": 4th moment of E (Expected = 3, Perfect Twin = 2):0|%5.3fx0|5:1,4:\n", maxres);
			printf(": 1st & 3rd moments of E (Expected = 0.798, 1.596, Perfect Twin = 0.886, 1.329):0|%5.3fx0|4:1,2,3:\n", maxres);
			//printf(": 6th & 8th moments of E (Expected = 15, 105, Perfect Twin = 6, 24):0|%5.3fx0|120:1,5,6:\n", maxres);
			printf("$$ 1/resol^2 <E> <E**3> <E**4> <E**6> <E**8> $$\n$$\n");
			
			for (int i=0; i<ncbins; i++) {
				double res = maxres * pow((double(i) + 0.5)/double(ncbins), 0.666666);
				double n = double(Nc[i]);
				if (n > 0 && I1c[i] > 0.0) {
					E1c[i] /= sqrt(I1c[i]*n);
					E3c[i] *= sqrt(n) / pow(I1c[i],1.5);
					I2c[i] *= n /( I1c[i]*I1c[i] );
					I3c[i] *= n*n / pow(I1c[i],3);
					I4c[i] *= pow(n,3) / pow(I1c[i],4);
				}
				printf("%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", res, E1c[i], E3c[i], I2c[i], I3c[i], I4c[i]);
			}
			printf("$$\n\n");
		}
		return;
	}
	
	void moments_E(clipper::HKL_data<clipper::data32::I_sigI>& isig, float maxres, int nbins, CCP4Program& prog)
	{
		// moments of E using clipper binning
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		const clipper::HKL_info& hklinf = isig.hkl_info();
		
		//acentrics
		clipper::HKL_data<clipper::data32::F_sigF> fs1( hklinf );
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			if ( !isig[ih].missing() && !ih.hkl_class().centric() ) {
				double I = isig[ih].I();
				double sigI = isig[ih].sigI();
				if ( I > 0.0 )
					fs1[ih] = clipper::data32::F_sigF( sqrt(I), 0.5*sigI/sqrt(I) );
			}
		}
		
		//printf("%d observations with I > 0\n\n", fs1.num_obs());
		
		int nprmk = nbins;
		std::vector<double> params_initk( nprmk, 1.0 );
		clipper::BasisFn_binner basis_fn1( fs1, nprmk, 1.5 );  // equal increments in invresolsq bins
		//clipper::BasisFn_binner basis_fn1( fs1, nprmk, 1.0 );  // equal volume bins
		
		
		clipper::TargetFn_meanFnth<clipper::data32::F_sigF> target_fn1( fs1, 1.0 );
		clipper::TargetFn_meanFnth<clipper::data32::F_sigF> target_fn2( fs1, 2.0 );
		clipper::TargetFn_meanFnth<clipper::data32::F_sigF> target_fn3( fs1, 3.0 );
		clipper::TargetFn_meanFnth<clipper::data32::F_sigF> target_fn4( fs1, 4.0 );
		
		clipper::ResolutionFn f1( hklinf, basis_fn1, target_fn1, params_initk );
		clipper::ResolutionFn f2( hklinf, basis_fn1, target_fn2, params_initk );
		clipper::ResolutionFn f3( hklinf, basis_fn1, target_fn3, params_initk );
		clipper::ResolutionFn f4( hklinf, basis_fn1, target_fn4, params_initk );
		
		//printf("$TABLE: Acentric moments of E for k=1,3,4,6,8:\n");
		printf("$TABLE: Acentric moments of E for k=1,3,4:\n");
		printf("$GRAPHS");
		printf(": 4th moment of E (Expected value = 2, Perfect Twin = 1.5):0|%5.3fx0|5:1,4:\n", maxres);
		printf(": 1st & 3rd moments of E (Expected values = 0.886, 1.329, Perfect twin = 0.94, 1.175):0|%5.3fx0|2:1,2,3:\n", maxres);
		//printf(": 6th & 8th moments of E (Expected value = 6, 24, Perfect Twin 3, 7.5):0|%5.3fx0|48:1,5,6:\n", maxres);
		
		//printf("$$ 1/resol^2 <E> <E**3> <E**4> <E**6> <E**8> $$\n$$\n");
		printf("$$ 1/resol^2 <E> <E**3> <E**4> $$\n$$\n");
		
		double mean1, mean3, mean4;
		mean1 = mean3 = mean4 = 0.0;
		for (int i=0; i != nbins; i++) {
			double res = double(i+1) * maxres / double(nbins);   // equal increments in invresolsq bins
			//double res = maxres * pow( double(i+1)/double(nbins), 0.666666 );  // equal volume bins
			double i1 = basis_fn1.f_s( res, f2.params() );
			printf("%10.6f %10.6f %10.6f %10.6f \n", res,
				   basis_fn1.f_s( res, f1.params() )/pow(i1,0.5),
				   basis_fn1.f_s( res, f3.params() )/pow(i1,1.5), 
				   basis_fn1.f_s( res, f4.params() )/pow(i1,2.0) ); 
			
			mean1 += basis_fn1.f_s( res, f1.params() )/pow(i1,0.5);
			mean3 += basis_fn1.f_s( res, f3.params() )/pow(i1,1.5); 
			mean4 += basis_fn1.f_s( res, f4.params() )/pow(i1,2.0);
		}
		printf("$$\n\n");
		
		mean1 /= double(nbins);
		mean3 /= double(nbins);
		mean4 /= double(nbins);
		
		prog.summary_beg();
		printf("\n\nMEAN ACENTRIC MOMENTS OF E:\n\n");
		//printf("mean1 = %6.3f %6.3f %6.3f %6.2f %6.2f\n",mean1,mean3,mean4,mean6,mean8);
		printf("<E> = %6.3f (Expected value = 0.886, Perfect Twin = 0.94)\n", mean1);
		printf("<E**3> = %6.3f (Expected value = 1.329, Perfect Twin = 1.175)\n", mean3);
		printf("<E**4> = %6.3f (Expected value = 2, Perfect Twin = 1.5)\n", mean4);
		if (mean4 < 2.0 && mean4 > 1.5) {
			double alpha = 0.5 - sqrt(0.5*mean4 - 0.75);
			printf("(equivalent to twin fraction of %6.3f)\n",alpha);
		}
		//printf("<E**6> = %6.2f (Expected value = 6, Perfect Twin = 3)\n", mean6);
		//printf("<E**8> = %6.2f (Expected value = 24, Perfect Twin = 7.5)\n", mean8);
		prog.summary_end();
		
		int ncentric = 0;
		
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			if ( !isig[ih].missing() )
				if (ih.hkl_class().centric()) ncentric += 1;
		}
		
		int ncbins = std::min(nbins, ncentric/10);
		
		nprmk = ncbins;
		//centrics
		if (ncentric) {
			clipper::HKL_data<clipper::data32::F_sigF> fs2( hklinf );
			for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
				if ( !isig[ih].missing() && ih.hkl_class().centric() ) {
					double I = isig[ih].I();
					double sigI = isig[ih].sigI();
					if ( I > 0.0 )
						fs2[ih] = clipper::data32::F_sigF( sqrt(I), 0.5*sigI/sqrt(I) );
				}
			}
			
			std::vector<double> params_initk2( nprmk, 1.0 );
			clipper::BasisFn_binner basis_fn2( fs2, nprmk, 1.5 );   // equal increments in invresolsq bins
			//clipper::BasisFn_binner basis_fn2( fs2, nprmk, 1.0 );   // equal volume bins
			
			clipper::TargetFn_meanFnth<clipper::data32::F_sigF> target_fn1c( fs2, 1.0 );
			clipper::TargetFn_meanFnth<clipper::data32::F_sigF> target_fn2c( fs2, 2.0 );
			clipper::TargetFn_meanFnth<clipper::data32::F_sigF> target_fn3c( fs2, 3.0 );
			clipper::TargetFn_meanFnth<clipper::data32::F_sigF> target_fn4c( fs2, 4.0 );
			
			clipper::ResolutionFn f1c( hklinf, basis_fn2, target_fn1c, params_initk2 );
			clipper::ResolutionFn f2c( hklinf, basis_fn2, target_fn2c, params_initk2 );
			clipper::ResolutionFn f3c( hklinf, basis_fn2, target_fn3c, params_initk2 );
			clipper::ResolutionFn f4c( hklinf, basis_fn2, target_fn4c, params_initk2 );
			
			//printf("$TABLE: Centric moments of E for k=1,3,4,6,8:\n");
			printf("$TABLE: Centric moments of E for k=1,3,4:\n");
			printf("$GRAPHS");
			printf(": 4th moment of E (Expected = 3, Perfect Twin = 2):0|%5.3fx0|5:1,4:\n", maxres);
			printf(": 1st & 3rd moments of E (Expected = 0.798, 1.596, Perfect Twin = 0.886, 1.329):0|%5.3fx0|4:1,2,3:\n", maxres);
			//printf(": 6th & 8th moments of E (Expected = 15, 105, Perfect Twin = 6, 24):0|%5.3fx0|120:1,5,6:\n", maxres);
			
			//printf("$$ 1/resol^2 <E> <E**3> <E**4> <E**6> <E**8> $$\n$$\n");
			printf("$$ 1/resol^2 <E> <E**3> <E**4> $$\n$$\n");
			
			
			for (int i=0; i != ncbins; i++) {
				double res = double(i+1) * maxres / double(ncbins);   // equal increments in invresolsq bins
				//double res = maxres * pow( double(i+1)/double(ncbins), 0.666666 );  // equal volume bins
				double i1 = basis_fn2.f_s( res, f2c.params() );
				printf("%10.6f %10.6f %10.6f %10.6f\n", res,
					   basis_fn2.f_s( res, f1c.params() )/pow(i1,0.5),
					   basis_fn2.f_s( res, f3c.params() )/pow(i1,1.5),
					   basis_fn2.f_s( res, f4c.params() )/pow(i1,2.0) ); 
				//basis_fn2.f_s( res, f6c.params() )/pow(i1,3.0), 
				//basis_fn2.f_s( res, f8c.params() )/pow(i1,4.0) ); 
			}
			
			printf("$$\n\n");
		}
	}		

}