//
//     CTRUNCATE
//     Copyright (C) 2006-2011 Norman Stein, Charles Ballard
//
//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Application.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//

#include "ctruncate_parity.h"

namespace ctruncate {
	
	int parity(clipper::HKL_data<clipper::data32::I_sigI>& isig, float maxres, int nbins ) 
	{
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		//printf("Starting parity group analysis:\n");
		
		//Parity group analysis
		
		printf( "Analysis of mean intensity by parity for reflection classes\n\n");
		printf("For each class, Mn(I/sig(I)) is given for even and odd parity with respect to the condition,\n");
		printf("eg group 1: h even & odd; group 7 h+k+l even & odd; group 8 h+k=2n & h+l=2n & k+l=2n or not\n\n");
		printf( " Range    Min_S    Dmax    Nref     1           2           3           4           5           6           7           8\n");
		printf( "                                    h           k           l          h+k         h+l         k+l        h+k+l    h+k,h+l,k+l\n");
		
		float Iparity[8][2][60], Itot[8][2];
		int Nparity[8][2][60], Ntot[8][2];
		
		for ( int i1=0; i1<8; i1++) { 
			for ( int i2=0; i2<2; i2++) {
				Itot[i1][i2] = 0.0;
				Ntot[i1][i2] = 0;
				for ( int i3=0; i3<nbins; i3++) {	
					Iparity[i1][i2][i3] = 0.0;
					Nparity[i1][i2][i3] = 0;
				}
			}
		}
		
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			if ( !isig[ih].missing() ) {
				clipper::HKL hkl = ih.hkl();
				int bin = int( double(nbins) * ih.invresolsq()/ double(maxres) - 0.5  );
				if (bin >= nbins || bin < 0) printf("Warning: (parity) illegal bin number %d\n", bin);
				//if ( ih.hkl_class().centric() ) printf("centric: %d %d %d\n", hkl.h(), hkl.k(), hkl.l() );
				int h = hkl.h();
				int k = hkl.k();
				int l = hkl.l();
				float I_over_sigma = 0.0;
				if ( isig[ih].sigI() > 0.0 ) I_over_sigma = isig[ih].I() / isig[ih].sigI();
				
				Iparity[0][abs(h%2)][bin] += I_over_sigma;
				Iparity[1][abs(k%2)][bin] += I_over_sigma;
				Iparity[2][abs(l%2)][bin] += I_over_sigma;
				Iparity[3][abs((h+k)%2)][bin] += I_over_sigma;
				Iparity[4][abs((h+l)%2)][bin] += I_over_sigma;
				Iparity[5][abs((k+l)%2)][bin] += I_over_sigma;
				Iparity[6][abs((h+k+l)%2)][bin] += I_over_sigma;
				
				Nparity[0][abs(h%2)][bin] ++;
				Nparity[1][abs(k%2)][bin] ++;
				Nparity[2][abs(l%2)][bin] ++;
				Nparity[3][abs((h+k)%2)][bin] ++;
				Nparity[4][abs((h+l)%2)][bin] ++;
				Nparity[5][abs((k+l)%2)][bin] ++;
				Nparity[6][abs((h+k+l)%2)][bin] ++;
				
				if ( (h+k)%2 == 0 && (h+l)%2 == 0 && (k+l)%2 == 0 ) {
					Iparity[7][0][bin] += I_over_sigma;
					Nparity[7][0][bin] ++;
				}
				else {
					Iparity[7][1][bin] += I_over_sigma;
					Nparity[7][1][bin] ++;
				}
			}
		}
		
		for (int i=0; i<60; i++) {
			for (int j=0; j<2; j++) {
				for (int k=0; k<8; k++) {
					Itot[k][j] += Iparity[k][j][i];
					Ntot[k][j] += Nparity[k][j][i];			  
				}
			}
		}
		
		for (int j=0; j<2; j++) {
			for (int k=0; k<8; k++) {
				if ( Ntot[k][j] > 0 ) Itot[k][j] /= Ntot[k][j];
			}
		}
		
		for (int i=0; i<60; i++) {
			double res = maxres * (double(i) + 0.5)/double(nbins);
			for (int j=0; j<2; j++) {
				for (int k=0; k<8; k++) {
					if (Nparity[k][j][i] > 0) {
						Iparity[k][j][i] /= Nparity[k][j][i];
					}
				}
			}
			printf( " %5d%10.5f%7.2f%8d%5.1f%5.1f  %5.1f%5.1f  %5.1f%5.1f  %5.1f%5.1f  %5.1f%5.1f  %5.1f%5.1f  %5.1f%5.1f  %5.1f%5.1f\n", 
				   i+1, res, 1.0/sqrt(res), Nparity[0][0][i]+Nparity[0][1][i], 
				   Iparity[0][0][i], Iparity[0][1][i], Iparity[1][0][i], Iparity[1][1][i], Iparity[2][0][i], Iparity[2][1][i],
				   Iparity[3][0][i], Iparity[3][1][i], Iparity[4][0][i], Iparity[4][1][i], Iparity[5][0][i], Iparity[5][1][i],
				   Iparity[6][0][i], Iparity[6][1][i], Iparity[7][0][i], Iparity[7][1][i] );
			
		}
		printf( "\nTotals:                %8d%5.1f%5.1f  %5.1f%5.1f  %5.1f%5.1f  %5.1f%5.1f  %5.1f%5.1f  %5.1f%5.1f  %5.1f%5.1f  %5.1f%5.1f\n", 
			   Ntot[0][0]+Ntot[0][1], 
			   Itot[0][0], Itot[0][1], Itot[1][0], Itot[1][1], Itot[2][0], Itot[2][1], Itot[3][0], Itot[3][1],
			   Itot[4][0], Itot[4][1], Itot[5][0], Itot[5][1], Itot[6][0], Itot[6][1], Itot[7][0], Itot[7][1] );
	}



}