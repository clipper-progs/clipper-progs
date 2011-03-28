//
//     CTRUNCATE
//     Copyright (C) 2006-2011 Norman Stein, Charles Ballard
//
//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Application.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//
 
#include "ctruncate_twin.h"

namespace ctruncate {
	
	void Htest( clipper::HKL_data<clipper::data32::I_sigI>& isig, clipper::Mat33<int>& twinop, int &itwin, int scalefac, 
			   clipper::String s, CCP4Program& prog, bool debug )
	{
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		double HT=0.0;
		double HT2=0.0;
		double NT=0.0;
		std::vector<int> pdf(50,0);
		std::vector<int> cdf(20,0);
		clipper::Vec3<int> jhkl, jhkl2;
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			if ( !isig[ih].missing() && !ih.hkl_class().centric() ) {
				clipper::HKL hkl = ih.hkl();
				jhkl[0] = hkl.h();
				jhkl[1] = hkl.k();
				jhkl[2] = hkl.l();
				clipper::HKL twin;
				jhkl2 = twinop*jhkl;
				twin.h() = jhkl2[0]/scalefac;
				twin.k() = jhkl2[1]/scalefac;
				twin.l() = jhkl2[2]/scalefac;
				//printf("%d %d %d %d %d %d\n",hkl.h(),hkl.k(),hkl.l(),twin.h(),twin.k(),twin.l());
				//printf("%d %d %d %d %d %d\n",jhkl[0],jhkl[1],jhkl[2],jhkl2[0],jhkl2[1],jhkl2[2]);
				if (!isig[twin].missing()) {
					double I1 = isig[ih].I();
					double I2 = isig[twin].I();
					//double weight = 1.0/(isig[ih].sigI() + isig[twin].sigI());
					double weight = 1.0;
					double H = 0.0;
					if ( I1 != 0.0 && I2 != 0.0) H = (I2-I1)/(I2+I1);
					if (fabs(H) < 1){
						HT += fabs(H)*weight;
						HT2 += H*H*weight;
						NT += weight;
						//fprintf(newfileout,"%10.4f %10.4f %10.4f %10.4f %8d\n",I1,I2,H,HT,NT);
						for (int i=0;i<20;i++) {
							if ( fabs(H) < (double(i+1))/20.0 ) cdf[i]++;
						}
					}
					// Britton test
					if ( I1 > 0.0 && I2 > 0.0 && I2 < I1 ) {
						double B = I2/(I1+I2);
						int bin = int(100.0*B);
						if (bin >= 0 && bin < 50) pdf[bin]++;
					}
				}
			}
		}
		double Hav = HT/NT;
		double H2av = HT2/NT;
		double alpha = 0.5-Hav;
		//printf("alpha = %f\n",alpha);
		printf("Applying the H test for twinning: (Yeates Acta Cryst. A44 142 (1980))\n");
		if (alpha > 0.05) {
			printf("H test suggests data is twinned\n");
			printf("Twinning fraction = %5.2f\n\n",alpha);
			itwin = 1;
			prog.summary_end();
			printf("$TABLE: H test for twinning (operator %s):\n", s.c_str() );
			printf("$GRAPHS");
			printf(": cumulative distribution function for |H|:0|1x0|1:1,2,3,4,5,6,7:\n");
			printf("$$ |H| 0.4 0.3 0.2 0.1 0.0 Observed $$\n$$\n");
			printf("0.000000 0.0 0.0 0.0 0.0 0.0 0.000000\n");
			
			for (int i=0;i<19;i++) {
				double x = (double(i+1))/20.0;
				//printf("%f %f %f %f %f %f %f\n", x, double(cdf[i])/NT, 5.0*x, 2.5*x, 1.667*x, 1.25*x, x  );
				printf("%f  -   -   -   -   -  %f\n", x, double(cdf[i])/NT);
			}
			printf("1.000000 5.0 2.5 1.667 1.25 1.0 1.0\n");
			printf("$$\n\n");
			prog.summary_beg();
		}
		else {
			printf("No twinning detected for this twinning operator\n\n");
		}
		alpha = 0.5*(1.0 - sqrt(3.0*H2av));
		//printf("alpha = %f\n",alpha);
		if (debug) {
			FILE *newfileout;
			newfileout=fopen("Htest.txt","w");
			for (int i=0;i<50;i++) {
				fprintf(newfileout,"%f %d\n", 0.01*double(i), pdf[i]);
			}
			fclose(newfileout);
		}
		return;
	}
	
}

