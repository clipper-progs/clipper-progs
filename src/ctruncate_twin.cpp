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
#include "ctruncate_utils.h"
#include "ccp4_fortran.h"
#include "csymlib.h"
#include "twinlaws.h"

namespace ctruncate {
	
	bool Ltest_driver(clipper::HKL_data<clipper::data32::I_sigI>& isig, CCP4Program& prog, bool debug)
	{
		// L test for twinning
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		
		bool itwin = false;
		
		double LT=0.0;
		double LT2=0.0;
		double NLT=0.0;
		std::vector<int> cdf(20,0);
		
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			if ( !isig[ih].missing() && !ih.hkl_class().centric() ) {
				clipper::HKL hkl = ih.hkl();
				int h = hkl.h();
				int k = hkl.k();
				int l = hkl.l();
				for ( int delta1 = -2; delta1 <= 2; delta1 += 2 ) {
					for ( int delta2 = -2; delta2 <= 2; delta2 += 2 ) {
						for ( int delta3 = -2; delta3 <= 2; delta3 += 2 ) {
							clipper::HKL hkl2;
							hkl2.h() = h+delta1;
							hkl2.k() = k+delta2;
							hkl2.l() = l+delta3;
							if ( !(delta1==0 && delta2==0 && delta3==0) ) {
								double I1 = isig[ih].I();
								double I2 = isig[hkl2].I();
								//double weight = 1.0/(isig[ih].sigI() + isig[jh].sigI());
								double weight = 1.0;
								double L = 0.0;
								//if ( I1 != 0.0 && I2 != 0.0 && I1/isig[ih].sigI() > 0.0 && I2/isig[hkl2].sigI() > 0.0 ) L = (I2-I1)/(I2+I1);
								if ( I1 != 0.0 && I2 != 0.0 ) L = (I2-I1)/(I2+I1);
								//printf("%f\n",L);
								if (fabs(L) < 1){
									LT += fabs(L)*weight;
									LT2 += L*L*weight;
									NLT += weight;
									for (int i=0;i<20;i++) {
										if ( fabs(L) < (double(i+1))/20.0 ) cdf[i]++;
									}
								}
							}
						}
					}
				}
			}
		}
		double Lav = LT/NLT;
		double L2av = LT2/NLT;
		//printf("Lav = %f  Untwinned 0.5 Perfect Twin 0.375\n",Lav);
		//printf("L2av = %f  Untwinned 0.333 Perfect Twin 0.200\n",L2av);
		if (Lav < 0.48) {
			printf("\nApplying the L test for twinning: (Padilla and Yeates Acta Cryst. D59 1124 (2003))\n");
			printf("L test suggests data is twinned\n");
			printf("L statistic = %6.3f  (untwinned 0.5 perfect twin 0.375)\n\n", Lav);
			itwin = true;
			printf("All data regardless of I/sigma(I) has been included in the L test\n");
			//printf("Data has been truncated at %6.2f A resolution\n",resopt);
		}
		
		if (!itwin) printf("No twinning detected\n\n");
		prog.summary_end();
		
		printf("$TABLE: L test for twinning:\n");
		printf("$GRAPHS");
		printf(": cumulative distribution function for |L|:0|1x0|1:1,2,3,4:\n");
		printf("$$ |L| Observed Expected_untwinned Expected_twinned $$\n$$\n");
		printf("0.000000 0.000000 0.000000 0.000000\n");
		
		for (int i=0;i<20;i++) {
			double x = (double(i+1))/20.0;
			printf("%f %f %f %f\n", x, double(cdf[i])/NLT, x, 0.5*x*(3.0-x*x)  );
		}
		printf("$$\n\n");
		return itwin;
	}
		
	bool Htest( clipper::HKL_data<clipper::data32::I_sigI>& isig, clipper::Mat33<int>& twinop, int scalefac, 
			   clipper::String s, CCP4Program& prog, bool debug )
	{
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		bool itwin = false;
		
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
			itwin = true;
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
		return itwin;
	}
	
	bool Htest_driver_fp(clipper::HKL_data<clipper::data32::I_sigI>& isig, CCP4Program& prog, bool debug)
	{
		bool itwin = false;
		
		clipper::Spacegroup spgr = isig.hkl_info().spacegroup();
		clipper::Cell      cell = isig.hkl_info().cell();
		
		const int scalefac = 12;           // scale factor for integer symops
		double sc_tol = 3.5;    // tolerance for determining candidate twinops, obliquity in degrees
		int lc = 48;             // maximum number of twinops that can be stored
		int nc = 0;
		int ivb = 0;
		int ierr = 0;
		int nc2;
		
		std::vector< clipper::Vec3<int> > trans;
		std::vector< clipper::Mat33<int> > rot;
		//std::vector< clipper::Mat33<int> > twin(lc);   // NB need to dimension these o/w output from fortran corrupted
		std::vector<double> score(lc);
		
		int nsymops = spgr.num_symops();
		//printf("nsymops = %d\n",nsymops);
		clipper::Grid g( 12, 12, 12 );
		for (int i=0; i!=nsymops; ++i) {
			clipper::Isymop isymop( spgr.symop(i), g );
			/*for (int j=0; j<3; j++) {
			 printf("%6d %6d %6d   %6d\n", isymop.rot()(j,0), isymop.rot()(j,1), isymop.rot()(j,2), isymop.trn()[j] );
			 }
			 printf("\n");*/
			trans.push_back( isymop.trn() );
			rot.push_back( isymop.rot() );
		}
		
		int vv[nsymops][3][3];
		int ww[nsymops][3];
		int uu_c[lc][3][3];
		double *sc_c = &score[0];
		double twin_cell[6];
		{ 
			twin_cell[0] = cell.a();
			twin_cell[1] = cell.b();
			twin_cell[2] = cell.c();
			twin_cell[3] = cell.alpha();
			twin_cell[4] = cell.beta();
			twin_cell[5] = cell.gamma();
			for (int i=0; i!=nsymops; ++i) for (int j=0 ; j != 3; ++j ) for (int k = 0 ; k != 3 ; ++k ) vv[i][j][k] = rot[i](j,k);
			for (int i=0; i!=nsymops; ++i) for (int j=0 ; j != 3; ++j ) ww[i][j] = trans[i][j];
		}
		
		yyy_cell2tg(twin_cell, sc_tol, nsymops, vv, ww, lc, nc, nc2, uu_c, sc_c, ivb, ierr);
		
		//printf("nc = %d  ierr = %d\n",nc,ierr);
		printf("First principles calculation of potential twinning operators using code by Andrey Lebedev:\n");
        if (!nc) 
            printf("First principles calculation has found no potential twinning operators\n\n");
        else
            printf("First principles calculation has found %d potential twinning operators\n\n", nc-1);
		
		
		if (nc > 1) {
			for (int k=1; k!=nc; ++k) {
				// Transpose matrix once because passed from Fortran to C, then a second time because Andrey's
				// convention is that (h,k,l) is postmultiplied by the twinning operator, whereas mine is that
				// (h,k,l) is premultiplied by the the twinning operator. Net result: don't do anything!
				clipper::Mat33<int> twinoper = clipper::Mat33<int>(uu_c[k][0][0],uu_c[k][0][1],uu_c[k][0][2],uu_c[k][1][0],uu_c[k][1][1]
																   ,uu_c[k][1][2],uu_c[k][2][0],uu_c[k][2][1],uu_c[k][2][2]);
				clipper::String s;
				MatrixToString(twinoper,s);
				std::cout << "Twinning operator: " << s << std::endl;
                std::cout << "      Obliquity of the twin operator is " << sc_c[k] << " degrees" << std::endl;
                std::cout << "      (above 2 degrees probably be ignored)" << std::endl;
				/*for (int i=0; i<3; i++) {
				 // Divide by 12 (scale factor for integer syops)
				 printf("%7.4f %7.4f %7.4f\n",double(twinoper(i,0))/12.0, double(twinoper(i,1))/12.0, double(twinoper(i,2))/12.0 );
				 }
				 printf("\n");*/
				itwin = Htest(isig, twinoper, scalefac, s, prog, debug);
			}
		}
		return itwin;
	}		
	
	bool Htest_driver_table(clipper::HKL_data<clipper::data32::I_sigI>& isig, CCP4Program& prog, bool debug)
	{
		bool itwin = false;
		
		clipper::Spacegroup spgr = isig.hkl_info().spacegroup();
		clipper::Cell      cell = isig.hkl_info().cell();
		
		int sg = spgr.spacegroup_number();
		
		CSym::CCP4SPG *spg1 = CSym::ccp4spg_load_by_ccp4_num(sg);
		char pointgroup[20];
		strcpy(pointgroup,spg1->point_group);
		
		printf("\n   Potential twinning operators found from tables:\n\n");
		clipper::Mat33<int> twinop(0,0,0,0,0,0,0,0,0);
		int scalefac = 1;
		clipper::String s;
		if ( (sg >= 75 && sg <= 80) || sg == 146 || (sg >= 168 && sg <= 173) || (sg >= 195 && sg <= 199) ) { 
			printf("Twinning operator k, h, -l\n");
			s = "k, h, -l";
			twinop(0,1) = 1;
			twinop(1,0) = 1;
			twinop(2,2) = -1;
			itwin = Htest ( isig, twinop, scalefac, s, prog, debug );
		}
		else if( sg >= 149 && sg <= 154 ) {
			printf("Twinning operator -h, -k, l\n");
			s = "-h, -k, l";
			twinop(0,0) = -1;
			twinop(1,1) = -1;
			twinop(2,2) = 1;
			itwin = Htest ( isig, twinop, scalefac, s, prog, debug );
		}
		else if( sg >= 143 && sg <= 145 ) {
			printf("Twinning operator k, h, -l\n");
			s = "k, h, -l";
			twinop(0,1) = 1;
			twinop(1,0) = 1;
			twinop(2,2) = -1;
			itwin = Htest ( isig, twinop, scalefac, s, prog, debug );
			
			printf("Twinning operator -k, -h, -l\n");
			s = "-k, -h, -l";
			twinop(0,1) = -1;
			twinop(1,0) = -1;
			itwin = Htest ( isig, twinop, scalefac, s, prog, debug );
			
			printf("Twinning operator -h, -k, l\n");
			s = "-h, -k, l";
			twinop(0,1) = 0;
			twinop(1,0) = 0;
			twinop(0,0) = -1;
			twinop(1,1) = -1;
			twinop(2,2) = 1;
			itwin = Htest ( isig, twinop, scalefac, s, prog, debug );
		}
		else if( !strcmp(pointgroup, "PG222") ) {
			//printf("PG222\n");
			// Can have pseudo-merohedral twinning in PG222 (orthorhombic) if a=b, b=c or c=a
			if ( fabs( 1.0 - cell.b()/cell.a() ) < 0.02 ) { 
				printf("Twinning operator k, h, -l\n");
				s = "k, h, -l";
				twinop(0,1) = 1;
				twinop(1,0) = 1;
				twinop(2,2) = -1;
				itwin = Htest ( isig, twinop, scalefac, s, prog, debug );
			}
			if ( fabs( 1.0 - cell.c()/cell.b() ) < 0.02 ) { 
				printf("Twinning operator -h, l, k\n");
				s = "-h, l, k";
				twinop(0,1) = 0;
				twinop(1,0) = 0;
				twinop(2,2) = 0;
				twinop(0,0) = -1;
				twinop(1,2) = 1;
				twinop(2,1) = 1;
				itwin = Htest ( isig, twinop, scalefac, s, prog, debug );
			}
			if ( fabs( 1.0 - cell.a()/cell.c() ) < 0.02 ) {
				printf("Twinning operator l, -k, h\n");
				s = "l, -k, h";
				twinop(0,0) = 0;
				twinop(1,2) = 0;
				twinop(2,1) = 0;
				twinop(1,1) = -1;
				twinop(0,2) = 1;
				twinop(2,0) = 1;
				itwin = Htest ( isig, twinop, scalefac, s, prog, debug );
			}
		}
		else if( !strcmp(pointgroup, "PG2") ) {
			// can have pseudomerohedral twinning in PG2 (monoclinic) if
			// beta = 90
			// a=c
			// cos(beta) = -a/2c, -c/a, -c/2a
			// sin(beta) = a/c
			if ( fabs( 1.0 - cell.a()/cell.c() ) < 0.02  || fabs( sin(cell.beta()) - cell.a()/cell.c() ) < 0.02 ) {
				printf("Twinning operator l, -k, h\n");
				s = "l, -k, h";
				twinop(1,1) = -1;
				twinop(0,2) = 1;
				twinop(2,0) = 1;
				itwin = Htest ( isig, twinop, scalefac, s, prog, debug );
			}
			
			if ( cell.beta() < clipper::Util::d2rad(93.0) ) {
				printf("Twinning operator -h, -k, l\n");
				s = "-h, -k, l";
				twinop(0,0) = -1;
				twinop(0,2) = 0;
				twinop(1,1) = -1;
				twinop(2,0) = 0;
				twinop(2,2) = 1;
				itwin = Htest ( isig, twinop, scalefac, s, prog, debug );
			}	 
			
			if ( fabs( cos(cell.beta()) + 0.5*cell.a()/cell.c() ) < 0.02 ) { 
				printf("Twinning operator -h, -k, h+l\n");
				s = "-h, -k, h+l";
				twinop(0,0) = -1;
				twinop(0,2) = 0;
				twinop(1,1) = -1;
				twinop(2,0) = 1;
				twinop(2,2) = 1;
				itwin = Htest ( isig, twinop, scalefac, s, prog, debug );
			}	  
			
			if ( fabs( cos(cell.beta()) + 0.5*cell.c()/cell.a() ) < 0.02 ) { 
				printf("Twinning operator h+l, -k, -l\n");
				s = "h+l, -k, -l";
				twinop(0,0) = 1;
				twinop(0,2) = 1;
				twinop(1,1) = -1;
				twinop(2,0) = 0;
				twinop(2,2) = -1;
				itwin = Htest ( isig, twinop, scalefac, s, prog, debug );
			}	  
			if ( fabs( cos(cell.beta()) + cell.c()/cell.a() ) < 0.02 ) { 
				printf("Twinning operator h+2l, -k, -l\n");
				s = "h+2l, -k, -l";
				twinop(0,0) = 1;
				twinop(0,2) = 2;
				twinop(1,1) = -1;
				twinop(2,0) = 0;
				twinop(2,2) = -1;
				itwin = Htest ( isig, twinop, scalefac, s, prog, debug );
			}
		}
		return itwin;
	}
	
}

