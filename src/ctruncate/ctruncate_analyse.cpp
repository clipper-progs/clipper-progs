//
//     CTRUNCATE
//     Copyright (C) 2006-2011 Norman Stein, Charles Ballard
//
//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Application.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//

#include "ctruncate_analyse.h"
#include "ctruncate_utils.h"
#include "intensity_target.h"
#include "ctruncate_wilson.h"
#include "cpsf_utils.h"
#include "best.h"
#include "best_rna.h"
#include <cmath>

#define ASSERT assert
#include <assert.h>

namespace ctruncate {
	
	int cumulative_plot(clipper::HKL_data<clipper::data32::I_sigI>& isig, clipper::HKL_data<clipper::data32::I_sigI>& Sigma)
	{
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		
		// construct cumulative distribution function for intensity (using Z rather than E)
		clipper::Range<double> intensity_range_centric;
		clipper::Range<double> intensity_range_acentric;
		// changed from jsig to isig
		
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			if ( !isig[ih].missing() ) {
				if ( ih.hkl_class().centric() ) intensity_range_centric.include( (isig[ih].I()/ih.hkl_class().epsilon() )/Sigma[ih].I() );
				else intensity_range_acentric.include( (isig[ih].I()/ih.hkl_class().epsilon() ) /Sigma[ih].I() );
			}
		}
		//printf("C2: %20.5f %20.5f\n",intensity_range_centric.max(),intensity_range_centric.min());
		//printf("A2: %20.5f %20.5f\n",intensity_range_acentric.max(),intensity_range_acentric.min());
		
		int ncentric = 0;
		clipper::Generic_ordinal intensity_ord_c;
		clipper::Generic_ordinal intensity_ord_a;
		intensity_ord_c.init( intensity_range_centric );
		intensity_ord_a.init( intensity_range_acentric );
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			if ( !isig[ih].missing() ) {
				if ( ih.hkl_class().centric() ) {
					intensity_ord_c.accumulate( ( isig[ih].I()/ih.hkl_class().epsilon() )/Sigma[ih].I()  );
					ncentric++;
				}
				else intensity_ord_a.accumulate( ( isig[ih].I()/ih.hkl_class().epsilon() )/Sigma[ih].I() );
			}
		}
		intensity_ord_c.prep_ordinal();
		intensity_ord_a.prep_ordinal();
		
		// theoretical values for cumulative intensity distribution
		double acen[51] = {0.0,
			0.0392106, 0.0768837, 0.1130796, 0.1478562, 0.1812692, 0.2133721, 0.2442163, 0.2738510, 0.3023237, 0.3296800,
			0.3559636, 0.3812166, 0.4054795, 0.4287909, 0.4511884, 0.4727076, 0.4933830, 0.5132477, 0.5323336, 0.5506710,
			0.5682895, 0.5852171, 0.6014810, 0.6171071, 0.6321206, 0.6465453, 0.6604045, 0.6737202, 0.6865138, 0.6988058,
			0.7106158, 0.7219627, 0.7328647, 0.7433392, 0.7534030, 0.7630722, 0.7723623, 0.7812881, 0.7898639, 0.7981035,
			0.8060200, 0.8136260, 0.8209339, 0.8279551, 0.8347011, 0.8411826, 0.8474099, 0.8533930, 0.8591416, 0.8646647};
		
		double cen[51] = {0.0,
			0.1585194, 0.2227026, 0.2709655, 0.3108435, 0.3452792, 0.3757939, 0.4032988, 0.4283924, 0.4514938, 0.4729107,
			0.4928775, 0.5115777, 0.5291583, 0.5457398, 0.5614220, 0.5762892, 0.5904133, 0.6038561, 0.6166715, 0.6289066,
			0.6406032, 0.6517983, 0.6625250, 0.6728131, 0.6826895, 0.6921785, 0.7013024, 0.7100815, 0.7185345, 0.7266783,
			0.7345289, 0.7421010, 0.7494079, 0.7564625, 0.7632764, 0.7698607, 0.7762255, 0.7823805, 0.7883348, 0.7940968,
			0.7996745, 0.8050755, 0.8103070, 0.8153755, 0.8202875, 0.8250491, 0.8296659, 0.8341433, 0.8384867, 0.8427008};
		
		double acen_twin[51] = {0.0,
			0.0030343, 0.0115132, 0.0245815, 0.0414833, 0.0615519, 0.0842006, 0.1089139, 0.1352404, 0.1627861, 0.1912079, 
			0.2202081, 0.2495299, 0.2789524, 0.3082868, 0.3373727, 0.3660750, 0.3942806, 0.4218963, 0.4488460, 0.4750691, 
			0.5005177, 0.5251562, 0.5489585, 0.5719077, 0.5939942, 0.6152149, 0.6355726, 0.6550744, 0.6737317, 0.6915590, 
			0.7085736, 0.7247951, 0.7402450, 0.7549459, 0.7689218, 0.7821971, 0.7947971, 0.8067470, 0.8180725, 0.8287987, 
			0.8389511, 0.8485543, 0.8576328, 0.8662106, 0.8743109, 0.8819565, 0.8891694, 0.8959710, 0.9023818, 0.9084218}; 
		
		
		printf("$TABLE: Cumulative intensity distribution:\n");
		printf("$GRAPHS");
		printf(": Cumulative intensity distribution (Acentric and centric):N:1,2,3,4,5,6:\n$$");
		printf(" Z Acent_theor Acent_twin Acent_obser Cent_theor Cent_obser $$\n$$\n");
		double x = 0.0;
		double deltax=0.04;
		for (int i=0; i<=50; i++) {
			if (ncentric) printf("%10.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n", x, acen[i], acen_twin[i], intensity_ord_a.ordinal(x), cen[i], intensity_ord_c.ordinal(x));
			else printf("%10.5f %8.5f %8.5f %8.5f %8.5f -\n", x, acen[i], acen_twin[i], intensity_ord_a.ordinal(x), cen[i]);
			x += deltax;
		}
		printf("$$\n\n");
		
		// count entries where acentric distribution lower than expected - sign of twinning
		x = 0.08;
		deltax = 0.12;
		int ntw = 0;
		for (int i=0;i<4; i++) {
			if ( (acen[3*i+2] - intensity_ord_a.ordinal(x))/acen[3*i+2] > 0.4 ) ntw ++;
			x += deltax;
		}
		return ntw;
	}
	
	int cumulative_plot(clipper::HKL_data<clipper::data32::I_sigI>& isig, clipper::ResolutionFn& Sigma)
	{
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		
		// construct cumulative distribution function for intensity (using Z rather than E)
		clipper::Range<double> intensity_range_centric;
		clipper::Range<double> intensity_range_acentric;
		// changed from jsig to isig
		
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			if ( !isig[ih].missing() ) {
				if ( ih.hkl_class().centric() ) intensity_range_centric.include( (isig[ih].I()/ih.hkl_class().epsilon() )/Sigma.f(ih) );
				else intensity_range_acentric.include( (isig[ih].I()/ih.hkl_class().epsilon() ) /Sigma.f(ih) );
			}
		}
		//printf("C2: %20.5f %20.5f\n",intensity_range_centric.max(),intensity_range_centric.min());
		//printf("A2: %20.5f %20.5f\n",intensity_range_acentric.max(),intensity_range_acentric.min());
		
		int ncentric = 0;
		clipper::Generic_ordinal intensity_ord_c;
		clipper::Generic_ordinal intensity_ord_a;
		intensity_ord_c.init( intensity_range_centric );
		intensity_ord_a.init( intensity_range_acentric );
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			if ( !isig[ih].missing() ) {
				if ( ih.hkl_class().centric() ) {
					intensity_ord_c.accumulate( ( isig[ih].I()/ih.hkl_class().epsilon() )/Sigma.f(ih)  );
					ncentric++;
				}
				else intensity_ord_a.accumulate( ( isig[ih].I()/ih.hkl_class().epsilon() )/Sigma.f(ih) );
			}
		}
		intensity_ord_c.prep_ordinal();
		intensity_ord_a.prep_ordinal();
		
		// theoretical values for cumulative intensity distribution
		double acen[51] = {0.0,
			0.0392106, 0.0768837, 0.1130796, 0.1478562, 0.1812692, 0.2133721, 0.2442163, 0.2738510, 0.3023237, 0.3296800,
			0.3559636, 0.3812166, 0.4054795, 0.4287909, 0.4511884, 0.4727076, 0.4933830, 0.5132477, 0.5323336, 0.5506710,
			0.5682895, 0.5852171, 0.6014810, 0.6171071, 0.6321206, 0.6465453, 0.6604045, 0.6737202, 0.6865138, 0.6988058,
			0.7106158, 0.7219627, 0.7328647, 0.7433392, 0.7534030, 0.7630722, 0.7723623, 0.7812881, 0.7898639, 0.7981035,
			0.8060200, 0.8136260, 0.8209339, 0.8279551, 0.8347011, 0.8411826, 0.8474099, 0.8533930, 0.8591416, 0.8646647};
		
		double cen[51] = {0.0,
			0.1585194, 0.2227026, 0.2709655, 0.3108435, 0.3452792, 0.3757939, 0.4032988, 0.4283924, 0.4514938, 0.4729107,
			0.4928775, 0.5115777, 0.5291583, 0.5457398, 0.5614220, 0.5762892, 0.5904133, 0.6038561, 0.6166715, 0.6289066,
			0.6406032, 0.6517983, 0.6625250, 0.6728131, 0.6826895, 0.6921785, 0.7013024, 0.7100815, 0.7185345, 0.7266783,
			0.7345289, 0.7421010, 0.7494079, 0.7564625, 0.7632764, 0.7698607, 0.7762255, 0.7823805, 0.7883348, 0.7940968,
			0.7996745, 0.8050755, 0.8103070, 0.8153755, 0.8202875, 0.8250491, 0.8296659, 0.8341433, 0.8384867, 0.8427008};
		
		double acen_twin[51] = {0.0,
			0.0030343, 0.0115132, 0.0245815, 0.0414833, 0.0615519, 0.0842006, 0.1089139, 0.1352404, 0.1627861, 0.1912079, 
			0.2202081, 0.2495299, 0.2789524, 0.3082868, 0.3373727, 0.3660750, 0.3942806, 0.4218963, 0.4488460, 0.4750691, 
			0.5005177, 0.5251562, 0.5489585, 0.5719077, 0.5939942, 0.6152149, 0.6355726, 0.6550744, 0.6737317, 0.6915590, 
			0.7085736, 0.7247951, 0.7402450, 0.7549459, 0.7689218, 0.7821971, 0.7947971, 0.8067470, 0.8180725, 0.8287987, 
			0.8389511, 0.8485543, 0.8576328, 0.8662106, 0.8743109, 0.8819565, 0.8891694, 0.8959710, 0.9023818, 0.9084218}; 
		
		
		printf("$TABLE: Cumulative intensity distribution:\n");
		printf("$GRAPHS");
		printf(": Cumulative intensity distribution (Acentric and centric):N:1,2,3,4,5,6:\n$$");
		printf(" Z Acent_theor Acent_twin Acent_obser Cent_theor Cent_obser $$\n$$\n");
		double x = 0.0;
		double deltax=0.04;
		for (int i=0; i<=50; i++) {
			if (ncentric) printf("%10.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n", x, acen[i], acen_twin[i], intensity_ord_a.ordinal(x), cen[i], intensity_ord_c.ordinal(x));
			else printf("%10.5f %8.5f %8.5f %8.5f %8.5f -\n", x, acen[i], acen_twin[i], intensity_ord_a.ordinal(x), cen[i]);
			x += deltax;
		}
		printf("$$\n\n");
		
		// count entries where acentric distribution lower than expected - sign of twinning
		x = 0.08;
		deltax = 0.12;
		int ntw = 0;
		for (int i=0;i<4; i++) {
			if ( (acen[3*i+2] - intensity_ord_a.ordinal(x))/acen[3*i+2] > 0.4 ) ntw ++;
			x += deltax;
		}
		return ntw;
	}
	
	void yorgo_modis_plot(clipper::HKL_data<clipper::data32::F_sigF>& fsig, float maxres, int nbins, CCP4Program& prog, clipper::U_aniso_orth uao )
	{
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		
        //get uao eigenvectors
		AnisoDirection<float> direct(uao);
		clipper::Mat33<float> e123;
			for (int i = 0 ; i !=3 ; ++i)
				for (int j = 0 ; j != 3 ; ++j)
            e123(i,j) = (direct.eigenVectors())[i][j];
		
		clipper::Cell cell = fsig.hkl_info().cell();
		clipper::Spacegroup spg = fsig.hkl_info().spacegroup();
		
		
		std::vector<float> somov(nbins,0.0);
		std::vector<float> somsdov(nbins,0.0);
		std::vector<int> numov(nbins,0);
		std::vector<float> enumov(nbins,0.0);
		
		float somdir[3][nbins];
		float somsddir[3][nbins];
		int numdir[3][nbins];
		float enumdir[3][nbins];
		
		for (int i=0;i<3;i++){
			for (int j=0;j<nbins;j++){
				somdir[i][j] = somsddir[i][j] = enumdir[i][j] = 0.0;
				numdir[i][j] = 0;
			}
		}
		
		int nzerosigma = 0;
		float cone = 30.0; //hardwired for now
		float ang;
		float cosang;
		
		for ( HRI ih = fsig.first(); !ih.last(); ih.next() ) {
			if ( !fsig[ih].missing() ) {
				// bin number different in C because arrays start at zero
				int bin = int( double(nbins) * ih.invresolsq() / maxres - 0.001);
				if (bin >= nbins || bin < 0) printf("Warning: (Modis) illegal bin number %d\n", bin);
				float epsiln = 1.0f/ih.hkl_class().epsilonc();
				
				for ( int jsym = 0; jsym != spg.num_primitive_symops() ; ++jsym ) {
					for (int friedal = 0 ; friedal != 2 ; ++friedal) {
						clipper::HKL ri = int(std::pow( -1.0f, float(friedal) ))*ih.hkl();
						clipper::HKL rj = ri.transform( spg.primitive_symop( jsym ) );
						
						clipper::Vec3<float> hc = e123*clipper::Vec3<float>(rj.coord_reci_orth(cell) );  //transpose into eigenspace
						
						for (int j=0;j!=3;++j) {
							cosang = fabs(hc[j])/sqrt(ih.invresolsq());
							// cosang can stray just past 1.0
							cosang = std::min(cosang, 1.0f);
							ang = acos(cosang);
							if ( ang < clipper::Util::d2rad(cone) ) {
								somdir[j][bin] += fsig[ih].f()*epsiln;
								if ( fsig[ih].sigf() > 0.0f ) somsddir[j][bin] += epsiln*fsig[ih].f()/fsig[ih].sigf();
								enumdir[j][bin] += epsiln;
								numdir[j][bin]++;
							}
						}
						somov[bin] += fsig[ih].f()*epsiln;
						if ( fsig[ih].sigf() > 0.0f ) somsdov[bin] += epsiln*fsig[ih].f()/fsig[ih].sigf();
						else nzerosigma++;
						enumov[bin] += epsiln;
						numov[bin]++;
					}
				}
			}
		}
		
		for (int i=0;i != nbins; ++i) {
			for (int j=0;j!=3;++j) {
				if (numdir[j][i] == 0) {
					somdir[j][i] = 0;
					somsddir[j][i] = 0;
				}
				else {
					somdir[j][i] /= enumdir[j][i];
					somsddir[j][i] /= enumdir[j][i];
				}
			}
			if (numov[i] == 0) {
				somov[i] = 0.0;
				somsdov[i] = 0.0;
			}
			else {
				somov[i] /= enumov[i];
				somsdov[i] /= enumov[i];
			}
		}
		
		if (nzerosigma > 0) {
			prog.summary_beg();
			printf("\nWARNING: ****  %d reflections have zero sigma ****\n\n", nzerosigma);
			prog.summary_end();
		}
		
		// calculate completeness
		std::vector<float> sumov(nbins,0.0);
		std::vector<float> summeas(nbins,0.0);
		std::vector<float> summeas1(nbins,0.0);
		std::vector<float> summeas2(nbins,0.0);
		std::vector<float> summeas3(nbins,0.0);
		std::vector<float> completeness(nbins,0.0);
		std::vector<float> completeness1(nbins,0.0);
		std::vector<float> completeness2(nbins,0.0);
		std::vector<float> completeness3(nbins,0.0);
		for ( HRI ih = fsig.first(); !ih.last(); ih.next() ) {
			// bin number different in C because arrays start at zero
			float mult = ih.hkl_class().epsilonc();
			int bin = int( double(nbins) * ih.invresolsq() / maxres - 0.001);
			//if (bin >= nbins || bin < 0) printf("Warning: (completeness) illegal bin number %d\n", bin);
			if ( bin < nbins && bin >= 0 ) sumov[bin] += mult;
			if ( !fsig[ih].missing() && bin < nbins && bin >= 0) {
				summeas[bin] += mult;
				float isigi = fsig[ih].f()/fsig[ih].sigf();
				if (isigi >= 1.0f ) summeas1[bin] += mult;
				if (isigi >= 2.0f ) summeas2[bin] += mult;
				if (isigi >= 3.0f ) summeas3[bin] += mult;
		}
		}
		for (int i=1; i!=nbins; ++i) {
			if (sumov[i] > 0.0) completeness[i] = summeas[i]/sumov[i];
			if (sumov[i] > 0.0) completeness1[i] = summeas1[i]/sumov[i];
			if (sumov[i] > 0.0) completeness2[i] = summeas2[i]/sumov[i];
			if (sumov[i] > 0.0) completeness3[i] = summeas3[i]/sumov[i];
		}
		
		
		printf("\n$TABLE: Structure amplitude statistics:\n");
		printf("$GRAPHS");
		printf(": Mn(F) v resolution:N:1,2,3,4,5:\n");
		printf(": Mn(F/sd) v resolution:N:1,6,7,8,9:\n");
		printf(": No. reflections v resolution:N:1,10,11,12,13:\n");
		printf(": Completeness v resolution:N:1,14,15,16,17:\n");
		printf("$$ 1/resol^2 Mn(F(d1)) Mn(F(d2)) Mn(F(d3)) Mn(F(ov) Mn(F/sd(d1)) Mn(F/sd(d2)) Mn(F/sd(d3)) Mn(F/sd(ov))");
		printf(" N(d1) N(d2) N(d3) N(ov) completeness sig1 sig2 sig3$$\n$$\n");
		
		
		for(int i=0;i<nbins;i++){
			double res = maxres*(double(i)+0.5)/double(nbins);
			printf("%10.6f %12.4e %12.4e %12.4e %12.4e ",res,somdir[0][i],somdir[1][i],somdir[2][i],somov[i]);
			printf("%12.4e %12.4e %12.4e %12.4e ",somsddir[0][i],somsddir[1][i],somsddir[2][i],somsdov[i]);
			printf("%8d %8d %8d %8d",numdir[0][i],numdir[1][i],numdir[2][i],numov[i]);
			printf("%8.4f %8.4f %8.4f %8.4f\n",completeness[i],completeness1[i],completeness2[i],completeness3[i]);
		}
		printf("$$\n\n");
	}
	
	void yorgo_modis_plot(clipper::HKL_data<clipper::data32::I_sigI>& isig, float maxres, int nbins, CCP4Program& prog, clipper::U_aniso_orth uao )
	{
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		
        //get uao eigenvectors
		AnisoDirection<float> direct(uao);
		clipper::Mat33<float> e123;
			for (int i = 0 ; i !=3 ; ++i)
				for (int j = 0 ; j != 3 ; ++j)
                e123(i,j) = (direct.eigenVectors())[i][j];
		
		clipper::Cell cell = isig.hkl_info().cell();
		clipper::Spacegroup spg = isig.hkl_info().spacegroup();
		
		
		std::vector<float> somov(nbins,0.0);
		std::vector<float> somsdov(nbins,0.0);
		std::vector<int> numov(nbins,0);
		std::vector<float> enumov(nbins,0.0);
		
		float somdir[3][nbins];
		float somsddir[3][nbins];
		int numdir[3][nbins];
		float enumdir[3][nbins];
		
		for (int i=0;i<3;i++){
			for (int j=0;j<nbins;j++){
				somdir[i][j] = somsddir[i][j] = enumdir[i][j] = 0.0;
				numdir[i][j] = 0;
			}
		}
		
		int nzerosigma = 0;
		float cone = 30.0; //hardwired for now
		float ang;
		float cosang;
		
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			if ( !isig[ih].missing() ) {
				// bin number different in C because arrays start at zero
				int bin = int( double(nbins) * ih.invresolsq() / maxres - 0.001);
				if (bin >= nbins || bin < 0) printf("Warning: (Modis) illegal bin number %d\n", bin);
				float epsiln = 1.0f/ih.hkl_class().epsilonc();
				
				for ( int jsym = 0; jsym != spg.num_primitive_symops() ; ++jsym ) {
					for (int friedal = 0 ; friedal != 2 ; ++friedal) {
						clipper::HKL ri = int(std::pow( -1.0f, float(friedal) ))*ih.hkl();
						clipper::HKL rj = ri.transform( spg.primitive_symop( jsym ) );
						
						clipper::Vec3<float> hc = e123*clipper::Vec3<float>(rj.coord_reci_orth(cell) );  //transpose into eigenspace
						
						for (int j=0;j!=3;++j) {
							cosang = fabs( hc[j] )/sqrt(ih.invresolsq());
							// cosang can stray just past 1.0
							cosang = std::min(cosang, 1.0f);
							ang = acos(cosang);
							if ( ang < clipper::Util::d2rad(cone) ) {
								somdir[j][bin] += isig[ih].I()*epsiln;
								if ( isig[ih].sigI() > 0.0f ) somsddir[j][bin] += epsiln*isig[ih].I()/isig[ih].sigI();
								enumdir[j][bin] += epsiln;
								numdir[j][bin]++;
							}
						}
						somov[bin] += isig[ih].I()*epsiln;
						if ( isig[ih].sigI() > 0.0f ) somsdov[bin] += epsiln*isig[ih].I()/isig[ih].sigI();
						else nzerosigma++;
						enumov[bin] += epsiln;
						numov[bin]++;
					}
				}
			}
		}
		
		for (int i=0;i != nbins; ++i) {
			for (int j=0;j!=3;++j) {
				if (numdir[j][i] == 0) {
					somdir[j][i] = 0;
					somsddir[j][i] = 0;
				}
				else {
					somdir[j][i] /= enumdir[j][i];
					somsddir[j][i] /= enumdir[j][i];
				}
			}
			if (numov[i] == 0) {
				somov[i] = 0.0;
				somsdov[i] = 0.0;
			}
			else {
				somov[i] /= enumov[i];
				somsdov[i] /= enumov[i];
			}
		}
		
		if (nzerosigma > 0) {
			prog.summary_beg();
			printf("\nWARNING: ****  %d reflections have zero sigma ****\n\n", nzerosigma);
			prog.summary_end();
		}
		
		// calculate completeness
		std::vector<float> sumov(nbins,0.0);
		std::vector<float> summeas(nbins,0.0);
		std::vector<float> summeas1(nbins,0.0);
		std::vector<float> summeas2(nbins,0.0);
		std::vector<float> summeas3(nbins,0.0);
		std::vector<float> completeness(nbins,0.0);
		std::vector<float> completeness1(nbins,0.0);
		std::vector<float> completeness2(nbins,0.0);
		std::vector<float> completeness3(nbins,0.0);
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			// bin number different in C because arrays start at zero
			float mult = ih.hkl_class().epsilonc();
			int bin = int( double(nbins) * ih.invresolsq() / maxres - 0.001);
			//if (bin >= nbins || bin < 0) printf("Warning: (completeness) illegal bin number %d\n", bin);
			if ( bin < nbins && bin >= 0 ) sumov[bin] += mult;
			if ( !isig[ih].missing() && bin < nbins && bin >= 0) {
				summeas[bin] += mult;
				float isigi = isig[ih].I()/isig[ih].sigI();
				if (isigi >= 1.0f ) summeas1[bin] += mult;
				if (isigi >= 2.0f ) summeas2[bin] += mult;
				if (isigi >= 3.0f ) summeas3[bin] += mult;
			}
		}
		for (int i=1; i!=nbins; ++i) {
			if (sumov[i] > 0.0) completeness[i] = summeas[i]/sumov[i];
			if (sumov[i] > 0.0) completeness1[i] = summeas1[i]/sumov[i];
			if (sumov[i] > 0.0) completeness2[i] = summeas2[i]/sumov[i];
			if (sumov[i] > 0.0) completeness3[i] = summeas3[i]/sumov[i];
		}
		
		
		printf("\n$TABLE: Intensity statistics:\n");
		printf("$GRAPHS");
		printf(": Mn(I) v resolution:N:1,2,3,4,5:\n");
		printf(": Mn(I/sd) v resolution:N:1,6,7,8,9:\n");
		printf(": No. reflections v resolution:N:1,10,11,12,13:\n");
		printf(": Completeness v resolution:N:1,14,15,16,17:\n");
		printf("$$ 1/resol^2 Mn(I(d1)) Mn(I(d2)) Mn(I(d3)) Mn(I(ov) Mn(I/sd(d1)) Mn(I/sd(d2)) Mn(I/sd(d3)) Mn(I/sd(ov))");
		printf(" N(d1) N(d2) N(d3) N(ov) completeness sig1 sig2 sig3$$\n$$\n");
		
		
		for(int i=0;i<nbins;i++){
			double res = maxres*(double(i)+0.5)/double(nbins);
			printf("%10.6f %12.4e %12.4e %12.4e %12.4e ",res,somdir[0][i],somdir[1][i],somdir[2][i],somov[i]);
			printf("%12.4e %12.4e %12.4e %12.4e ",somsddir[0][i],somsddir[1][i],somsddir[2][i],somsdov[i]);
			printf("%8d %8d %8d %8d",numdir[0][i],numdir[1][i],numdir[2][i],numov[i]);
			printf("%8.4f %8.4f %8.4f %8.4f\n",completeness[i],completeness1[i],completeness2[i],completeness3[i]);
		}
		printf("$$\n\n");
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
	
	float PattPeak::operator()(clipper::HKL_data<clipper::data32::I_sigI>& Sigma)
	{
		int nprm = 60;
		const clipper::HKL_info& hklinf = Sigma.hkl_info();
		std::vector<double> params( nprm, 1.0 );
		clipper::BasisFn_spline basis_fo( Sigma, nprm, 2.0 );
		TargetFn_meanInth<clipper::data32::I_sigI> target_fo( Sigma, 1);
		clipper::ResolutionFn resol_fn( hklinf, basis_fo, target_fo, params );
		
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
		std::vector<clipper::ftype> weights(_nbins,1.0f);
		std::vector<clipper::ftype> x(_nbins);
		std::vector<clipper::ftype> y(_nbins);
		
		for( int i = 0 ; i != _nbins ; ++i) {
			float dist = (float(i)+0.5f)*_width/float(_nbins);
			x[i] = 0.25*dist*dist;
			y[i] = std::log(1.0f/_patterson[i]);
		}
		
        clipper::ftype a, b, siga, sigb;
		
		straight_line_fit(x,y,weights,_nbins,a,b,siga,sigb);
		
		_P000 = std::exp(-b);
		
		float btemp = 0.25*a;
		
		_sigma = std::sqrt(1.0f/std::abs(2.0f*btemp) );
		
		return;	
	}
	
	//******YorgoModis***************************************************************************
	
	template <class D> void YorgoModis<D>::operator() (clipper::HKL_data<D>& isig, clipper::Resolution reso )
	{
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		
		_t = type(isig[isig.first()]); //set our type
		
		//if resolution not set use observed
		_reso = reso;
        if (_reso.is_null() ) _reso = isig.hkl_info().resolution();
		
        //get uao eigenvectors
		AnisoDirection<float> direct(_uao);
			for (int i = 0 ; i !=3 ; ++i)
				for (int j = 0 ; j != 3 ; ++j)
                _e123(i,j) = (direct.eigenVectors())[i][j];
		
		clipper::Cell cell = isig.hkl_info().cell();
		clipper::Spacegroup spg = isig.hkl_info().spacegroup();
		
		int _nzerosigma = 0;
		float cone = 30.0; //hardwired for now
		float ang;
		float cosang;
		
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			if ( !isig[ih].missing() ) {
				// bin number different in C because arrays start at zero
				int bin = int( double(_nbins) * ih.invresolsq() / _reso.invresolsq_limit() - 0.001);
				if (bin >= _nbins || bin < 0) printf("Warning: (Modis) illegal bin number %d\n", bin);
				float epsiln = 1.0f/ih.hkl_class().epsilonc();
				clipper::ftype Ival=obs(isig[ih]);
				clipper::ftype Isig=sigobs(isig[ih]);
				for ( int jsym = 0; jsym != spg.num_primitive_symops() ; ++jsym ) {
					for (int friedal = 0 ; friedal != 2 ; ++friedal) {
						clipper::HKL ri = int(std::pow( -1.0f, float(friedal) ))*ih.hkl();
						clipper::HKL rj = ri.transform( spg.primitive_symop( jsym ) );
						
						clipper::Vec3<clipper::ftype> hc = _e123*clipper::Vec3<clipper::ftype>(rj.coord_reci_orth(cell) );  //transpose into eigenspace

						
						for (int j=0;j!=3;++j) {
							int jn = j*_nbins+bin;
							cosang = fabs( hc[j] )/sqrt(ih.invresolsq());
							// cosang can stray just past 1.0
							cosang = std::min(cosang, 1.0f);
							ang = acos(cosang);
							if ( ang < clipper::Util::d2rad(cone) ) {
								_somdir[jn] += Ival*epsiln;
								if ( Isig > 0.0f ) _somsddir[jn] += epsiln*Ival/Isig;
								_enumdir[jn] += epsiln;
								++_numdir[jn];
							}
						}
						_somov[bin] += Ival*epsiln;
						if ( Isig > 0.0f ) _somsdov[bin] += epsiln*Ival/Isig;
						else _nzerosigma++;
						_enumov[bin] += epsiln;
						_numov[bin]++;
					}
				}
			}
		}
		
		for (int i=0;i != _nbins; ++i) {
			for (int j=0;j!=3;++j) {
				int jn = j*_nbins+i;
				if (_numdir[jn] == 0) {
					_somdir[jn] = 0.0;
					_somsddir[jn] = 0.0;
				}
				else {
					_somdir[jn] /= _enumdir[jn];
					_somsddir[jn] /= _enumdir[jn];
				}
			}
			if (_numov[i] == 0) {
				_somov[i] = 0.0;
				_somsdov[i] = 0.0;
			}
			else {
				_somov[i] /= _enumov[i];
				_somsdov[i] /= _enumov[i];
			}
		}
		return;
	}
	
	template<class D> void YorgoModis<D>::plot() 
	{
		if (_t == I ) {
			printf("\n$TABLE: Intensity statistics:\n");
			printf("$GRAPHS");
			printf(": Mn(I) v resolution:N:1,2,3,4,5:\n");
			printf(": Mn(I/sd) v resolution:N:1,6,7,8,9:\n");
			printf(": No. reflections v resolution:N:1,10,11,12,13:\n");
			printf("$$ 1/resol^2 Mn(I(d1)) Mn(I(d2)) Mn(I(d3)) Mn(I(ov) Mn(I/sd(d1)) Mn(I/sd(d2)) Mn(I/sd(d3)) Mn(I/sd(ov))");
			printf(" N(d1) N(d2) N(d3) N(ov)$$\n$$\n");
		} else {
			printf("\n$TABLE: Anisotropy analysis (Yorgo Modis):\n");
			printf("$GRAPHS");
			printf(": Mn(F) v resolution:N:1,2,3,4,5:\n");
			printf(": Mn(F/sd) v resolution:N:1,6,7,8,9:\n");
			printf(": No. reflections v resolution:N:1,10,11,12,13:\n");
			printf("$$ 1/resol^2 Mn(F(d1)) Mn(F(d2)) Mn(F(d3)) Mn(F(ov) Mn(F/sd(d1)) Mn(F/sd(d2)) Mn(F/sd(d3)) Mn(F/sd(ov))");
			printf(" N(d1) N(d2) N(d3) N(ov)$$\n$$\n");
		}
		
		for(int i=0;i!=_nbins;++i){
			double res = _reso.invresolsq_limit()*(double(i)+0.5)/double(_nbins);
			printf("%10.6f %12.4e %12.4e %12.4e %12.4e ",res,_somdir[i],_somdir[_nbins+i],_somdir[2*_nbins+i],_somov[i]);
			printf("%12.4e %12.4e %12.4e %12.4e ",_somsddir[i],_somsddir[_nbins+i],_somsddir[2*_nbins+i],_somsdov[i]);
			printf("%8d %8d %8d %8d\n",_numdir[i],_numdir[_nbins+i],_numdir[2*_nbins+i],_numov[i]);
		}
		printf("$$\n\n");
	}
	
	//******Completeness***************************************************************************
	
	template<class D> void Completeness<D>::operator() (clipper::HKL_data<D>& isig, clipper::Resolution reso )
	{
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		
		_t = type(isig[isig.first()]); //set our type
		
		//if resolution not set use observed
		_reso = reso;
        if (_reso.is_null() ) _reso = isig.hkl_info().resolution();
		
		// calculate completeness
		std::vector<clipper::ftype> sumov(_nbins,0.0);
		std::vector<clipper::ftype> summeas(_nbins,0.0);
		
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			// bin number different in C because arrays start at zero
			float mult = ih.hkl_class().epsilonc();
			int bin = int( double(_nbins) * ih.invresolsq() / _reso.invresolsq_limit() - 0.001);
			//if (bin >= _nbins || bin < 0) printf("Warning: (completeness) illegal bin number %d\n", bin);
			if ( bin < _nbins && bin >= 0 ) {
				sumov[bin] += mult;
				if ( !isig[ih].missing() ||
				!clipper::Util::is_nan(obs(isig[ih])) || !clipper::Util::is_nan(sigobs(isig[ih])) ) {
					_compsig[bin] += mult;
					clipper::ftype Ival=obs(isig[ih]);
					clipper::ftype Isig=sigobs(isig[ih]);
					clipper::ftype isigi = Ival/Isig;
					if (isigi >= 1.0f ) _compsig1[bin] += mult;
					if (isigi >= 2.0f ) _compsig2[bin] += mult;
					if (isigi >= 3.0f ) _compsig3[bin] += mult;
					_standard[bin] += (isigi > 0.01 ) ? mult/(isigi) : mult/100.0 ;
                                        summeas[bin] += mult;
				}
			}
		}
		for (int i=0; i !=_nbins; ++i) {
			if (sumov[i] > 0.0) {
				_compsig[i] /= sumov[i];
				_compsig1[i] /= sumov[i];
				_compsig2[i] /= sumov[i];
				_compsig3[i] /= sumov[i];
				_standard[i] /= summeas[i];
			}
		}
	}
		
	template<class D> void Completeness<D>::plot() 
	{
		if (_t == I ) {
			printf("\n$TABLE: Intensity Completeness:\n");
			printf("$GRAPHS");
			printf(": Completeness v resolution:N:1,2,3,4,5,6:\n");
			printf("$$ 1/resol^2 completeness sig1 sig2 sig3 standard$$\n$$\n");
		} else {
			printf("\n$TABLE: Structure factor completeness:\n");
			printf("$GRAPHS");
			printf(": Completeness v resolution:N:1,2,3,4,5,6:\n");
			printf("$$ 1/resol^2 completeness sig1 sig2 sig3 standard$$\n$$\n");
		}
		
		for(int i=0;i!=_nbins;++i){
			double res = _reso.invresolsq_limit()*(double(i)+0.5)/double(_nbins);
			printf("%10.6f %8.4f %8.4f %8.4f %8.4f %8.4f\n",res, _compsig[i],_compsig1[i],_compsig2[i],_compsig3[i]
				   ,_standard[i]);
		}
		printf("$$\n\n");
	}
	
	template <class T> const std::vector<clipper::Coord_frac>& tNCS<T>::operator() (clipper::HKL_data<clipper::datatypes::I_sigI<T> >& I, clipper::Resolution r)
	{
		intensity = &I;
		reso = r;
		
		clipper::HKL_info hklinf(I.hkl_info());
		clipper::Cell cell = hklinf.cell();
		clipper::Spacegroup spgr = hklinf.spacegroup();
		
		// check for pseudo translation (taken from cpatterson)
		// get Patterson spacegroup
		clipper::HKL_info hklp;
		clipper::Spacegroup
		pspgr( clipper::Spgr_descr( spgr.generator_ops().patterson_ops() ) );
		hklp.init( pspgr, cell, reso, true );
		
		// make patterson coeffs
		clipper::HKL_data<clipper::datatypes::F_phi<T> > fphi( hklp );
		for ( clipper::HKL_data_base::HKL_reference_index  ih = fphi.first(); !ih.last(); ih.next() ) {
			clipper::datatypes::I_sigI<T> i = I[ih.hkl()];
			if ( !i.missing() ) {
				fphi[ih].f() = i.I();
				fphi[ih].phi() = 0.0 ;
			}
		}
		
		// make grid if necessary
		clipper::Grid_sampling grid( pspgr, cell, reso );
		
		// make xmap
		clipper::Xmap<float> patterson( pspgr, cell, grid );
		patterson.fft_from( fphi );
		
		
		// use Charles's stuff to find peaks
		PeakSearch pksch;                      // peak search object
		PeakInterp pkinterp;                   // peak interpolation methods
		
		int npeak = 5;
		
		std::vector<int> ppks = pksch( patterson );
		
		T top_peak = patterson.get_data( ppks[0] );
		int i = 0;
		T pval(0.0);
		
		do {
			T next_peak = patterson.get_data( ppks[++i] );
			clipper::Coord_frac c0 = patterson.coord_of( ppks[i] ).coord_frac(grid);
			T ratio = next_peak/top_peak;
			T dist2 = std::sqrt(c0.lengthsq(cell) );
			// look for peaks > 20% of origin peak and at least 14A distant from origin
			// precentage estimate is Zwartz CCP4 Newsletter 42
                        if (dist2 > 14.0 ) {
			const T aval = 0.0679;
			const T bval = 3.56;
			pval = (1.0 - std::exp(-std::pow(ratio/(aval*(T(1.0)-ratio)),-bval)) )*100.0;
			if (pval < 1.0) {
				peaks.push_back(c0);
				peak_prob.push_back(pval);
				peak_height.push_back(ratio);
			}
			
		} while ( pval < 1.0 );
	
		return peaks;
	}
	
	template <class T> void tNCS<T>::summary() {
		printf("\n\nTRANSLATIONAL NCS:\n");
		if ( peaks.size() && peak_prob[0] < 1.0 ) { 
			clipper::Coord_frac c0 = peaks[0];
			printf("Translational NCS has been detected at (%6.3f, %6.3f, %6.3f).\n  The probability based on peak ratio is %5.2f%% \n that this is by chance (with resolution limited to %5.2f A). \n", c0[0],c0[1],c0[2],peak_prob[0],reso.limit() );
			printf("This will have a major impact on the twinning estimates and effectiveness of the truncate procedure\n");
			for (int i = 0; i != peaks.size() ; ++i )
				printf("Peak %d Ratio = %5.2f with Peak Vector = (%6.3f, %6.3f, %6.3f)\n",i+1,peak_height[i],peaks[i][0],peaks[i][1],peaks[i][2]);
		}
		else {
			printf("No translational NCS detected (with resolution limited to %5.2f A)\n", reso.limit() );
		}
	}
	
	
	template <class T> AnisoPlot<T>::AnisoPlot(clipper::U_aniso_orth& uao)
	{
		int steps = 60;
		clipper::ftype lev[] = { 0.5, 0.25, 0.125};
		std::vector<clipper::ftype> levels(lev,lev+sizeof(lev)/sizeof(clipper::ftype));
		// get eigenvalues of aniso_U
			{
				clipper::Matrix<T> m(3,3);
				for (int i = 0 ; i !=3 ; ++i)
					for (int j = 0 ; j != 3 ; ++j)
						m(i,j) = uao(i,j);
				_eigen = m.eigen();
				for (int i = 0 ; i !=3 ; ++i)
					_eigen[i] = 2.0/(clipper::Util::twopi2()*_eigen[i]);
				for (int i = 0 ; i !=3 ; ++i)
					for (int j = 0 ; j != 3 ; ++j)
						_e123(i,j) = m(i,j);
			}
			clipper::ftype angle, sigmau, sigmav;
			ellipse(_eigen[0],_eigen[1],0.0,angle,sigmau,sigmav);
			for (int l=0 ; l != levels.size() ; ++l ) 
				_isoline1.push_back(isoline(0.0,0.0,sigmau,sigmav, angle,levels[l], steps));
			ellipse(_eigen[0],_eigen[2],0.0,angle,sigmau,sigmav);
			for (int l=0 ; l != levels.size() ; ++l ) 
				_isoline2.push_back(isoline(0.0,0.0,sigmau,sigmav, angle,levels[l], steps));
			ellipse(_eigen[1],_eigen[2],0.0,angle,sigmau,sigmav);
			for (int l=0 ; l != levels.size() ; ++l ) 
				_isoline3.push_back(isoline(0.0,0.0,sigmau,sigmav, angle,levels[l], steps));
		}
	
	template <class T> AnisoPlot<T>::AnisoPlot(clipper::ftype scale, clipper::U_aniso_orth& uao)
	{
		int steps = 60;
		clipper::ftype lev[] = { 3, 2, 1};
		std::vector<clipper::ftype> levels(lev,lev+sizeof(lev)/sizeof(clipper::ftype));
		for (int i = 0; i != levels.size() ; ++i ) levels[i] /= scale;
		// get eigenvalues of aniso_U
		{
			clipper::Matrix<T> m(3,3);
			for (int i = 0 ; i !=3 ; ++i)
				for (int j = 0 ; j != 3 ; ++j)
					m(i,j) = uao(i,j);
			_eigen = m.eigen();
			for (int i = 0 ; i !=3 ; ++i)
				_eigen[i] = 2.0/(clipper::Util::twopi2()*_eigen[i]);
			for (int i = 0 ; i !=3 ; ++i)
				for (int j = 0 ; j != 3 ; ++j)
					_e123(i,j) = m(i,j);
		}
		clipper::ftype angle, sigmau, sigmav;
		ellipse(_eigen[0],_eigen[1],0.0,angle,sigmau,sigmav);
		for (int l=0 ; l != levels.size() ; ++l ) 
			_isoline1.push_back(isoline(0.0,0.0,sigmau,sigmav, angle,levels[l], steps));
		ellipse(_eigen[0],_eigen[2],0.0,angle,sigmau,sigmav);
		for (int l=0 ; l != levels.size() ; ++l ) 
			_isoline2.push_back(isoline(0.0,0.0,sigmau,sigmav, angle,levels[l], steps));
		ellipse(_eigen[1],_eigen[2],0.0,angle,sigmau,sigmav);
		for (int l=0 ; l != levels.size() ; ++l ) 
			_isoline3.push_back(isoline(0.0,0.0,sigmau,sigmav, angle,levels[l], steps));
	}
	
		/* ellipse */
		template <class T> void AnisoPlot<T>::loggraph() {
			std::stringstream x1,x2,x3;
			x1 << "("  << _e123(0,0) << "h," << _e123(1,0) << "k," << _e123(2,0) << "l)";
			x2 << "("  << _e123(0,1) << "h," << _e123(1,1) << "k," << _e123(2,1) << "l)";
			x3 << "("  << _e123(0,2) << "h," << _e123(1,2) << "k," << _e123(2,2) << "l)";
			
			clipper::ftype maxv(-99.0), minv(99);
			// get plot extremes
			{
			int l = _isoline1.size()-1;
			for(int i=0;i!=_isoline1[l].size();++i) {
				maxv = std::max(maxv,std::max(_isoline1[l][i][0],_isoline1[l][i][1]) );
				minv = std::min(minv,std::min(_isoline1[l][i][0],_isoline1[l][i][1]) );
			}
			l = _isoline2.size()-1;
			for(int i=0;i!=_isoline2[l].size();++i) {
				maxv = std::max(maxv,std::max(_isoline2[l][i][0],_isoline2[l][i][1]) );
				minv = std::min(minv,std::min(_isoline2[l][i][0],_isoline2[l][i][1]) );
			}
			l = _isoline3.size()-1;
			for(int i=0;i!=_isoline1[l].size();++i) {
				maxv = std::max(maxv,std::max(_isoline3[l][i][0],_isoline3[l][i][1]) );
				minv = std::min(minv,std::min(_isoline3[l][i][0],_isoline3[l][i][1]) );
			}
			}
			printf("\n$TABLE: Anisotropy:\n");
			printf("$SCATTER");
			printf(": Falloff plane 1:%f|%fx%f|%f:1,2:\n",minv,maxv,minv,maxv);
			printf(": Falloff plane 2:%f|%fx%f|%f:3,4:\n",minv,maxv,minv,maxv);
			printf(": Falloff plane 3:%f|%fx%f|%f:5,6:\n",minv,maxv,minv,maxv);
			printf("$$ %s %s %s %s %s %s$$\n$$\n",x1.str().c_str(),x2.str().c_str(),x1.str().c_str(),
				   x3.str().c_str(),x2.str().c_str(),x3.str().c_str());	
			
			// ellipse level and step
			for(int l=0;l!=_isoline1.size(); ++l)
				for(int i=0;i!=_isoline1[l].size();++i)
					printf("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", 
						   _isoline1[l][i][0],_isoline1[l][i][1],
						   _isoline2[l][i][0],_isoline2[l][i][1],
						   _isoline3[l][i][0],_isoline3[l][i][1]);
			
			printf("$$\n\n");
		}
		
		
		//want a point, but use Coord_orth for now
		template <class T> clipper::Coord_orth AnisoPlot<T>::point(clipper::ftype offsetx, clipper::ftype offsety, clipper::ftype sigmau, 
																clipper::ftype sigmav, clipper::ftype angleuv, clipper::ftype theta, clipper::ftype frac)
		{
			clipper::ftype f = std::sqrt(2*std::log(1.0/frac));
			clipper::ftype sinuv = std::sin(angleuv);
			clipper::ftype cosuv = std::cos(angleuv);

			clipper::ftype s = std::sin(theta);
			clipper::ftype c = std::cos(theta);		 
			// theta equations
			clipper::ftype x = offsetx + f*(sigmau * c * cosuv - sigmav * s * sinuv);
			clipper::ftype y= offsety + f*(sigmau * c * sinuv + sigmav * s * cosuv);
			return clipper::Coord_orth(x,y,0.0);
		}
		
		// plot at frac of height
	template <class T> std::vector<clipper::Coord_orth> AnisoPlot<T>::isoline(clipper::ftype offsetx, clipper::ftype offsety, clipper::ftype sigmau, 
													  clipper::ftype sigmav, clipper::ftype angleuv, clipper::ftype frac, int steps)
		{
			std::vector<clipper::Coord_orth> points(steps);
			clipper::ftype f = std::sqrt(2*std::log(1.0/frac));
			clipper::ftype sinuv = std::sin(angleuv);
			clipper::ftype cosuv = std::cos(angleuv);
			
			for (int i = 0; i != steps ;++i) {
				clipper::ftype angle = 2.0*clipper::Util::pi()*clipper::ftype(i)/clipper::ftype(steps) ;
				clipper::ftype s = std::sin(angle);
				clipper::ftype c = std::cos(angle);		
				// theta equations
				points[i] = clipper::Coord_orth(offsetx + f*(sigmau * c * cosuv - sigmav * s * sinuv),
												offsety + f*(sigmau * c * sinuv + sigmav * s * cosuv),
												0.0);
			}
			
			return points;
		}
		
		// calculate the ellipse parameters
		template <class T> void AnisoPlot<T>::ellipse(clipper::ftype sig1, clipper::ftype sig2, clipper::ftype cov, 
													  clipper::ftype& angle, clipper::ftype& sigmau, clipper::ftype& sigmav)		 
		{
			angle = ( std::abs(sig1 - sig2 ) < 1.0e-5) ? 0.0 : std::atan(2.0*cov*sig1*sig2/(sig1*sig1-sig2*sig2));
			//tan(2alpha) = 2*cov(x,y)*sig(x)*sig(y)/(sig(x)**2-sig(y)**2)
			
			clipper::ftype t1 = sig1*sig1*sig2*sig2*(1.0-cov*cov);
			clipper::ftype s = std::sin(angle);
			clipper::ftype c = std::cos(angle);
			clipper::ftype t2 = 2*cov*sig1*sig2*s*c;
			
			sigmav = std::sqrt(t1/(sig2*sig2*c*c-t2+sig1*sig1*s*s));
			sigmau = std::sqrt(t1/(sig2*sig2*s*s+t2+sig1*sig1*c*c));
			return;
		}

    //---------Calculate anisotropy-------------------------------------


    template <class SCALER, class DATA, class T> void AnisoCorr<SCALER,DATA,T>::calc(clipper::HKL_data<DATA>& isig, bool protein, bool rna)
    {
        const clipper::HKL_info& hklinf = isig.hkl_info();
        clipper::Cell cell(hklinf.cell());
        clipper::Spacegroup spgr(hklinf.spacegroup());
		clipper::HKL_data<DATA> Ibest(hklinf);
		
        //generate reference scattering curve (initially protein only)
		if (protein) {
        Scattering scat;
        clipper::ftype totalscatter = spgr.num_symops()*scat(cell, spgr);
        
		for ( clipper::HKL_data_base::HKL_reference_index ih = Ibest.first(); !ih.last(); ih.next() ) {
            T reso = ih.invresolsq();
            Ibest[ih] = DATA(ih.hkl_class().epsilon()*totalscatter*ctruncate::BEST(reso), 1.0f);
		} // scale against BEST
        } else {
			Scattering scat(Scattering::NUCLEIC);
			clipper::ftype totalscatter = spgr.num_symops()*scat(cell, spgr);
			
			for ( clipper::HKL_data_base::HKL_reference_index ih = Ibest.first(); !ih.last(); ih.next() ) {
				T reso = ih.invresolsq();
				Ibest[ih] = DATA(ih.hkl_class().epsilon()*totalscatter*ctruncate::BEST_rna(reso), 1.0f);
			} // scale against BEST
			
		}
        _iscale( Ibest, isig );

        return ;
    }
    
    //! calculate anisotropy correction
    template <class SCALER, class DATA, class T> 
    const clipper::U_aniso_orth& AnisoCorr<SCALER,DATA,T>::operator()(clipper::HKL_data<DATA>& observed, bool protein, bool rna )
    {
        calc(observed,protein,rna);
        return _iscale.u_aniso_orth(Scaling::F);
    }
    
    template <class SCALER, class DATA, class T> 
    const clipper::U_aniso_orth& AnisoCorr<SCALER,DATA,T>::u_aniso_orth( Scaling::TYPE t ) const
    {
        if ( t == Scaling::I ) return _iscale.u_aniso_orth(Scaling::I);
        else return _iscale.u_aniso_orth(Scaling::F);
    }
    
    //---------Calculate anisotropy eigenvalues------------------------

    //! Calculate eigenvalues
    
    template <class T> AnisoDirection<T>::AnisoDirection(clipper::U_aniso_orth& uao)
    {
        // Eigenvalue calculation
        _uao = &uao;
        clipper::Matrix<T> mat( 3, 3, 0.0 );
        for (int i=0; i !=3; ++i) {
            for (int j=0; j !=3 ; ++j) {
                mat(i,j) = uao(i,j);
            }
        }
        std::vector<T> v = mat.eigen( true );
        
        _max = 0.0;
        for (int i=0 ; i!=3 ; ++i) 
            if ( v[i] > _max ) _max = v[i];
        
        std::vector<int> close(3,-1);
        
        for (int i=0 ; i!=3 ; ++i) { // loop a*, b*, c*
            clipper::ftype max = 0.0;
            int jj = -1;
            for (int j=0; j!=3; ++j) { // loop vectors
                if (close[0] == j || close[1] == j ) continue;
                if (std::fabs(mat(j,i) > max) ) {
                    max = std::fabs(mat(j,i));
                    jj = j;
                }
            }
            close[i] = jj;
        }
        // close[i] for i -> a*, b*, c* is closest vector
        for (int i=0; i!=3 ;++i) { // loop a*, b*, c*
            _eigenvalues.push_back(v[close[i]]);
            _eigenvectors.push_back(clipper::Vec3<T>(mat(0,close[i]),mat(1,close[i]),mat(2,close[i])) );
        }
        ASSERT( _eigenvalues.size() == 3 );
        return ;
    }
    
    // run on scaled on merged data, however the best stats will be from the unmerged set.
    // get from aimless
    template<class T> AnomStats<T>::AnomStats(clipper::HKL_data<clipper::datatypes::J_sigJ_ano<T> >& isig_ano, int nbins)
    {
        typedef clipper::HKL_data_base::HKL_reference_index HRI;
        
        _nbins = nbins;
        std::vector<int> sumov(nbins,0.0);
		std::vector<int> summeas(nbins,0.0);
        std::vector<clipper::ftype> meandI(nbins,0.0);
        std::vector<clipper::ftype> meandIsigdI(nbins,0.0);
        std::vector<clipper::ftype> meanI(nbins,0.0);
        
        
        clipper::ftype maxres = isig_ano.hkl_info().resolution().invresolsq_limit();
        
		for (HRI ih=isig_ano.first() ; !ih.last() ; ih.next() ) {
            int eps = ih.hkl_class().epsilonc();
            int bin = int( double(nbins) * ih.invresolsq() / maxres - 0.5);
			sumov[bin] += eps;
            if ( ih.hkl_class().centric() ) {
                //no anomalous information
                if ( !clipper::Util::is_nan(obs_pl(isig_ano[ih]) )  &&  !clipper::Util::is_nan(obs_mi(isig_ano[ih]) ) ) {
                    meanI[bin] += 0.5*eps*(obs_pl(isig_ano[ih])+obs_mi(isig_ano[ih]));
                } else if ( !clipper::Util::is_nan(obs_pl(isig_ano[ih]) ) ) {
                    meanI[bin] += eps*obs_pl(isig_ano[ih]);
                } else if ( !clipper::Util::is_nan(obs_mi(isig_ano[ih]) ) ) {
                    meanI[bin] += eps*obs_mi(isig_ano[ih]);
                }
            } else {
				if ( !clipper::Util::is_nan(obs_pl(isig_ano[ih]) )  &&  !clipper::Util::is_nan(obs_mi(isig_ano[ih]) )  ) {
					clipper::ftype Ip = obs_pl(isig_ano[ih]);
					clipper::ftype sp = sigobs_pl(isig_ano[ih]);
					clipper::ftype Im = obs_mi(isig_ano[ih]);
					clipper::ftype sm = sigobs_mi(isig_ano[ih]);
					clipper::ftype dI = std::fabs(Ip - Im);
					clipper::ftype ds = std::sqrt(sp*sp+sm*sm);
					if ( dI/ds >= 3.0 && Ip/sp >= 3.0 && Im/sm > 3.0 ) {
						summeas[bin] += eps;
					}
					meandI[bin] += eps*dI;
					meandIsigdI[bin] += eps*dI/ds;
					meanI[bin] += 0.5*eps*(obs_pl(isig_ano[ih])+obs_mi(isig_ano[ih]));
				} else if ( !clipper::Util::is_nan(obs_pl(isig_ano[ih]) ) ) {
                    meanI[bin] += eps*obs_pl(isig_ano[ih]);
				} else if ( !clipper::Util::is_nan(obs_mi(isig_ano[ih]) ) ) {
                    meanI[bin] += eps*obs_mi(isig_ano[ih]);
                }
            }
        }
        
        for (int i= 0; i != nbins ; ++i ) {
			meanI[i] /= (float) sumov[i];
            meandI[i] /= (float) sumov[i];
            meandIsigdI[i] /= (float) sumov[i];
        }
        
        //assume values decrease monatomically.  Cut at measurability 5%
        //Dauter Acta D62 (2006) 867
        //Zwart Acta D61 (2005) 1437
        clipper::Range<float> meas_limit;
        for (int i1 = 0; i1 != nbins ; ++i1 ) {
            if ( float(summeas[i1])/float(sumov[i1]) > 0.05 )
				meas_limit.include(maxres*(float(i1)+0.5)/float(nbins));
        } 
        
        //assume values decrease monatomically.  Cut at DeltaAnom at 1.3
        //Dauter Acta D62 (2006) 867
        //Schneider Acta D58 (2002) 1772
		clipper::Range<float> anom_limit;
        for (int i1 = 0  ; i1 != nbins ; ++i1 ) {
            if ( meandIsigdI[i1] > 1.3 ) 
				anom_limit.include(maxres*(float(i1)+0.5)/float(nbins));
        }
        
        //assume values decrease monatomically.  Cut at deltaI/I 0.6%
        //Zwart Acta D61 (2005) 1437
        //Wang Methods Enzymol 115 (1985) 90
        clipper::Range<float> wang_limit;
        for (int i1 = 0  ; i1 != nbins ; ++i1 ) {
            if ( meandI[i1]/meanI[i1] > 0.006 ) 
				wang_limit.include(maxres*(float(i1)+0.5)/float(nbins));;
        }
        
        std::cout << "Estimated limits of anomalous signal" << std::endl;
        
        std::cout << "      Wang limit (deltaI/I) > 0.6% : " << 1.0/std::sqrt(wang_limit.max() ) << " A " << std::endl;
        std::cout << "      anomalous limit (deltaI/sig) > 1.3 : " << 1.0/std::sqrt(anom_limit.max() ) << " A " << std::endl;
        std::cout << "      measurability limit (Nanon/Nov) > 5% : " << 1.0/std::sqrt(meas_limit.max() ) << " A " << std::endl;
        std::cout << "  These calculations are performed using scaled and merged data.  More accurate estimates of the limit of the anomalous signal can be obtained using scaled and unmerged data in the half dataset correlation calculation of aimless. " << std::endl;
        
        printf("\n$TABLE: Intensity anomalous analysis:\n");
		printf("$GRAPHS");
		printf(": Mn(dI) v resolution:N:1,2:\n");
        printf(": Mn(dI/sigdI) v resolution:N:1,3:\n");
		printf(": Mn(dI/I) v resolution:N:1,4:\n");
		printf(": Mesurability v resolution:N:1,5:\n");
		printf("$$ 1/resol^2 Mn(dI) Mn(dI/sigdI)) Mn(dI/I) measurability$$\n$$\n");
		for(int i=0;i!=nbins;++i){
			double res = maxres*(double(i)+0.5)/double(nbins);
			printf("%10.6f %12.4e %12.4e %12.4e %12.4e\n",res,meandI[i],meandIsigdI[i],meandI[i]/meanI[i],float(summeas[i])/float(sumov[i]));
		}
		printf("$$\n\n");
        
    }
    
    template<class T> AnomStats<T>::AnomStats(clipper::HKL_data<clipper::datatypes::G_sigG_ano<T> >& isig_ano, int nbins)
    {
        typedef clipper::HKL_data_base::HKL_reference_index HRI;
        
        _nbins = nbins;
        std::vector<int> sumov(nbins,0.0);
		std::vector<int> summeas(nbins,0.0);
        std::vector<clipper::ftype> meandI(nbins,0.0);
        std::vector<clipper::ftype> meandIsigdI(nbins,0.0);
        std::vector<clipper::ftype> meanI(nbins,0.0);
        
        
        clipper::ftype maxres = isig_ano.hkl_info().resolution().invresolsq_limit();
        
		for (HRI ih=isig_ano.first() ; !ih.last() ; ih.next() ) {
            int eps = ih.hkl_class().epsilon();
            int bin = int( double(nbins) * ih.invresolsq() / maxres - 0.5);
			sumov[bin] += eps;
            if ( ih.hkl_class().centric() ) {
                //no anomalous information
                if ( !clipper::Util::is_nan(obs_pl(isig_ano[ih]) )  &&  !clipper::Util::is_nan(obs_mi(isig_ano[ih]) ) ) {
                    meanI[bin] += 0.5*(obs_pl(isig_ano[ih])+obs_mi(isig_ano[ih]));
                } else if ( !clipper::Util::is_nan(obs_pl(isig_ano[ih]) ) ) {
                    meanI[bin] += obs_pl(isig_ano[ih]);
                } else if ( !clipper::Util::is_nan(obs_mi(isig_ano[ih]) ) ) {
                    meanI[bin] += obs_mi(isig_ano[ih]);
                }
            } else {
				if ( !clipper::Util::is_nan(obs_pl(isig_ano[ih]) )  &&  !clipper::Util::is_nan(obs_mi(isig_ano[ih]) )  ) {
					clipper::ftype Ip = obs_pl(isig_ano[ih]);
					clipper::ftype sp = sigobs_pl(isig_ano[ih]);
					clipper::ftype Im = obs_mi(isig_ano[ih]);
					clipper::ftype sm = sigobs_mi(isig_ano[ih]);
					clipper::ftype dI = std::fabs(Ip - Im);
					clipper::ftype ds = std::sqrt(sp*sp+sm*sm);
					if ( dI/ds >= 3.0 && Ip/sp >= 3.0 && Im/sm > 3.0 ) {
						summeas[bin] += eps;
					}
					meandI[bin] += dI;
					meandIsigdI[bin] += dI/ds;
					meanI[bin] += 0.5*(obs_pl(isig_ano[ih])+obs_mi(isig_ano[ih]));
				} else if ( !clipper::Util::is_nan(obs_pl(isig_ano[ih]) ) ) {
                    meanI[bin] += obs_pl(isig_ano[ih]);
				} else if ( !clipper::Util::is_nan(obs_mi(isig_ano[ih]) ) ) {
                    meanI[bin] += obs_mi(isig_ano[ih]);
                }
            }
        }		
        
        for (int i= 0; i != nbins ; ++i ) {
            meandI[i] /= (float) sumov[i];
            meandIsigdI[i] /= (float) sumov[i];
        }
        
        //assume values decrease monatomically.  Cut at measurability 5%
        //Dauter Acta D62 (2006) 867
        //Zwart Acta D61 (2005) 1437
        clipper::ftype meas_limit;
        int i1 = 0;
        for (  ; i1 != nbins ; ++i1 ) {
            if ( float(summeas[i1])/float(sumov[i1]) < 0.05 ) break;
        }
        if ( i1 == nbins ) {
            meas_limit = maxres;
        } else {
            meas_limit = maxres*(float(i1)+0.5)/float(nbins);
        }
        
        //assume values decrease monatomically.  Cut at DeltaAnom at 1.3
        //Dauter Acta D62 (2006) 867
        //Schneider Acta D58 (2002) 1772
        clipper::ftype anom_limit;
        i1 = 0;
        for (  ; i1 != nbins ; ++i1 ) {
            if ( meandIsigdI[i1] < 1.3 ) break;
        }
        if ( i1 == nbins ) {
            anom_limit = maxres;
        } else {
            anom_limit = maxres*(float(i1)+0.5)/float(nbins);
        }
        
        //assume values decrease monatomically.  Cut at deltaI/I 0.6%
        //Zwart Acta D61 (2005) 1437
        //Wang Methods Enzymol 115 (1985) 90
        clipper::ftype wang_limit;
        i1 = 0;
        for (  ; i1 != nbins ; ++i1 ) {
            if ( meandI[i1]/meanI[i1] < 0.006 ) break;
        }
        if ( i1 == nbins ) {
            wang_limit = maxres;
        } else {
            wang_limit = maxres*(float(i1)+0.5)/float(nbins);
        }
        
        std::cout << "Estimated limits of anomalous signal" << std::endl;
        
        std::cout << "      Wang limit (deltaF/F) > 0.6% : " << 1.0/std::sqrt(wang_limit) << " A " << std::endl;
        std::cout << "      anomalous limit (deltaF/sig) > 1.3 : " << 1.0/std::sqrt(anom_limit) << " A " << std::endl;
        std::cout << "      measurability limit (Nanon/Nov) > 5% : " << 1.0/std::sqrt(meas_limit) << " A " << std::endl;
        std::cout << "  These calculations are performed using scaled and merged data.  More accurate estimates of the limit of the anomalous signal can be obtained using scaled and unmerged data in the half dataset correlation calculation of aimless. " << std::endl;
        
        printf("\n$TABLE: Structure factor anomalous analysis:\n");
		printf("$GRAPHS");
		printf(": Mn(dF) v resolution:N:1,2:\n");
        printf(": Mn(dF/sigdF) v resolution:N:1,3:\n");
		printf(": Mn(dF/F) v resolution:N:1,4:\n");
		printf(": Mesurability v resolution:N:1,5:\n");
		printf("$$ 1/resol^2 Mn(dI) Mn(dI/sigdI)) Mn(dI/I) measurability$$\n$$\n");
		for(int i=0;i!=nbins;++i){
			double res = maxres*(double(i)+0.5)/double(nbins);
			printf("%10.6f %12.4e %12.4e %12.4e %12.4e\n",res,meandI[i],meandIsigdI[i],meandI[i]/meanI[i],float(summeas[i])/float(sumov[i]));
		}
		printf("$$\n\n");
    }
	
	//******ResoCorrels***************************************************************************
	
	template<class D1, class D2> void ResoCorrel<D1,D2>::operator() (clipper::HKL_data<D1>& isig1, clipper::HKL_data<D2>& isig2, 
																	 clipper::Range<clipper::ftype> reso )
	{
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		
		_t = type(isig1[isig1.first()]); //set our type
		
		//if resolution not set use observed
		_reso = reso;		
		if (reso.range() < 0.0 ) {
			_reso = isig2.hkl_info().invresolsq_range();
		}
		
		clipper::ftype num_symops = isig1.hkl_info().spacegroup().num_symops();
		
		// calculate completeness
		clipper::ftype eof(0.0), ecf(0.0), eoof(0.0), eccf(0.0), eocf(0.0), nf(0.0);
		std::vector<clipper::ftype> v_eof(_nbins,0.0);
		std::vector<clipper::ftype> v_ecf(_nbins,0.0);
		std::vector<clipper::ftype> v_eoof(_nbins,0.0);
		std::vector<clipper::ftype> v_eccf(_nbins,0.0);
		std::vector<clipper::ftype> v_eocf(_nbins,0.0);
		std::vector<clipper::ftype> v_nf(_nbins,0.0);
		
		for ( HRI ih = isig1.first(); !ih.last(); ih.next() ) {
			// bin number different in C because arrays start at zero
			if ( _reso.contains(ih.invresolsq() )) {
				clipper::ftype mult = num_symops/ih.hkl_class().epsilonc();
				
				int bin = int( double(_nbins) * (ih.invresolsq()- _reso.min()) / _reso.range());
				if ( bin < _nbins && bin >= 0 ) {
  					if ( ih.hkl() != clipper::HKL::zero() &&
						(!clipper::Util::is_nan(obs(isig1[ih])) || !clipper::Util::is_nan(sigobs(isig1[ih]) ) ) &&
						(!clipper::Util::is_nan(obs(isig2[ih.hkl()])) || !clipper::Util::is_nan(sigobs(isig2[ih.hkl()]) ) ) ) {
						clipper::ftype ff = obs(isig1[ih]);
						clipper::ftype foof = obs(isig2[ih.hkl()]);
						nf += mult;
						eof += foof*mult;
						eoof += foof * foof*mult;
						ecf += ff*mult;
						eocf += foof* ff*mult;
						eccf += ff * ff*mult;
						v_nf[bin] += mult;
						v_eof[bin] += foof*mult;
						v_eoof[bin] += foof * foof*mult;
						v_ecf[bin] += ff*mult;
						v_eocf[bin] += foof* ff*mult;
						v_eccf[bin] += ff * ff*mult;
					}
                }
			}
		}
 		for (int i=0; i !=_nbins; ++i) {
			_comp[i] = ( v_nf[i] * v_eocf[i] - v_eof[i] * v_ecf[i] ) 
					/ sqrt ( ( v_nf[i] * v_eccf[i] - v_ecf[i] * v_ecf[i] ) * ( v_nf[i] * v_eoof[i] - v_eof[i] * v_eof[i] ) );
		}
		_cc = ( nf * eocf - eof * ecf ) / sqrt ( ( nf * eccf - ecf * ecf ) * ( nf * eoof - eof * eof ) );
	}
	
	template<class D1,class D2> void ResoCorrel<D1,D2>::plot() 
	{
		printf("\n$TABLE: Resolution dependent correlation coefficient (CC=%5.4f):\n",_cc);
		printf("$GRAPHS");
		printf(": CC v resolution:N:1,2:\n");
		printf("$$ 1/resol^2 CC$$\n$$\n");
		for(int i=0;i!=_nbins;++i){
			double res = _reso.range()*(double(i)+0.5)/double(_nbins)+_reso.min();
			printf("%10.6f %8.4f \n",res, ( clipper::Util::is_nan(_comp[i]) ) ? 0.0 : _comp[i]);
		}
		printf("$$\n\n");
	}
	
    //---------Instantiate templates------------------------------------
	template class AnisoPlot<clipper::ftype32>;
	template class AnisoPlot<clipper::ftype64>;
    
    template class AnisoCorr<ctruncate::Iscale_logLikeAniso<clipper::ftype32>, clipper::datatypes::F_sigF<clipper::ftype32>,clipper::ftype32>;
    template class AnisoCorr<ctruncate::Iscale_logLikeAniso<clipper::ftype64>,clipper::datatypes::F_sigF<clipper::ftype64>,clipper::ftype64>;
    template class AnisoCorr<ctruncate::Iscale_logLikeAniso<clipper::ftype32>,clipper::datatypes::I_sigI<clipper::ftype32>,clipper::ftype32>;
    template class AnisoCorr<ctruncate::Iscale_logLikeAniso<clipper::ftype64>,clipper::datatypes::I_sigI<clipper::ftype64>,clipper::ftype64>;
    template class AnisoCorr<ctruncate::Iscale_wilsonAniso<clipper::ftype32>,clipper::datatypes::F_sigF<clipper::ftype32>,clipper::ftype32>;
    template class AnisoCorr<ctruncate::Iscale_wilsonAniso<clipper::ftype64>,clipper::datatypes::F_sigF<clipper::ftype64>,clipper::ftype64>;
    template class AnisoCorr<ctruncate::Iscale_wilsonAniso<clipper::ftype32>,clipper::datatypes::I_sigI<clipper::ftype32>,clipper::ftype32>;
    template class AnisoCorr<ctruncate::Iscale_wilsonAniso<clipper::ftype64>,clipper::datatypes::I_sigI<clipper::ftype64>,clipper::ftype64>;

    template class AnisoDirection<clipper::ftype32>;
    template class AnisoDirection<clipper::ftype64>;
    
    template class AnomStats<clipper::ftype32>;
    template class AnomStats<clipper::ftype64>;
	
	template class Completeness<clipper::datatypes::F_sigF<clipper::ftype32> >;
	template class Completeness<clipper::datatypes::I_sigI<clipper::ftype32> >;
	
	template class tNCS<clipper::ftype32>;
	template class tNCS<clipper::ftype64>;
	
	template class YorgoModis<clipper::datatypes::F_sigF<clipper::ftype32> >;
	template class YorgoModis<clipper::datatypes::I_sigI<clipper::ftype32> >;
	
	template class ResoCorrel<clipper::datatypes::I_sigI<clipper::ftype32>,clipper::datatypes::I_sigI<clipper::ftype32> >;
	template class ResoCorrel<clipper::datatypes::F_sigF<clipper::ftype32>,clipper::datatypes::F_sigF<clipper::ftype32> >;
	}
