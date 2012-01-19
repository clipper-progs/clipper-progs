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
	
	void yorgo_modis_plot(clipper::HKL_data<clipper::data32::F_sigF>& fsig, float maxres, int nbins, CCP4Program& prog)
	{
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		
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
						
						clipper::Coord_reci_orth hc(rj.coord_reci_orth(cell) );
						
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
		std::vector<float> completeness(nbins,0.0);
		for ( HRI ih = fsig.first(); !ih.last(); ih.next() ) {
			// bin number different in C because arrays start at zero
			int bin = int( double(nbins) * ih.invresolsq() / maxres - 0.001);
			//if (bin >= nbins || bin < 0) printf("Warning: (completeness) illegal bin number %d\n", bin);
			if ( bin < nbins && bin >= 0 ) sumov[bin] += 1.0;
			if ( !fsig[ih].missing() && bin < nbins && bin >= 0) summeas[bin] += 1.0;
		}
		for (int i=1; i<nbins; i++) {
			if (sumov[i] > 0.0) completeness[i] = summeas[i]/sumov[i];
		}
		
		
		printf("\n$TABLE: Anisotropy analysis (Yorgo Modis):\n");
		printf("$GRAPHS");
		printf(": Mn(F) v resolution:N:1,2,3,4,5:\n");
		printf(": Mn(F/sd) v resolution:N:1,6,7,8,9:\n");
		printf(": No. reflections v resolution:N:1,10,11,12,13:\n");
		printf(": Completeness v resolution:N:1,14:\n");
		printf("$$ 1/resol^2 Mn(F(d1)) Mn(F(d2)) Mn(F(d3)) Mn(F(ov) Mn(F/sd(d1)) Mn(F/sd(d2)) Mn(F/sd(d3)) Mn(F/sd(ov))");
		printf(" N(d1) N(d2) N(d3) N(ov) completeness$$\n$$\n");
		
		
		for(int i=0;i<nbins;i++){
			double res = maxres*(double(i)+0.5)/double(nbins);
			printf("%10.6f %12.4e %12.4e %12.4e %12.4e ",res,somdir[0][i],somdir[1][i],somdir[2][i],somov[i]);
			printf("%12.4e %12.4e %12.4e %12.4e ",somsddir[0][i],somsddir[1][i],somsddir[2][i],somsdov[i]);
			printf("%8d %8d %8d %8d",numdir[0][i],numdir[1][i],numdir[2][i],numov[i]);
			printf("%8.4f\n",completeness[i]);
		}
		printf("$$\n\n");
	}
	
	
	void yorgo_modis_plot(clipper::HKL_data<clipper::data32::I_sigI>& isig, float maxres, int nbins, CCP4Program& prog)
	{
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		
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
						
						clipper::Coord_reci_orth hc(rj.coord_reci_orth(cell) );
						
						for (int j=0;j!=3;++j) {
							cosang = fabs(hc[j])/sqrt(ih.invresolsq());
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
		std::vector<float> weights(_nbins,1.0f);
		std::vector<float> x(_nbins);
		std::vector<float> y(_nbins);
		
		for( int i = 0 ; i != _nbins ; ++i) {
			float dist = (float(i)+0.5f)*_width/float(_nbins);
			x[i] = 0.25*dist*dist;
			y[i] = std::log(1.0f/_patterson[i]);
		}
		
		float a, b, siga, sigb;
		
		straight_line_fit(x,y,weights,_nbins,a,b,siga,sigb);
		
		_P000 = std::exp(-a);
		
		float btemp = 0.25*b;
		
		_sigma = std::sqrt(1.0f/std::abs(2.0f*btemp) );
		
		return;	
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
            int eps = ih.hkl_class().epsilon();
            int bin = int( double(nbins) * ih.invresolsq() / maxres - 0.5);
            if ( ih.hkl_class().centric() ) {
                //no anomalous information
                sumov[bin] += eps;
                if ( !clipper::Util::is_nan(obs_pl(isig_ano[ih]) )  &&  !clipper::Util::is_nan(obs_mi(isig_ano[ih]) ) ) {
                    meanI[bin] += 0.5*(obs_pl(isig_ano[ih])+obs_mi(isig_ano[ih]));
                } else if ( !clipper::Util::is_nan(obs_pl(isig_ano[ih]) ) ) {
                    meanI[bin] += obs_pl(isig_ano[ih]);
                } else {
                    meanI[bin] += obs_mi(isig_ano[ih]);
                }
            } else if ( !clipper::Util::is_nan(obs_pl(isig_ano[ih]) )  &&  !clipper::Util::is_nan(obs_pl(isig_ano[ih]) ) ) {
                clipper::ftype Ip = obs_pl(isig_ano[ih]);
                clipper::ftype sp = sigobs_pl(isig_ano[ih]);
                clipper::ftype Im = obs_mi(isig_ano[ih]);
                clipper::ftype sm = sigobs_mi(isig_ano[ih]);
                clipper::ftype dI = Ip - Im;
                clipper::ftype ds = std::sqrt(sp*sp+sm*sm);
                if ( dI/ds >= 3.0 && Ip/sp >= 3.0 && Im/sm > 3.0 ) {
                    summeas[bin] += eps;
                }
                meandI[bin] += dI;
                meandIsigdI[bin] += dI/ds;
                meanI[bin] += 0.5*(obs_pl(isig_ano[ih])+obs_mi(isig_ano[ih]));
            } else  {
                sumov[bin] += eps;
                if ( !clipper::Util::is_nan(obs_pl(isig_ano[ih]) ) ) {
                    meanI[bin] += obs_pl(isig_ano[ih]);
                } else {
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
 
        std::cout << "      Wang limit (deltaI/I) > 0.6% : " << 1.0/std::sqrt(wang_limit) << " A " << std::endl;
        std::cout << "      anomalous limit (deltaI/sig) > 1.3 : " << 1.0/std::sqrt(anom_limit) << " A " << std::endl;
        std::cout << "      measurability limit (Nanon/Nov) > 5% : " << 1.0/std::sqrt(meas_limit) << " A " << std::endl;
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
			printf("%10.6f %12.4e %12.4e %12.4e %12.4e\n",res,meandI[i],meandIsigdI[i],meandI[i]/meanI[i],float(summeas[i1])/float(sumov[i1]));
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
            if ( ih.hkl_class().centric() ) {
                //no anomalous information
                sumov[bin] += eps;
                if ( !clipper::Util::is_nan(obs_pl(isig_ano[ih]) )  &&  !clipper::Util::is_nan(obs_mi(isig_ano[ih]) ) ) {
                    meanI[bin] += 0.5*(obs_pl(isig_ano[ih])+obs_mi(isig_ano[ih]));
                } else if ( !clipper::Util::is_nan(obs_pl(isig_ano[ih]) ) ) {
                    meanI[bin] += obs_pl(isig_ano[ih]);
                } else {
                    meanI[bin] += obs_mi(isig_ano[ih]);
                }
            } else if ( !clipper::Util::is_nan(obs_pl(isig_ano[ih]) )  &&  !clipper::Util::is_nan(obs_pl(isig_ano[ih]) ) ) {
                clipper::ftype Ip = obs_pl(isig_ano[ih]);
                clipper::ftype sp = sigobs_pl(isig_ano[ih]);
                clipper::ftype Im = obs_mi(isig_ano[ih]);
                clipper::ftype sm = sigobs_mi(isig_ano[ih]);
                clipper::ftype dI = Ip - Im;
                clipper::ftype ds = std::sqrt(sp*sp+sm*sm);
                if ( dI/ds >= 3.0 && Ip/sp >= 3.0 && Im/sm > 3.0 ) {
                    summeas[bin] += eps;
                }
                meandI[bin] += dI;
                meandIsigdI[bin] += dI/ds;
                meanI[bin] += 0.5*(obs_pl(isig_ano[ih])+obs_mi(isig_ano[ih]));
            } else  {
                sumov[bin] += eps;
                if ( !clipper::Util::is_nan(obs_pl(isig_ano[ih]) ) ) {
                    meanI[bin] += obs_pl(isig_ano[ih]);
                } else {
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
			printf("%10.6f %12.4e %12.4e %12.4e %12.4e\n",res,meandI[i],meandIsigdI[i],meandI[i]/meanI[i],float(summeas[i1])/float(sumov[i1]));
		}
		printf("$$\n\n");
    }

    
            template class AnomStats<clipper::ftype32>;
            template class AnomStats<clipper::ftype64>;
}