//
//     CTRUNCATE
//     Copyright (C) 2006-2011 Norman Stein, Charles Ballard
//
//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Application.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//

#include "ctruncate_aniso.h"
#include "ctruncate_utils.h"
#include "intensity_target.h"
#include "ctruncate_wilson.h"
#include "cpsf_utils.h"
#include "best.h"
#include "best_rna.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

#define ASSERT assert
#include <assert.h>

namespace ctruncate {

	void yorgo_modis_plot(clipper::HKL_data<clipper::data32::F_sigF>& fsig, float maxres, int nbins, CCP4Program& prog, clipper::U_aniso_orth uao )
	{
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		
        //get uao eigenvectors
		AnisoDirection direct(uao);
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
		printf("$$ 1/resol^2 Mn(F(1)) Mn(F(2)) Mn(F(3)) Mn(F) Mn(F/s(1)) Mn(F/s(2)) Mn(F/s(3)) Mn(F/s)");
		printf(" N(1) N(2) N(3) N completeness sig1 sig2 sig3$$\n$$\n");
		
		
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
		AnisoDirection direct(uao);
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
		printf("$$ 1/resol^2 Mn(I(1)) Mn(I(2)) Mn(I(3)) Mn(I) Mn(I/s(1)) Mn(I/s(2)) Mn(I/s(3)) Mn(I/s)");
		printf(" N(1) N(2) N(3) N completeness sig1 sig2 sig3$$\n$$\n");
		
		
		for(int i=0;i<nbins;i++){
			double res = maxres*(double(i)+0.5)/double(nbins);
			printf("%10.6f %12.4e %12.4e %12.4e %12.4e ",res,somdir[0][i],somdir[1][i],somdir[2][i],somov[i]);
			printf("%12.4e %12.4e %12.4e %12.4e ",somsddir[0][i],somsddir[1][i],somsddir[2][i],somsdov[i]);
			printf("%8d %8d %8d %8d",numdir[0][i],numdir[1][i],numdir[2][i],numov[i]);
			printf("%8.4f %8.4f %8.4f %8.4f\n",completeness[i],completeness1[i],completeness2[i],completeness3[i]);
		}
		printf("$$\n\n");
	}


	//******YorgoModis***************************************************************************
	
	/*template <class D> void YorgoModis<D>::operator() (clipper::HKL_data<D>& isig, clipper::Resolution reso )
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
				float epsiln = (_t == I) ? 1.0f/ih.hkl_class().epsilonc() : std::sqrt(1.0f/ih.hkl_class().epsilonc() );
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
	} */
	
	void YorgoModis::output() const
	{
		if (is_intensity() ) {
			printf("\n$TABLE: Intensity statistics:\n");
			printf("$GRAPHS");
			printf(": Mn(I) v resolution:N:1,2,3,4,5:\n");
			printf(": Mn(I/sd) v resolution:N:1,6,7,8,9:\n");
			printf(": No. reflections v resolution:N:1,10,11,12,13:\n");
			printf("$$ 1/resol^2 Mn(I(1)) Mn(I(2)) Mn(I(3)) Mn(I) Mn(I/s(1)) Mn(I/s(2)) Mn(I/s(3)) Mn(I/s)");
			printf("     N(1)     N(2)     N(3)     N$$\n$$\n");
		} else {
			printf("\n$TABLE: Anisotropy analysis (Yorgo Modis):\n");
			printf("$GRAPHS");
			printf(": Mn(F) v resolution:N:1,2,3,4,5:\n");
			printf(": Mn(F/sd) v resolution:N:1,6,7,8,9:\n");
			printf(": No. reflections v resolution:N:1,10,11,12,13:\n");
			printf("$$ 1/resol^2 Mn(F(1)) Mn(F(2)) Mn(F(3)) Mn(F) Mn(F/s(1)) Mn(F/s(2)) Mn(F/s(3)) Mn(F/s)");
			printf("     N(1)     N(2)     N(3)     N$$\n$$\n");
		}
		
		int bins=this->size();
		for(int i=0;i!=bins;++i){
			printf("   %6.4f %8.2f %8.2f %8.2f %8.2f ",_b_reso[i],_somdir[i],_somdir[bins+i],_somdir[2*bins+i],_somov[i]);
			printf("%8.2f   %8.2f   %8.2f   %8.2f ",_somsddir[i],_somsddir[bins+i],_somsddir[2*bins+i],_somsdov[i]);
			printf("%8.1f %8.1f %8.1f %8.1f\n",_numdir[i],_numdir[bins+i],_numdir[2*bins+i],_numov[i]);
		}
		printf("$$\n\n");
		printf("The directional plots are along the directions of the moments of the anisotropy temperature matrix.  These are ordered such that direction 1 has maximum alignment with a*, directions 2 with b*, etc.\n");
	}

	std::stringstream& YorgoModis::xml_output(std::stringstream& ss) const
	{
		std::string is = (is_intensity() ) ? "Intensities" : "Amplitude" ;
		ss << "<CCP4Table groupID=\"graph\" id=\"AnisotropyAnalysis\" title=\""<< is << " Anisotropy Analysis\">" << std::endl;
		ss << "<plot>" << std::endl;
		if (is_intensity() ) ss << "<title>Directional plots for mean(I)</title>" << std::endl;
		else ss << "<title>Directorional plots for mean(F)</title>" << std::endl;
        ss << "<xscale>oneoversqrt</xscale>" << std::endl;
		ss << "<yrange min=\"0\" max=\"None\"/>" << std::endl;
		ss << "<plotline xcol=\"1\" ycol=\"5\" >" <<  std::endl;
		ss << "<symbolsize>  0</symbolsize>" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>yellow</colour>" << std::endl;
		ss << "</plotline>" << std::endl;
		ss << "<plotline xcol=\"1\" ycol=\"4\" >" <<  std::endl;
		ss << "<symbolsize>  0</symbolsize>" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>red</colour>" << std::endl;
		ss << "</plotline>" << std::endl;
		ss << "<plotline xcol=\"1\" ycol=\"3\" >" << std::endl;
		ss << "<symbolsize>  0</symbolsize>" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>green</colour>" << std::endl;
		ss << "</plotline>" << std::endl;
		ss << "<plotline xcol=\"  1\" ycol=\"  2\" >" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>black</colour>" << std::endl;
		ss << "</plotline>" << std::endl;
		ss << "</plot>" << std::endl;
		ss << "<plot>" << std::endl;
		if (is_intensity() ) ss << "<title>Directorional plots for mean(I)/mean(Sigma)</title>" << std::endl;
		else ss << "<title>Directorional plots for mean(F)/mean(Sigma)</title>" << std::endl;
		ss << "<xscale>oneoversqrt</xscale>" << std::endl;
		ss << "<yrange min=\"0\" max=\"None\"/>" << std::endl;
		ss << "<plotline xcol=\"1\" ycol=\"9\" >" <<  std::endl;
		ss << "<symbolsize>  0</symbolsize>" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>yellow</colour>" << std::endl;
		ss << "</plotline>" << std::endl;		
		ss << "<plotline xcol=\"1\" ycol=\"8\" >" <<  std::endl;
		ss << "<symbolsize>  0</symbolsize>" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>red</colour>" << std::endl;
		ss << "</plotline>" << std::endl;
		ss << "<plotline xcol=\"1\" ycol=\"7\" >" << std::endl;
		ss << "<symbolsize>  0</symbolsize>" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>green</colour>" << std::endl;
		ss << "</plotline>" << std::endl;
		ss << "<plotline xcol=\"  1\" ycol=\"6\" >" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>black</colour>" << std::endl;
		ss << "</plotline>" << std::endl;
		ss << "</plot>" << std::endl;
		ss << "<plot>" << std::endl;
		ss << "<title>Directorional plots for number of reflections</title>" << std::endl;
        ss << "<xscale>oneoversqrt</xscale>" << std::endl;
		ss << "<yrange min=\"0\" max=\"None\"/>" << std::endl;
		ss << "<plotline xcol=\"1\" ycol=\"13\" >" <<  std::endl;
		ss << "<symbolsize>  0</symbolsize>" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>yellow</colour>" << std::endl;
		ss << "</plotline>" << std::endl;		
		ss << "<plotline xcol=\"1\" ycol=\"12\" >" <<  std::endl;
		ss << "<symbolsize>  0</symbolsize>" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>red</colour>" << std::endl;
		ss << "</plotline>" << std::endl;
		ss << "<plotline xcol=\"1\" ycol=\"11\" >" << std::endl;
		ss << "<symbolsize>  0</symbolsize>" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>green</colour>" << std::endl;
		ss << "</plotline>" << std::endl;
		ss << "<plotline xcol=\"  1\" ycol=\"10\" >" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>black</colour>" << std::endl;
		ss << "</plotline>" << std::endl;
		ss << "</plot>" << std::endl;
		if (this->is_intensity()) ss << "<headers separator=\" \">\n 1/resol^2 Mn(I(d1)) Mn(I(d2)) Mn(I(d3)) Mn(I(ov) Mn(I/sd(d1)) Mn(I/sd(d2)) Mn(I/sd(d3)) Mn(I/sd(ov))";
		else ss << "<headers separator=\" \">\n  1/resol^2 Mn(F(d1)) Mn(F(d2)) Mn(F(d3)) Mn(F) Mn(F/sd(d1)) Mn(F/sd(d2)) Mn(F/sd(d3)) Mn(F/sd)";
		ss << " N(d1) N(d2) N(d3) N \n </headers>" << std::endl;
		ss << "<data>" << std::endl;;
		int bins=this->size();
		for(int i=0;i!=bins;++i){
			ss << " " << std::fixed << std::setw(6) << std::setprecision(4)  << _b_reso[i];
			ss << std::fixed << std::setprecision(2) << " " << std::setw(8) << _somdir[i] << " " << std::setw(8) << _somdir[bins+i] << " " << std::setw(8) << _somdir[2*bins+i] << " " << std::setw(8) << _somov[i];
			ss << std::fixed << std::setprecision(2) << " " << std::setw(8) << _somsddir[i] << " " << std::setw(8) << _somsddir[bins+i] << " " << std::setw(8) << _somsddir[2*bins+i] << " " << std::setw(8) << _somsdov[i];
			ss << std::fixed << std::setprecision(1) << " " << std::setw(8) << _numdir[i] << " " << std::setw(8) << _numdir[bins+i] << " " << std::setw(8) << _numdir[2*bins+i] << " " << std::setw(8) << _numov[i] << std::endl;
		}
		ss << "</data>" << std::endl;
		ss << "</CCP4Table>" << std::endl;
        ss << "<AnisotropyAnalysis id='YorgoModis'>" << std::endl;
        ss << "<Comment id='YorgoModisPlot'>" << std::endl;
        ss << "The directional plots are along the directions of the moments of the anisotropy temperature matrix.  These are ordered such that d1 has maximum alignment with a*, d2 with b*, etc." << std::endl;
        ss << "</Comment>" << std::endl;
        ss << "</AnisotropyAnalysis>" << std::endl;
		return ss;
	}

	//******AnisoPlot***************************************************************************

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
	
	//---------Calculate anisotropy plot------------------------
	
	
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
	
	
    template <class SCALER, class DATA, class T> void AnisoCorr<SCALER,DATA,T>::calc(clipper::HKL_data<DATA>& isig)
    {
		const clipper::HKL_info& hklinf = isig.hkl_info();
		if (_range.range() < 0.0 ) {
			_range = isig.hkl_info().invresolsq_range();
		}
		
		if ( _is_protein || _is_nucl ) {
			clipper::Cell cell(hklinf.cell());
			clipper::Spacegroup spgr(hklinf.spacegroup());
			clipper::HKL_data<DATA> Ibest(hklinf);
			
			//generate reference scattering curve (initially protein only)
			if (_is_protein) {
				Scattering scat;
				clipper::ftype totalscatter = spgr.num_symops()*scat(cell, spgr);
				
				for ( clipper::HKL_data_base::HKL_reference_index ih = Ibest.first(); !ih.last(); ih.next() ) {
					T reso = ih.invresolsq();
					Ibest[ih] = DATA(ih.hkl_class().epsilon()*totalscatter*ctruncate::Best::value(reso), 1.0f);
				} // scale against BEST
			} else {
				Scattering scat(Scattering::NUCLEIC);
				clipper::ftype totalscatter = spgr.num_symops()*scat(cell, spgr);
				
				for ( clipper::HKL_data_base::HKL_reference_index ih = Ibest.first(); !ih.last(); ih.next() ) {
					T reso = ih.invresolsq();
                    Ibest[ih] = DATA(ih.hkl_class().epsilon()*totalscatter*ctruncate::Best_rna::value(reso), 1.0f);
				} // scale against BEST
				
			}
			_iscale( Ibest, isig );
		} else {
			_iscale( isig, _range, 12);
		}
		
        return ;
    }
    
    //! calculate anisotropy correction
    template <class SCALER, class DATA, class T> 
    const clipper::U_aniso_orth& AnisoCorr<SCALER,DATA,T>::operator()(clipper::HKL_data<DATA>& observed, bool protein, bool rna, clipper::Range<clipper::ftype> reso )
    {
		_is_protein = protein;
		_is_nucl = rna;
		_range = reso;
        calc(observed);
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
    
	AnisoDirection::AnisoDirection(const clipper::U_aniso_orth& uao) {
		this->operator()(uao);
		return;
	}
	
	const std::vector<clipper::ftype>& AnisoDirection::operator()(const clipper::U_aniso_orth& uao)
    {
        // Eigenvalue calculation
        _uao = &uao;
        clipper::Matrix<clipper::ftype> mat( 3, 3, 0.0 );
        for (int i=0; i !=3; ++i) {
            for (int j=0; j !=3 ; ++j) {
                mat(i,j) = uao(i,j);
            }
        }
        std::vector<clipper::ftype> v = mat.eigen( true );
        
        _max = -999.0;
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
            _eigenvectors.push_back(clipper::Vec3<clipper::ftype>(mat(0,close[i]),mat(1,close[i]),mat(2,close[i])) );
        }
        ASSERT( _eigenvalues.size() == 3 );
        return _eigenvalues;
    }
	
	//---------AnisoAnalysis--------------------------------------------
	bool AnisoAnalysis::allowed_by_symmetry(const clipper::Spacegroup& spgr, const clipper::Cell& cell)
	{
		clipper::U_aniso_orth uao_sum(1,3,5,7,11,13);
		for (int i = 1; i != spgr.num_symops() ; ++i ) {
			clipper::U_aniso_orth uao(1,3,5,7,11,13);
			clipper::U_aniso_orth tmp = uao.transform(spgr.symop(i).rtop_orth(cell) );
			uao_sum = uao_sum + tmp;
		}
		return (std::fabs(uao_sum.mat00() - uao_sum.mat11() ) > 0.5 || std::fabs(uao_sum.mat00() - uao_sum.mat22()) > 0.5 );
	}
	
	bool AnisoAnalysis::is_anisotropic() const {
        if ( !allowed_by_symmetry() ) return false;
		clipper::ftype vm = std::exp(_dm.max() );
		const std::vector<clipper::ftype> v = _dm.eigenValues();
		clipper::ftype v1(std::exp(v[0])/vm), v2(std::exp(v[1])/vm), v3(std::exp(v[2])/vm);
		clipper::ftype ratio = std::min(v1,std::min(v2,v3) );
		return ratio < 0.95;
	}
	
	void AnisoAnalysis::output() const
	{
		if (_anisobysymm) {
			
			std::vector<clipper::ftype> v = _dm.eigenValues();
            std::vector<clipper::Vec3<clipper::ftype> > ev = _dm.eigenVectors();
			clipper::ftype max = std::exp(_dm.max() );
			clipper::ftype v1(std::exp(v[0])/max), v2(std::exp(v[1])/max), v3(std::exp(v[2])/max);
			clipper::ftype ratio = std::min(v1,std::min(v2,v3) );
			clipper::U_aniso_frac uaf = _u_aniso_orth_i.u_aniso_frac( _base->hkl_info().cell() );
			
			printf("\nANISOTROPY ANALYSIS:\n");
            printf("\nAnalysis using data from %6.2fA to %6.2fA.\n",1.0/std::sqrt(_range.min()),1.0/std::sqrt(_range.max() ) );
			printf("\nEigenvalues: %8.4f %8.4f %8.4f\n", clipper::Util::u2b(v[0]/2.0),clipper::Util::u2b(v[1]/2.0),clipper::Util::u2b(v[2]/2.0) );
			printf("Eigenvalue ratios: %8.4f %8.4f %8.4f\n", v1, v2, v3);
            std::cout << std::endl;
			if (ratio < 0.5 ) std::cout << "Warning: data is highly anisotropic." << std::endl;
			else if (ratio < 0.9) std::cout << "Some anisotropy detetect.  This may have an effect on statistics." << std::endl;
			else std::cout << "Little or no anisotropy detected. " << std::endl;
            std::cout << "The presence of anisotropy may indicate that the crystal is poorly ordered along one of the axes." << std::endl;
			
            printf("\nAnisotropic B scaling (orthogonal coords):\n\n");
			
			printf("| %10.6f %10.6f %10.6f |\n",clipper::Util::u2b( _u_aniso_orth_f(0,0) ), clipper::Util::u2b( _u_aniso_orth_f(0,1) ), clipper::Util::u2b( _u_aniso_orth_f(0,2) ) );
			printf("| %10.6f %10.6f %10.6f |\n",clipper::Util::u2b( _u_aniso_orth_f(1,0) ), clipper::Util::u2b( _u_aniso_orth_f(1,1) ), clipper::Util::u2b( _u_aniso_orth_f(1,2) ) );
			printf("| %10.6f %10.6f %10.6f |\n",clipper::Util::u2b( _u_aniso_orth_f(2,0) ), clipper::Util::u2b( _u_aniso_orth_f(2,1) ), clipper::Util::u2b( _u_aniso_orth_f(2,2) ) );
            
			printf("\nAnisotropic U (orthogonal coords):\n\n");
			printf("| %10.6f %10.6f %10.6f |\n", _u_aniso_orth_i(0,0) ,  _u_aniso_orth_i(0,1) ,  _u_aniso_orth_i(0,2)  );
			printf("| %10.6f %10.6f %10.6f |\n", _u_aniso_orth_i(1,0) ,  _u_aniso_orth_i(1,1) ,  _u_aniso_orth_i(1,2)  );
			printf("| %10.6f %10.6f %10.6f |\n", _u_aniso_orth_i(2,0) ,  _u_aniso_orth_i(2,1) ,  _u_aniso_orth_i(2,2)  );
			
			printf("\nEigenvector breakdown:\n\n");
            printf("Eigenvalue  Eigenvector(a*,b*,c*)\n");
			
			printf("%10.6f ( %10.6f %10.6f %10.6f )\n", v[0], ev[0][0] ,  ev[0][1] ,  ev[0][2]  );
			printf("%10.6f ( %10.6f %10.6f %10.6f )\n", v[1], ev[1][0] ,  ev[1][1] ,  ev[1][2]  );
			printf("%10.6f ( %10.6f %10.6f %10.6f )\n", v[2], ev[2][0] ,  ev[2][1] ,  ev[2][2]  );
			
            printf("\nAnisotropic correction (orthogonal coords):\n\n");
			printf("| %10.6f %10.6f %10.6f |\n", _u_aniso_orth_corr_f(0,0) ,  _u_aniso_orth_corr_f(0,1) ,  _u_aniso_orth_corr_f(0,2)  );
			printf("| %10.6f %10.6f %10.6f |\n", _u_aniso_orth_corr_f(1,0) ,  _u_aniso_orth_corr_f(1,1) ,  _u_aniso_orth_corr_f(1,2)  );
			printf("| %10.6f %10.6f %10.6f |\n", _u_aniso_orth_corr_f(2,0) ,  _u_aniso_orth_corr_f(2,1) ,  _u_aniso_orth_corr_f(2,2)  );
            
            _ym.output();
		} else { 
			printf("\nNo anisotropy by symmetry. \n");
		}
		
        std::cout << std::endl;
        
		return;
	}
	
	std::stringstream& AnisoAnalysis::xml_output(std::stringstream& ss) const
	{
		ss << "<AnisotropyAnalysis>" << std::endl;
		if (_anisobysymm) {
            std::vector<clipper::ftype> v = _dm.eigenValues();
            std::vector<clipper::Vec3<clipper::ftype> > ev = _dm.eigenVectors();
            clipper::U_aniso_frac uai = _u_aniso_orth_i.u_aniso_frac( _base->hkl_info().cell() );
            clipper::U_aniso_frac uaf = _u_aniso_orth_f.u_aniso_frac( _base->hkl_info().cell() );
            clipper::U_aniso_frac uac = _u_aniso_orth_corr_f.u_aniso_frac( _base->hkl_info().cell() );
            clipper::ftype max = std::exp(_dm.max() );
            clipper::ftype v1(std::exp(v[0])/max), v2(std::exp(v[1])/max), v3(std::exp(v[2])/max);
            clipper::ftype ratio = std::min(v1,std::min(v2,v3) );
            
			ss << "  <UAniso id=\"obs\" type=\"intensity\">" << std::endl;
			ss << "    <Orthogonal>" <<std::endl;
			ss << std::fixed << std::setw(10) << std::setprecision(6) << _u_aniso_orth_i(0,0) << " " << std::setw(10) << _u_aniso_orth_i(0,1) << " " << std::setw(10) << _u_aniso_orth_i(0,2) << std::endl;
			ss << std::fixed << std::setw(10) << std::setprecision(6) << _u_aniso_orth_i(1,0) << " " << std::setw(10) << _u_aniso_orth_i(1,1) << " " << std::setw(10) << _u_aniso_orth_i(1,2) << std::endl;
			ss << std::fixed << std::setw(10) << std::setprecision(6) << _u_aniso_orth_i(2,0) << " " << std::setw(10) << _u_aniso_orth_i(2,1) << " " << std::setw(10) << _u_aniso_orth_i(2,2) << std::endl;
			ss << "    </Orthogonal>" << std::endl;
			
			ss << "    <Fractional>" <<std::endl;
			ss << std::fixed << std::setw(10) << std::setprecision(6) << uai(0,0) << " " << uai(0,1) << " " << uai(0,2) << std::endl;
			ss << std::fixed << std::setw(10) << std::setprecision(6) << uai(1,0) << " " << uai(1,1) << " " << uai(1,2) << std::endl;
			ss << std::fixed << std::setw(10) << std::setprecision(6) << uai(2,0) << " " << uai(2,1) << " " << uai(2,2) << std::endl;
			ss << "    </Fractional>" << std::endl;
            ss << "  </UAniso>" << std::endl;
			ss << "   <Eigenvalue id=\"0\">" << std::fixed << std::setw(8) << std::setprecision(4) << v[0] << "</Eigenvalue>" << std::endl;
            ss << "  <EigenVector id=\"0\">" << std::fixed << std::setw(10) << std::setprecision(6) << ev[0][0] << " " << ev[0][1] << " " << ev[0][2] << "</EigenVector>" << std::endl;
            ss << "   <Eigenvalue id=\"1\">" << std::fixed << std::setw(8) << std::setprecision(4) << v[1] << "</Eigenvalue>" << std::endl;
            ss << "  <EigenVector id=\"1\">" <<std::fixed << std::setw(10) << std::setprecision(6) << ev[1][0] << " " << ev[1][1] << " " << ev[1][2] << "</EigenVector>" << std::endl;
            ss << "   <Eigenvalue id=\"2\">" << std::fixed << std::setw(8) << std::setprecision(4) << v[2] << "</Eigenvalue>" << std::endl;
            ss << "  <EigenVector id=\"2\">" <<std::fixed << std::setw(10) << std::setprecision(6) << ev[2][0] << " " << ev[2][1] << " " << ev[2][2] << "</EigenVector>" << std::endl;
            ss << "  <BAniso id=\"obs\" type=\"amplitude\">" << std::endl;
			ss << "    <Orthogonal>" <<std::endl;
			ss << std::fixed << std::setw(10) << std::setprecision(6) << clipper::Util::u2b(_u_aniso_orth_f(0,0) ) << " " << clipper::Util::u2b(_u_aniso_orth_f(0,1) ) << " " << clipper::Util::u2b(_u_aniso_orth_f(0,2) ) << std::endl;
			ss << std::fixed << std::setw(10) << std::setprecision(6) << clipper::Util::u2b(_u_aniso_orth_f(1,0) ) << " " << clipper::Util::u2b(_u_aniso_orth_f(1,1) ) << " " << clipper::Util::u2b(_u_aniso_orth_f(1,2) ) << std::endl;
			ss << std::fixed << std::setw(10) << std::setprecision(6) << clipper::Util::u2b(_u_aniso_orth_f(2,0) ) << " " << clipper::Util::u2b(_u_aniso_orth_f(2,1) ) << " " << clipper::Util::u2b(_u_aniso_orth_f(2,2) ) << std::endl;
			ss << "    </Orthogonal>" << std::endl;
			ss << "    <Fractional>" <<std::endl;
			ss << std::fixed << std::setw(10) << std::setprecision(6) << clipper::Util::u2b(uaf(0,0) ) << " " << clipper::Util::u2b(uaf(0,1) ) << " " << clipper::Util::u2b(uaf(0,2) ) << std::endl;
			ss << std::fixed << std::setw(10) << std::setprecision(6) << clipper::Util::u2b(uaf(1,0) ) << " " << clipper::Util::u2b(uaf(1,1) ) << " " << clipper::Util::u2b(uaf(1,2) ) << std::endl;
			ss << std::fixed << std::setw(10) << std::setprecision(6) << clipper::Util::u2b(uaf(2,0) ) << " " << clipper::Util::u2b(uaf(2,1) ) << " " << clipper::Util::u2b(uaf(2,2) ) << std::endl;
			ss << "    </Fractional>" << std::endl;
            ss << "  </BAniso>" << std::endl;
            ss << "  <UAniso id=\"correction\" type=\"amplitude\">" << std::endl;
			ss << "    <Orthogonal>" <<std::endl;
			ss << std::fixed << std::setw(10) << std::setprecision(6) << _u_aniso_orth_corr_f(0,0) << " " << _u_aniso_orth_corr_f(0,1) << " " << _u_aniso_orth_corr_f(0,2) << std::endl;
			ss << std::fixed << std::setw(10) << std::setprecision(6) << _u_aniso_orth_corr_f(1,0) << " " << _u_aniso_orth_corr_f(1,1) << " " << _u_aniso_orth_corr_f(1,2) << std::endl;
			ss << std::fixed << std::setw(10) << std::setprecision(6) << _u_aniso_orth_corr_f(2,0) << " " << _u_aniso_orth_corr_f(2,1) << " " << _u_aniso_orth_corr_f(2,2) << std::endl;
			ss << "    </Orthogonal>" << std::endl;
			ss << "    <Fractional>" <<std::endl;
			ss << std::fixed << std::setw(10) << std::setprecision(6) << uac(0,0) << " " <<  uac(0,1) << " " << uac(0,2) << std::endl;
			ss << std::fixed << std::setw(10) << std::setprecision(6) << uac(1,0) << " " <<  uac(1,1) << " " << uac(1,2) << std::endl;
			ss << std::fixed << std::setw(10) << std::setprecision(6) << uac(2,0) << " " <<  uac(2,1) << " " << uac(2,2) << std::endl;
			ss << "    </Fractional>" << std::endl;
            ss << "  </UAniso>" << std::endl;
			ss << "  <Comment id=\"AnisotropyResult\">" << std::endl;
			if (ratio < 0.5 ) ss << "Warning: data is highly anisotropic." << std::endl;
			else if (ratio < 0.9) ss << "Some anisotropy detetect.  This may have an effect on statistics." << std::endl;
			else ss << "Little or no anisotropy detected. " << std::endl;
            ss << "The presence of anisotropy may indicate that the crystal is poorly ordered along one of the axes." << std::endl;
			ss << "  </Comment>" << std::endl;
            ss << "</AnisotropyAnalysis>" << std::endl;
			_ym.xml_output(ss); 
		} else {
			ss << "  <Comment id=\"Anisotropy Analysis\">" << std::endl;
			ss << "No anisotropy by symmetry." << std::endl;
			ss << "  </Comment>" << std::endl;
            ss << "</AnisotropyAnalysis>" << std::endl;
		}
		return ss;
	}
	
    //******HKLStats_d_completeness***************************************************************
	
	HKLStats_d_completeness::HKLStats_d_completeness(const clipper::HKL_data_base& hkldata, clipper::U_aniso_orth& uao, clipper::ftype val)
	{
		_val=val;
        //get uao eigenvectors
		AnisoDirection direct(uao);
		for (int i = 0 ; i !=3 ; ++i)
			for (int j = 0 ; j != 3 ; ++j)
                _e123(i,j) = (direct.eigenVectors())[i][j];
		init(hkldata);
		calc(hkldata);
	}
	
	HKLStats_d_completeness::HKLStats_d_completeness(const ResolStats_base& base, clipper::U_aniso_orth& uao, clipper::ftype val)
	{
		_val=val;//get uao eigenvectors
		AnisoDirection direct(uao);
		for (int i = 0 ; i !=3 ; ++i)
			for (int j = 0 ; j != 3 ; ++j)
                _e123(i,j) = (direct.eigenVectors())[i][j];
		init(base);
		calc(*(this->parent() ) );
	}
	
	void HKLStats_d_completeness::calc(const clipper::HKL_data_base& hkldata)
	{
		int nbins = this->size();
        int nfriedal(0);
		
        clipper::ftype cone(30.0),ang, cosan; //hardwired for now
        
		_d1_completeness.resize(nbins,0.0);
        _d2_completeness.resize(nbins,0.0);
        _d3_completeness.resize(nbins,0.0);
		
		std::vector<clipper::ftype> d1_sumov(nbins,0.0);
        std::vector<clipper::ftype> d2_sumov(nbins,0.0);
        std::vector<clipper::ftype> d3_sumov(nbins,0.0);
        
        clipper::Spacegroup spg = hkldata.hkl_info().spacegroup();
        clipper::Cell cell = hkldata.hkl_info().cell();
		
		clipper::xtype working[hkldata.data_size()];
		for ( clipper::HKL_data_base::HKL_reference_index ih = hkldata.first(); !ih.last(); ih.next() ) {
			clipper::ftype eps = (this->is_intensity() ) ? 1.0/ih.hkl_class().epsilonc() : 1.0/std::sqrt(ih.hkl_class().epsilonc());
            clipper::HKL ri = ih.hkl();
            nfriedal = (ih.hkl_class().centric() ) ? 1 : 2;
            int bin = this->ResolStats_base::operator()(ih.invresolsq() );
            if (!hkldata.missing(ih.index() ) ) {
                clipper::ftype I(0.0);
                clipper::ftype sig(1.0);
                clipper::ftype s = ih.invresolsq();
                hkldata.data_export(ih.hkl(),working);
                if (this->is_anomalous() ) {
                    clipper::ftype Ip(working[0]);
                    clipper::ftype sp(working[1]);
                    clipper::ftype Im(working[2]);
                    clipper::ftype sm(working[3]);
                    if (!clipper::Util::is_nan(Ip) && !clipper::Util::is_nan(Im) && sp > 0.0 && sm > 0.0 ) {
                        for ( int jsym = 0; jsym != spg.num_primitive_symops() ; ++jsym ) {
                            clipper::HKL rj = ri.transform( spg.primitive_symop( jsym ) );
                            clipper::Vec3<clipper::ftype> hc = _e123*clipper::Vec3<clipper::ftype>(rj.coord_reci_orth(cell) );  //transpose into eigenspace
                            if ( acos(std::min((fabs( hc[0] )/sqrt(ih.invresolsq())),1.0)) < clipper::Util::d2rad(cone) ) {
                                if ( Ip/sp >= _val ) _d1_completeness[bin] += eps;
                                if ( Im/sm >= _val ) _d1_completeness[bin] += eps;
                            }
                            if ( acos(std::min((fabs( hc[1] )/sqrt(ih.invresolsq())),1.0)) < clipper::Util::d2rad(cone) ) {
                                if ( Ip/sp >= _val ) _d2_completeness[bin] += eps;
                                if ( Im/sm >= _val ) _d2_completeness[bin] += eps;
                            }
                            if ( acos(std::min((fabs( hc[2] )/sqrt(ih.invresolsq())),1.0)) < clipper::Util::d2rad(cone) ) {
                                if ( Ip/sp >= _val ) _d3_completeness[bin] += eps;
                                if ( Im/sm >= _val ) _d3_completeness[bin] += eps;
                            }
                        }
                    } else if (!clipper::Util::is_nan(Ip) && sp > 0.0 ) {                            for ( int jsym = 0; jsym != spg.num_primitive_symops() ; ++jsym ) {
                            clipper::HKL rj = ri.transform( spg.primitive_symop( jsym ) );
                            clipper::Vec3<clipper::ftype> hc = _e123*clipper::Vec3<clipper::ftype>(rj.coord_reci_orth(cell) );  //transpose into eigenspace
                            if ( acos(std::min((fabs( hc[0] )/sqrt(ih.invresolsq())),1.0)) < clipper::Util::d2rad(cone) ) {
                                if ( Ip/sp >= _val ) _d1_completeness[bin] += eps;
                            }
                            if ( acos(std::min((fabs( hc[1] )/sqrt(ih.invresolsq())),1.0)) < clipper::Util::d2rad(cone) ) {
                                if ( Ip/sp >= _val ) _d2_completeness[bin] += eps;
                            }
                            if ( acos(std::min((fabs( hc[2] )/sqrt(ih.invresolsq())),1.0)) < clipper::Util::d2rad(cone) ) {
                                if ( Ip/sp >= _val ) _d3_completeness[bin] += eps;
                            }
                        }
                    } else if (!clipper::Util::is_nan(Im) && sm > 0.0 ) {
                        for ( int jsym = 0; jsym != spg.num_primitive_symops() ; ++jsym ) {
                            clipper::HKL rj = ri.transform( spg.primitive_symop( jsym ) );
                            clipper::Vec3<clipper::ftype> hc = _e123*clipper::Vec3<clipper::ftype>(rj.coord_reci_orth(cell) );  //transpose into eigenspace
                            if ( acos(std::min((fabs( hc[0] )/sqrt(ih.invresolsq())),1.0)) < clipper::Util::d2rad(cone) ) {
                                if ( Im/sm >= _val ) _d1_completeness[bin] += eps;
                            }
                            if ( acos(std::min((fabs( hc[1] )/sqrt(ih.invresolsq())),1.0)) < clipper::Util::d2rad(cone) ) {
                                if ( Im/sm >= _val ) _d2_completeness[bin] += eps;
                            }
                            if ( acos(std::min((fabs( hc[2] )/sqrt(ih.invresolsq())),1.0)) < clipper::Util::d2rad(cone) ) {
                                if ( Im/sm >= _val ) _d3_completeness[bin] += eps;
                            }
                        }
                    }
                } else {
                    clipper::ftype Ip(working[0]);
                    clipper::ftype sp(working[1]);
                    if (!clipper::Util::is_nan(Ip) && sp > 0.0 ) {
                            for ( int jsym = 0; jsym != spg.num_primitive_symops() ; ++jsym ) {
                            clipper::HKL rj = ri.transform( spg.primitive_symop( jsym ) );
                            clipper::Vec3<clipper::ftype> hc = _e123*clipper::Vec3<clipper::ftype>(rj.coord_reci_orth(cell) );  //transpose into eigenspace
                            if ( acos(std::min((fabs( hc[0] )/sqrt(ih.invresolsq())),1.0)) < clipper::Util::d2rad(cone) ) {
                                if ( Ip/sp >= _val ) _d1_completeness[bin] += eps;
                            }
                            if ( acos(std::min((fabs( hc[1] )/sqrt(ih.invresolsq())),1.0)) < clipper::Util::d2rad(cone) ) {
                                if ( Ip/sp >= _val ) _d2_completeness[bin] += eps;
                            }
                            if ( acos(std::min((fabs( hc[2] )/sqrt(ih.invresolsq())),1.0)) < clipper::Util::d2rad(cone) ) {
                                if ( Ip/sp >= _val ) _d3_completeness[bin] += eps;
                            }
                        }
                    }
                }
            }
            for ( int jsym = 0; jsym != spg.num_primitive_symops() ; ++jsym ) {
                clipper::HKL rj = ri.transform( spg.primitive_symop( jsym ) );
                clipper::Vec3<clipper::ftype> hc = _e123*clipper::Vec3<clipper::ftype>(rj.coord_reci_orth(cell) );  //transpose into eigenspace
                if ( acos(std::min((fabs( hc[0] )/sqrt(ih.invresolsq())),1.0)) < clipper::Util::d2rad(cone) ) {
                    d1_sumov[bin] += spg.num_primitive_symops()*float(nfriedal)*eps;
                }
                if ( acos(std::min((fabs( hc[1] )/sqrt(ih.invresolsq())),1.0)) < clipper::Util::d2rad(cone) ) {
                    d2_sumov[bin] += spg.num_primitive_symops()*float(nfriedal)*eps;
                }
                if ( acos(std::min((fabs( hc[2] )/sqrt(ih.invresolsq())),1.0)) < clipper::Util::d2rad(cone) ) {
                    d3_sumov[bin] += spg.num_primitive_symops()*float(nfriedal)*eps;
                }
            }
         }
		
		for (int i=0 ; i != nbins ; ++i) _d1_completeness[i] /= d1_sumov[i];
        for (int i=0 ; i != nbins ; ++i) _d2_completeness[i] /= d2_sumov[i];
        for (int i=0 ; i != nbins ; ++i) _d3_completeness[i] /= d3_sumov[i];
	}
	
	clipper::ftype HKLStats_d_completeness::operator[](const int index) const
	{
		int nbins = this->size();
		int bin = clipper::Util::bound( 0, index, 2 );
		return (_d1_completeness[bin]+_d2_completeness[bin]+_d3_completeness[bin])/3.0;
	}
    
    clipper::ftype HKLStats_d_completeness::operator()(const int dir, const int index) const
	{
		int nbins = this->size();
		int bin = clipper::Util::bound( 0, index, nbins-1 );
        clipper::ftype comp;
        switch (dir) {
            case 1:
                comp = _d1_completeness[bin];
            case 2:
                comp = _d2_completeness[bin];
            case 3:
                comp = _d3_completeness[bin];
            default:
                comp = 0.0;
        }
		return comp;
	}
	
	HKLStats_d_completeness& HKLStats_d_completeness::operator=(const HKLStats_d_completeness& orig)
	{
		init(orig);
		_val = orig._val;
		_d1_completeness = orig._d1_completeness;
        _d2_completeness = orig._d2_completeness;
        _d3_completeness = orig._d3_completeness;
        _e123 = orig._e123;
		return *this;
	}
    
    clipper::Vec3<clipper::ftype> HKLStats_d_completeness::direction(const int index) const
	{
        int bin = clipper::Util::bound( 0, index, 2 );
		return clipper::Vec3<clipper::ftype>(_e123(bin,0),_e123(bin,1),_e123(bin,2));
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
	
 
}

