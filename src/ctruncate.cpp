//
//     CTRUNCATE
//     Copyright (C) 2006-2008 Norman Stein
//
//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Application.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//

#include "clipper/clipper.h"
#include "clipper/clipper-contrib.h"
#include "clipper/clipper-ccp4.h"
#include "clipper/clipper-mmdb.h"
#include "clipper/core/clipper_util.h"
#include "clipper/core/atomsf.h"
#include "clipper/core/coords.h"
#include "ccp4_general.h"
#include "clipper/ccp4/ccp4_mtz_io.h"
#include "cmtzlib.h"
#include "csymlib.h"
#include <iostream>
#include <math.h>
#include "intensity_target.h"  // contains additions to resol_targetfn.h

#include "cpsf_utils.h"


using namespace clipper;

// replacement for Wilson/Truncate

int truncate(  HKL_data<data32::I_sigI> isig,   HKL_data<data32::I_sigI>& jsig,   HKL_data<data32::F_sigF>& fsig,
			   clipper::ResolutionFn Sigma, float scalef, CSym::CCP4SPG *spg1);
int truncate_acentric(float I, float sigma, float S, float& J, float& sigJ, float& F, float& sigF);
int truncate_centric(float I, float sigma, float S, float& J, float& sigJ, float& F, float& sigF);
void straight_line_fit(std::vector<float> x, std::vector<float> y, std::vector<float> w, int n, float &a, float &b, float &siga, float &sigb);
void tricart(Cell cell, Mat33<float>& transf);
void Htest( HKL_data<data32::I_sigI> isig, Mat33<int> twinop, int &itwin, bool debug );

int main(int argc, char **argv)
{
  CCP4Program prog( "ctruncate", "0.1.02", "$Date: 2008/01/07" );
  
  // defaults
  clipper::String outfile = "ctruncate_out.mtz";
  clipper::String outcol = "F";
  clipper::String meancol = "IMEAN";
  clipper::String pluscol = "I(+)";
  clipper::String minuscol = "I(-)";
  clipper::String ipfile = "NONE";
  bool aniso = true;
  bool debug = false;

  int mtzinarg = 0;
  int anomalous = 0;

  // clipper seems to use its own column labels, then append yours

  CCP4MTZfile mtzfile, mtzout;
  HKL_info hklinf, hklp;

  // command input
  printf("\nUSER SUPPLIED INPUT:\n");
  CCP4CommandInput args( argc, argv, true );  
  int arg = 0;
  while ( ++arg < args.size() ) {
    if ( args[arg] == "-mtzin" || args[arg] == "-hklin") {
		if ( ++arg < args.size() ) {
			ipfile = args[arg];
			mtzinarg = arg;
		}
    } else if ( args[arg] == "-mtzout" || args[arg] == "-hklout") {
      if ( ++arg < args.size() ) outfile = args[arg];
    } else if ( args[arg] == "-colin" ) {
      if ( ++arg < args.size() ) meancol = args[arg];
    } else if ( args[arg] == "-colplus" ) {
      if ( ++arg < args.size() ) pluscol = args[arg];
	  anomalous = 1;
    } else if ( args[arg] == "-colminus" ) {
      if ( ++arg < args.size() ) minuscol = args[arg];
	  anomalous = 1;
    } else if ( args[arg] == "-colout" ) {
      if ( ++arg < args.size() ) outcol = args[arg];
    } else if ( args[arg] == "-no-aniso" ) {
      aniso = false;
	} else if ( args[arg] == "-debug" ) {
      debug = true;
	} else {
	  printf("Unrecognised argument\n");
	  return(0);
	}

  }
  if ( args.size() <= 1 ) {
	  CCP4::ccperror(1,"Usage: ctruncate -mtzin <filename>  -mtzout <filename>  -colin <colpath> -colplus <colpath> -colminus <colpath>");
  }

  if (mtzinarg == 0) CCP4::ccperror(1, "No input mtz file");

  typedef clipper::HKL_data_base::HKL_reference_index HRI;

  //mtzfile.open_read( argv[1] );
  mtzfile.open_read( ipfile );
  mtzfile.import_hkl_info( hklinf );  
  // allocate memory to isig by reading in hklinf before declaring isig
  HKL_data<data32::I_sigI> isig(hklinf);   // raw I and sigma
  HKL_data<data32::I_sigI> jsig(hklinf);   // post-truncate I and sigma
  HKL_data<data32::F_sigF> fsig(hklinf);   // post-truncate F and sigma 
  HKL_data<data32::I_sigI> isig_plus(hklinf);   // raw I(+) and sigma
  HKL_data<data32::I_sigI> jsig_plus(hklinf);   // post-truncate I and sigma
  HKL_data<data32::F_sigF> fsig_plus(hklinf);   // post-truncate F and sigma 
  HKL_data<data32::I_sigI> isig_minus(hklinf);   // raw I(-) and sigma
  HKL_data<data32::I_sigI> jsig_minus(hklinf);   // post-truncate I and sigma
  HKL_data<data32::F_sigF> fsig_minus(hklinf);   // post-truncate F and sigma 
  HKL_data<data32::F_sigF> Dano(hklinf);   // anomalous difference and sigma 
  HKL_data<data32::I_sigI> ianiso(hklinf);   // anisotropy corrected I and sigma
  // column labels originally hard wired  (aucn_mrg.mtz from $CEXAM used for Testing) 
  //mtzfile.import_hkl_data( isig, "/*/*/[I(+),SIGI(+)]" );
  meancol = "/*/*/[" + meancol + ",SIG" + meancol + "]";
  //mtzfile.import_hkl_data( isig, "/*/*/[IMEAN,SIGIMEAN]" );

//  clipper::MTZcrystal cxtl;
//  mtzfile.import_crystal( cxtl, meancol );
//  clipper::HKL_data<clipper::data32::F_sigF> faniso( hklinf, cxtl );  // don't seem to need crystal info
  clipper::HKL_data<clipper::data32::F_sigF> faniso( hklinf );

  mtzfile.import_hkl_data( isig, meancol );

  if (anomalous) {
      pluscol = "/*/*/[" + pluscol + ",SIG" + pluscol + "]";
	  mtzfile.import_hkl_data( isig_plus, pluscol );
      minuscol = "/*/*/[" + minuscol + ",SIG" + minuscol + "]";
	  mtzfile.import_hkl_data( isig_minus, minuscol );
  }

  //mtzfile.import_hkl_data( isig, "/*/*/[I,SIGI]" );
  // need this mumbo-jumbo in order to write to output file
  if ( outcol[0] != '/' ) outcol = mtzfile.assigned_paths()[0].notail()+"/"+outcol;

  mtzfile.close_read();



  // check for pseudo translation (taken from cpatterson)
  clipper::Spacegroup spgr = mtzfile.spacegroup();
  clipper::Cell       cell1 = mtzfile.cell();
  clipper::Resolution reso;
  clipper::Grid_sampling grid;
  clipper::String opfile = "patterson.map";
  reso = mtzfile.resolution();

  // get Patterson spacegroup
  clipper::Spacegroup
    pspgr( clipper::Spgr_descr( spgr.generator_ops().patterson_ops() ) );
  hklp.init( pspgr, cell1, reso, true );

  // make patterson coeffs
  clipper::HKL_data<clipper::data32::F_phi> fphi( hklp );
  for ( HRI ih = fphi.first(); !ih.last(); ih.next() ) {
    clipper::data32::I_sigI i = isig[ih.hkl()];
    if ( !i.missing() ) {
      fphi[ih].f() = i.I();
      fphi[ih].phi() = 0.0 ;
    }
  }

  // make grid if necessary
  if ( grid.is_null() ) grid.init( pspgr, cell1, reso );

  // make xmap
  clipper::Xmap<float> patterson( pspgr, cell1, grid );
  patterson.fft_from( fphi );

  // write map
  clipper::CCP4MAPfile mapout;
  mapout.open_write( opfile );
  mapout.export_xmap( patterson );
  mapout.close_write();


  // use Charles's stuff to find peaks
  PeakSearch pksch;                      // peak search object
  PeakInterp pkinterp;                   // peak interpolation methods

  //std::vector<clipper::Vec3<float> > cvec;
  		//harker hs(spgr);
		//clipper::Resolution p_reso( sreso.limit() );
		//clipper::Grid_sampling tgrid( pattspg, cell, p_reso);
		//clipper::Xmap<float> patterson( pattspg, cell, tgrid );
		//clipper::HKL_info hklp( pattspg, cell, reso, true);
  //std::vector<clipper::Coord_frac> origins = alt_origin(spgr);

  int npeak = 5;

  std::vector<int> ppks = pksch( patterson );

  float top_peak = patterson.get_data( ppks[0] );
  float next_peak = patterson.get_data( ppks[1] );
  clipper::Coord_frac c0 = patterson.coord_of( ppks[1] ).coord_frac(grid);
  float ratio = next_peak/top_peak;
  float dist2 = pow(c0[0], 2.0) + pow(c0[1], 2.0) + pow(c0[2], 2.0);
  // look for peaks > 20% of origin peak and at least 0.1 distant from origin
  printf("\n\nTRANSLATIONAL NCS:\n\n");
  if ( debug || (ratio > 0.2 && dist2 > 0.01) ) { 
      printf("\n\nNCS:\n\n");
	  printf("Translational NCS has been detected\n");
      printf("Ratio = %f\n",ratio);
      printf("Vector = (%6.3f, %6.3f, %6.3f)\n",c0[0],c0[1],c0[2]);
  }
  else {
	  printf("No translational NCS detected\n");
  }


/*  		clipper::Map_stats pattstats( patterson );
		float mean = pattstats.mean();
	    float std_dev = pattstats.std_dev();
		int i = 0;
		while ( (patterson.get_data( ppks[i] )-mean)/std_dev > 1.0f && cvec.size() < 4*npeak ) {
			clipper::Coord_frac c0 = patterson.coord_of( ppks[i] ).coord_frac(grid); */
			/* must screen for peaks on harker sections */
			// interpolate if requested
/*#ifdef INTERP
			clipper::Coord_map tmp = pkinterp( patterson , c0.coord_map(grid) );
			c0 = tmp.coord_frac(grid);
#endif	
			if ( c0.lengthsq(cell1) > std::pow(reso.limit()/2.0,2.0) && !hs.is_harker(c0) ) {


				bool ih = false;
				for ( int iih = 1 ; iih != origins.size() ; ++iih ) {
					clipper::Coord_frac t(c0 + origins[iih] );
					if ( hs.is_harker( t ) ) ih = true;
				}
				if (!ih) cvec.push_back(clipper::Vec3<float>(c0[0],c0[1],c0[2]) );
			}
			++i;
		}
		for (int i = 0 ; i != cvec.size() ; ++i ) {
			clipper::Coord_frac t(cvec[i][0],cvec[i][1],cvec[i][2]);
			std::cout << i << t.format() << patterson.interp<clipper::Interp_cubic>( t ) << std::endl;
		}*/



  // anisotropy correction
  if (aniso)  {
	  printf("\n\nANISOTROPY CORRECTION:\n");

      //clipper::HKL_data<clipper::data32::F_sigF> faniso( hklinf );
	  double Itotal = 0.0;
      for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
		  if ( !isig[ih].missing() ) {
	          double I = isig[ih].I();
	          double sigI = isig[ih].sigI();
			  if ( I > 0.0 ) {
			      Itotal += I;
	              faniso[ih] = clipper::data32::F_sigF( sqrt(I), 0.5*sigI/sqrt(I) );
			  }
          }
	  }

      // scale structure factors
      clipper::SFscale_aniso<float> sfscl( 3.0 );
      sfscl( faniso );
      //std::cout << "\nAnisotropic scaling:\n" << sfscl.u_aniso_orth().format() << "\n";
	  printf("\nAnisotropic scaling (orthogonal coords):\n\n");

	  printf("|%8.4f %8.4f %8.4f |\n", sfscl.u_aniso_orth()(0,0), sfscl.u_aniso_orth()(0,1), sfscl.u_aniso_orth()(0,2) );
	  printf("|%8.4f %8.4f %8.4f |\n", sfscl.u_aniso_orth()(1,0), sfscl.u_aniso_orth()(1,1), sfscl.u_aniso_orth()(1,2) );
	  printf("|%8.4f %8.4f %8.4f |\n", sfscl.u_aniso_orth()(2,0), sfscl.u_aniso_orth()(2,1), sfscl.u_aniso_orth()(2,2) );

	  clipper::U_aniso_orth uao = sfscl.u_aniso_orth();
	  clipper::U_aniso_frac uaf = uao.u_aniso_frac( cell1 );

      printf("\nAnisotropic U scaling (fractional coords):\n\n"); 

	  printf("| %11.3e %11.3e %11.3e |\n", uaf(0,0) ,  uaf(0,1) ,  uaf(0,2)  );
	  printf("| %11.3e %11.3e %11.3e |\n", uaf(1,0) ,  uaf(1,1) ,  uaf(1,2)  );
	  printf("| %11.3e %11.3e %11.3e |\n", uaf(2,0) ,  uaf(2,1) ,  uaf(2,2)  );

      printf("\nAnisotropic B scaling (fractional coords):\n\n"); 

	  printf("| %11.3e %11.3e %11.3e |\n",clipper::Util::u2b( uaf(0,0) ), clipper::Util::u2b( uaf(0,1) ), clipper::Util::u2b( uaf(0,2) ) );
	  printf("| %11.3e %11.3e %11.3e |\n",clipper::Util::u2b( uaf(1,0) ), clipper::Util::u2b( uaf(1,1) ), clipper::Util::u2b( uaf(1,2) ) );
	  printf("| %11.3e %11.3e %11.3e |\n",clipper::Util::u2b( uaf(2,0) ), clipper::Util::u2b( uaf(2,1) ), clipper::Util::u2b( uaf(2,2) ) );

	  double FFtotal = 0.0;
      for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
		  if ( !isig[ih].missing() ) {
			  double I = isig[ih].I();
			  double sigI = isig[ih].sigI();
	          double F = faniso[ih].f();
	          double sigF = faniso[ih].sigf();
			  if ( I > 0.0 ) {
		          FFtotal += F*F;
	              ianiso[ih] = clipper::data32::I_sigI( F*F, 2.0*F*sigF );
			  }
			  else {
                  ianiso[ih] = clipper::data32::I_sigI( I, sigI );
			  }
		  }
      }
	  double scalefac = Itotal/FFtotal;
	  if (debug) printf("\nscalefactor = %6.3f %6.3f %6.3f\n",scalefac,Itotal,FFtotal);
	  for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
		  if ( !isig[ih].missing() ) {
		      ianiso[ih].I() *= scalefac;
		      ianiso[ih].sigI() *=scalefac;
		  }
	  }
  }
  else {  // copy isig to aniso - is there a slicker way?
	  for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
	  	  double I = isig[ih].I();
	      double sigI = isig[ih].sigI();
     	  ianiso[ih] = clipper::data32::I_sigI( I, sigI );
	  }
  }



  // calculate Sigma (mean intensity in resolution shell)
  const int nprm = 12;

  std::vector<double> params_init( nprm, 1.0 );
  clipper::BasisFn_spline basis_fo( ianiso, nprm );
  TargetFn_meanInth<clipper::data32::I_sigI> target_fo( ianiso, 1 );
  clipper::ResolutionFn Sigma( hklinf, basis_fo, target_fo, params_init );


  // can't seem to get max resolution from clipper, so use CCP4 methods
  CMtz::MTZ *mtz1=NULL;
  int read_refs=1;  // not sure what read_refs actually does
  float minres,maxres;
  mtz1 = CMtz::MtzGet(argv[mtzinarg], read_refs);
  CMtz::MtzResLimits(mtz1,&minres,&maxres);
  printf("\n\nMinimum resolution = %7.3f A\nMaximum resolution = %7.3f A\n\n",1.0/sqrt(minres),1.0/sqrt(maxres));
  if (debug) printf("Minimum resolution = %f \nMaximum resolution = %f \n\n",minres,maxres);
  CSym::CCP4SPG *spg1 = CSym::ccp4spg_load_by_standard_num(CMtz::MtzSpacegroupNumber(mtz1));

  std::cout << "\nSpacegroup: " << spgr.symbol_hm() << " (number " << spgr.descr().spacegroup_number() << ")" << std::endl;

  if (debug) {
      FILE *ftestfile;
      ftestfile = fopen("sigma.txt","w");
      for (int i=0; i<60; i++) {
          double res = maxres * pow( double(i+1)/60.0, 0.666666 );
	      fprintf(ftestfile,"%10.6f %10.6f \n", res, basis_fo.f_s( res, Sigma.params() ));
      }
      fclose(ftestfile);
  }

  char pointgroup[20];
  strcpy(pointgroup,spg1->point_group);
  //pointgroup = spg1->point_group;
  printf("\npointgroup = %s\n\n",pointgroup);

  FILE *floggraph, *flogfile;
  floggraph = fopen("ctruncate.loggraph","w");
  flogfile = fopen("ctruncate.log","w");

  // calculate moments of Z using truncate methods

  int nbins = 60;

  std::vector<double> Na(nbins,0.0),  Nc(nbins,0.0);
  std::vector<double> E1a(nbins,0.0), E1c(nbins,0.0);
  std::vector<double> E3a(nbins,0.0), E3c(nbins,0.0);
  std::vector<double> I1a(nbins,0.0), I1c(nbins,0.0);
  std::vector<double> I2a(nbins,0.0), I2c(nbins,0.0);
  std::vector<double> I3a(nbins,0.0), I3c(nbins,0.0);
  std::vector<double> I4a(nbins,0.0), I4c(nbins,0.0);

  for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
	  if ( !isig[ih].missing() ) {
		  //if( ih.hkl_class().epsilon() > 1.01) printf("%d %d %d epsilon = %f\n",ih.hkl().h(), ih.hkl().k(), ih.hkl().l(), ih.hkl_class().epsilon());
		  double I = isig[ih].I() / ih.hkl_class().epsilon();
		  int bin = int( nbins * pow( ih.invresolsq() / double(maxres), 1.0 ) - 0.5  );
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
		  else {
			  Nc[bin]++;
		      if (I > 0.0) {
			      E1c[bin] += sqrt(I);
			      I1c[bin] += I;
			      E3c[bin] += pow(I,1.5);
			      I2c[bin] += I*I;
			      I3c[bin] += I*I*I;
			      I4c[bin] += I*I*I*I;
			  }
		  }
	  }
  }

  printf("$TABLE: Acentric moments of E using Truncate method:\n");
  printf("$GRAPHS");
  printf(": 1st & 3rd moments of E (Expected values = 0.886, 1.329, Perfect twin = 0.94, 1.175):0|%5.3fx0|2:1,2,3:\n", maxres);
  printf(": 4th moment of E (Expected value = 2, Perfect Twin = 1.5):0|%5.3fx0|5:1,4:\n", maxres);
  //printf(": 6th & 8th moments of E (Expected value = 6, 24, Perfect Twin 3, 7.5):0|%5.3fx0|48:1,5,6:\n", maxres);
  printf("$$ 1/resol^2 <E> <E**3> <E**4> <E**6> <E**8> $$\n$$\n");

  for (int i=0; i<60; i++) {
	  double res = maxres * pow((double(i) + 0.5)/double(nbins), 1.00000);
	  if (Na[i] > 0 && I1a[i] > 0.0) {
		  E1a[i] /= sqrt(I1a[i]*Na[i]);
		  E3a[i] *= sqrt(Na[i]) / pow(I1a[i],1.5);
		  I2a[i] *= Na[i] /( I1a[i]*I1a[i] );
		  I3a[i] *= Na[i]*Na[i] / pow(I1a[i],3);
		  I4a[i] *= pow(Na[i],3) / pow(I1a[i],4);
	  }
	  printf("%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", res, E1a[i], E3a[i], I2a[i], I3a[i], I4a[i]);
  }
  printf("$$\n\n");

  printf("$TABLE: Centric moments of E using Truncate method:\n");
  printf("$GRAPHS");
  printf(": 1st & 3rd moments of E (Expected = 0.798, 1.596, Perfect Twin = 0.886, 1.329):0|%5.3fx0|4:1,2,3:\n", maxres);
  printf(": 4th moment of E (Expected = 3, Perfect Twin = 2):0|%5.3fx0|5:1,4:\n", maxres);
  //printf(": 6th & 8th moments of E (Expected = 15, 105, Perfect Twin = 6, 24):0|%5.3fx0|120:1,5,6:\n", maxres);
  printf("$$ 1/resol^2 <E> <E**3> <E**4> <E**6> <E**8> $$\n$$\n");

  for (int i=0; i<60; i++) {
	  double res = maxres * pow((double(i) + 0.5)/double(nbins), 0.666666);
	  if (Nc[i] > 0 && I1c[i] > 0.0) {
		  E1c[i] /= sqrt(I1c[i]*Nc[i]);
		  E3c[i] *= sqrt(Nc[i]) / pow(I1c[i],1.5);
		  I2c[i] *= Nc[i] /( I1c[i]*I1c[i] );
		  I3c[i] *= Nc[i]*Nc[i] / pow(I1c[i],3);
		  I4c[i] *= pow(Nc[i],3) / pow(I1c[i],4);
	  }
	  printf("%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", res, E1c[i], E3c[i], I2c[i], I3c[i], I4c[i]);
  }
  printf("$$\n\n");

  printf("\nTWINNING ANALYSIS:\n\n");
  int itwin = 0;

  // H test for twinning

  Cell cell = hklinf.cell();
  Mat33<int> twinop(0,0,0,0,0,0,0,0,0);
  int sg = CMtz::MtzSpacegroupNumber(mtz1);
  if ( (sg >= 75 && sg <= 80) || sg == 146 || (sg >= 168 && sg <= 173) || (sg >= 195 && sg <= 199) ) { 
	  printf("twinning operator k, h, -l\n");
	  //printf("twinning operator h, -h-k, -l\n");
	  twinop(0,1) = 1;
	  twinop(1,0) = 1;
	  twinop(2,2) = -1;
	  //twinop(0,0) = 1;
	  //twinop(1,0) = -1;
	  //twinop(1,1) = -1;
	  Htest(isig,twinop,itwin,debug);
  }
  else if( sg >= 149 && sg <= 154 ) {
	  printf("twinning operator -h, -k, l\n");
	  twinop(0,0) = -1;
	  twinop(1,1) = -1;
	  twinop(2,2) = 1;
	  Htest(isig,twinop,itwin,debug);
  }
  else if( sg >= 143 && sg <= 145 ) {
	  printf("twinning operator k, h, -l\n");
	  twinop(0,1) = 1;
	  twinop(1,0) = 1;
	  twinop(2,2) = -1;
	  Htest(isig,twinop,itwin,debug);

	  printf("twinning operator -k, -h, -l\n");
	  twinop(0,1) = -1;
	  twinop(1,0) = -1;
	  Htest(isig,twinop,itwin,debug);

	  printf("twinning operator -h, -k, l\n");
      twinop(0,1) = 0;
	  twinop(1,0) = 0;
	  twinop(0,0) = -1;
	  twinop(1,1) = -1;
	  twinop(2,2) = 1;
	  Htest(isig,twinop,itwin,debug);
  }
  else if( !strcmp(pointgroup, "PG222") ) {
	  //printf("PG222\n");
	  // Can have pseudo-merohedral twinning in PG222 (orthorhombic) if a=b, b=c or c=a
	  if ( fabs( 1.0 - cell.b()/cell.a() ) < 0.02 ) { 
	      printf("twinning operator k, h, -l\n");
	      twinop(0,1) = 1;
	      twinop(1,0) = 1;
	      twinop(2,2) = -1;
	      Htest(isig,twinop,itwin,debug);
	  }
	  if ( fabs( 1.0 - cell.c()/cell.b() ) < 0.02 ) { 
	      printf("twinning operator -h, l, k\n");
	      twinop(0,1) = 0;
	      twinop(1,0) = 0;
	      twinop(2,2) = 0;
		  twinop(0,0) = -1;
		  twinop(1,2) = 1;
		  twinop(2,1) = 1;
	      Htest(isig,twinop,itwin,debug);
	  }
	  if ( fabs( 1.0 - cell.a()/cell.c() ) < 0.02 ) {
	      printf("twinning operator l, -k, h\n");
	      twinop(0,0) = 0;
	      twinop(1,2) = 0;
	      twinop(2,1) = 0;
		  twinop(1,1) = -1;
		  twinop(0,2) = 1;
		  twinop(2,0) = 1;
	      Htest(isig,twinop,itwin,debug);
	  }
  }
  else if( !strcmp(pointgroup, "PG2") ) {
      // can have pseudomerohedral twinning in PG2 (monoclinic) if
	  // beta = 90
	  // a=c
	  // cos(beta) = -a/2c, -c/a, -c/2a
	  // sin(beta) = a/c
	  if ( fabs( 1.0 - cell.a()/cell.c() ) < 0.02  || fabs( sin(cell.beta()) - cell.a()/cell.c() ) < 0.02 ) {
	      printf("twinning operator l, -k, h\n");
		  twinop(1,1) = -1;
		  twinop(0,2) = 1;
		  twinop(2,0) = 1;
	      Htest(isig,twinop,itwin,debug);
	  }

	  if ( cell.beta() < Util::d2rad(93.0) ) {
		  printf("twinning operator h, -k, l\n");
		  twinop(0,0) = 1;
		  twinop(0,2) = 0;
		  twinop(1,1) = -1;
		  twinop(2,0) = 0;
		  twinop(2,2) = 1;
		  Htest(isig,twinop,itwin,debug);
	  }	 

	  if ( fabs( cos(cell.beta()) + 0.5*cell.a()/cell.c() ) < 0.02 ) { 
		  printf("twinning operator -h, -k, h+l\n");
		  twinop(0,0) = -1;
		  twinop(0,2) = 0;
		  twinop(1,1) = -1;
		  twinop(2,0) = 1;
		  twinop(2,2) = 1;
		  Htest(isig,twinop,itwin,debug);
	  }	  

	  if ( fabs( cos(cell.beta()) + 0.5*cell.c()/cell.a() ) < 0.02 ) { 
		  printf("twinning operator h+l, -k, -l\n");
		  twinop(0,0) = 1;
		  twinop(0,2) = 1;
		  twinop(1,1) = -1;
		  twinop(2,0) = 0;
		  twinop(2,2) = -1;
		  Htest(isig,twinop,itwin,debug);
	  }	  
	  if ( fabs( cos(cell.beta()) + cell.c()/cell.a() ) < 0.02 ) { 
		  printf("twinning operator h+2l, -k, -l\n");
		  twinop(0,0) = 1;
		  twinop(0,2) = 2;
		  twinop(1,1) = -1;
		  twinop(2,0) = 0;
		  twinop(2,2) = -1;
		  Htest(isig,twinop,itwin,debug);
	  }
  }

  // L test for twinning

  double LT=0.0;
  double LT2=0.0;
  double NLT=0.0;
  std::vector<int> cdf(20,0);

  for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
	  if ( !isig[ih].missing() && !ih.hkl_class().centric() ) {
		  HKL hkl = ih.hkl();
		  int h = hkl.h();
		  int k = hkl.k();
		  int l = hkl.l();
		  for ( int delta1 = -2; delta1 <= 2; delta1 += 2 ) {
			  for ( int delta2 = -2; delta2 <= 2; delta2 += 2 ) {
				  for ( int delta3 = -2; delta3 <= 2; delta3 += 2 ) {
					  HKL hkl2;
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
		  HKL hkl = ih.hkl();
			              if (fabs(L) < 1){
			                  LT += fabs(L)*weight;
			                  LT2 += L*L*weight;
			                  NLT += weight;
							  for (int i=0;i<20;i++) {
								  if ( fabs(L) < (double(i+1))/20.0 ) cdf[i]++;
					  HKL hkl2;
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
	  printf("\nApplying the L test for twinning: (Padilla and Yeates Acta D59 1124 (2003))\n");
	  printf("L test suggests data is twinned\n");
	  printf("L statistic = %6.3f  (untwinned 0.5 perfect twin 0.375)\n\n", Lav);
	  itwin = 1;
	  printf("All data regardless of I/sigma(I) has been included in the L test\n");
	  printf("Anisotropy correction has not been applied before calculating L\n\n");
  }

  if (!itwin) printf("No twinning detected\n\n");

  printf("$TABLE: L test for twinning:\n");
  printf("$GRAPHS");
  printf(": cumulative distribution function for |L|:0|1x0|1:1,2,3,4:\n");
  printf("$$ |L| untwinned perfect_twin data $$\n$$\n");
  printf("0.000000 0.000000 0.000000 0.000000\n");

  for (int i=0;i<20;i++) {
	  double x = (double(i+1))/20.0;
	  printf("%f %f %f %f\n", x, x, 0.5*x*(3.0-x*x),  double(cdf[i])/NLT);
  }
  printf("$$\n\n");
   

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
		  HKL hkl = ih.hkl();
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
		  HKL hkl = ih.hkl();
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
  fclose(flogfile);


  //Wilson pre
  
  printf("\nWILSON SCALING:\n\n");
  int nsym = spg1->nsymop;
  int nresidues = int(0.5*hklinf.cell().volume()/(nsym*157));
  //nresidues /= 2.0;
  printf("Estimated number of residues = %d\n",nresidues);

  std::string name[4] = { "C", "N", "O", "H" };
  int numatoms[4]; 
  numatoms[0] = 5*nresidues;
  numatoms[1] = int(1.35*nresidues);
  numatoms[2] = int(1.5*nresidues);
  numatoms[3] = 8*nresidues;


  // Wilson plot
  
  std::vector<float> xi, yi, wi; 
  std::vector<int> flags(nbins,0);
  float totalscat; 
  float minres_scaling = 0.0625;   // 4 Angstroms

  printf("$TABLE: Wilson plot:\n");
  printf("$SCATTER");
  //printf(": Wilson plot:0|0.1111x-7|-5:1,2:\n$$");  // limits hardwired
  printf(": Wilson plot:A:1,2:\n$$");  
  printf(" 1/resol^2 ln(I/I_th) $$\n$$\n");

  for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
	  if ( !isig[ih].missing() ) {
		  float lnS = -log(Sigma.f(ih));
		  float res = ih.invresolsq();

		  totalscat = 0;
		  for (int i=0;i<4;i++) {
			  Atom atom;
              atom.set_occupancy(1.0);
              atom.set_element(name[i]);
              atom.set_u_iso(0.0);
              atom.set_u_aniso_orth( U_aniso_orth( U_aniso_orth::null() ) ); // need this o/w next line hangs
              AtomShapeFn sf(atom);
			  float scat = sf.f(res);
			  totalscat +=  float( nsym * numatoms[i] ) * scat * scat;
		  }
		  lnS += log(totalscat);

		  int bin = int( nbins * pow( ih.invresolsq() / double(maxres), 1.5 ) - 0.5  );
		  if (bin >= nbins || bin < 0) printf("Warning: (Wilson) illegal bin number %d\n", bin);
	      if (flags[bin] != 1) {
		      printf("%10.5f %10.5f\n", res, -lnS);
		      flags[bin] = 1;
	      }
		  if (res > minres_scaling) {  
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
  printf("$$\n\n");

  int nobs = xi.size();
  //printf("%d %d %d\n", xi.size(), yi.size(), wi.size());
  float a,b,siga,sigb,a1,b1;
  b = 0.0;

  if ( wi.size() > 200 && maxres > 0.0816) {               // 3.5 Angstroms
      straight_line_fit(xi,yi,wi,nobs,a,b,siga,sigb);
	  printf("\nResults from Clipper style Wilson plot:\n");
      a *= 2.0;
      printf ("B = %6.3f intercept = %6.3f siga = %6.3f sigb = %6.3f\n",a,b,siga,sigb);
      printf("scale factor on intensity = %10.4f\n\n",(exp(b)));
  }
  else {
	  printf("Too few high resolution points to determine B factor and Wilson scale factor\n");
  }

  // Truncate style Wilson plot

  std::vector<int> N_all(nbins,0);
  std::vector<int> N_obs(nbins,0);
  std::vector<float> I_obs(nbins,0.0);

  std::vector<float> xtr, ytr, wtr; 

  printf("$TABLE: Truncate style Wilson plot:\n");
  printf("$GRAPHS");
  printf(": Wilson plot:A:1,2,3:\n$$");  
  //printf(": Wilson plot:0|0.1111x-8|-5.5:1,2,3:\n$$");  // limits hardwired
  printf(" 1/resol^2 obs all $$\n$$\n");

  for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
      int bin = int( float(nbins) * ih.invresolsq() / maxres - 0.5 );
	  if (bin >= nbins || bin < 0) printf("Warning: (Wilson 2) illegal bin number %d\n", bin);
	  N_all[bin]++;
	  if ( !isig[ih].missing() ) {
		  I_obs[bin] += isig[ih].I();
		  N_obs[bin]++;
	  }
  }

  for ( int j=0; j<nbins; j++ ) {
	  float res = maxres*(float(j)+0.5)/float(nbins);
	  totalscat = 0;
	  for (int i=0;i<4;i++) {
		  Atom atom;
          atom.set_occupancy(1.0);
          atom.set_element(name[i]);
          atom.set_u_iso(0.0);
		  atom.set_u_aniso_orth( U_aniso_orth( U_aniso_orth::null() ) ); // need this o/w next line hangs
          AtomShapeFn sf(atom);
	      float scat = sf.f(res);
		  totalscat +=  float( nsym * numatoms[i] ) * scat * scat;
	  }
	  float w1 = log( I_obs[j] / (float(N_obs[j]) * totalscat) );
	  float w2 = log( I_obs[j] / (float(N_all[j]) * totalscat) );
	  printf( "%10.6f %10.6f %10.6f\n", res, w1, w2);
	  if (res > minres_scaling) {  
		  xtr.push_back(res);
		  ytr.push_back(w1);
		  wtr.push_back(1.0);
	  }
  }
  printf("$$\n\n");

  nobs = xtr.size();
  if ( wi.size() > 200 && maxres > 0.0816) {
      straight_line_fit(xtr,ytr,wtr,nobs,a1,b1,siga,sigb);
      printf("\nresults from fitting Truncate style Wilson plot\n");
      printf ("B = %6.3f intercept = %6.3f siga = %6.3f sigb = %6.3f\n",-2.0*a1,-b1,siga,sigb);
      printf("scale factor on intensity = %10.4f\n\n", exp(-b1));
  }
 
  // apply the Truncate procedure

  float scalef = sqrt(exp(b));
  //float scalef = 4.66047;  //hardwired (for now) scalefactor

  if (anomalous) {
	  truncate( isig_plus, jsig_plus, fsig_plus, Sigma, scalef, spg1 );
      truncate( isig_minus, jsig_minus, fsig_minus, Sigma, scalef, spg1 );
	  int iwarn = 0;
	  for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
	       if ( !isig_plus[ih].missing() && !isig_minus[ih].missing() ) {
			   fsig[ih].f() = 0.5 * ( fsig_plus[ih].f() + fsig_minus[ih].f() );
			   fsig[ih].sigf() = 0.5 * sqrt( pow( fsig_plus[ih].sigf(), 2 ) + pow( fsig_minus[ih].sigf(), 2 ) );
			   Dano[ih].f() = fsig_plus[ih].f() - fsig_minus[ih].f();
			   Dano[ih].sigf() = 2.0 * fsig[ih].sigf();
		   }
		   else if ( !isig_plus[ih].missing() ) {
			   fsig[ih].f() = fsig_plus[ih].f();
			   fsig[ih].f() = fsig_plus[ih].sigf();
		   }
		   else if ( !isig_minus[ih].missing() ) {
			   fsig[ih].f() = fsig_minus[ih].f();
			   fsig[ih].f() = fsig_minus[ih].sigf();
		   }
		   else if ( !isig[ih].missing() && iwarn != 1 ) {
			   printf("\nWARNING: Imean exists but I(+), I(-) do not\n\n");
			   iwarn = 1;
		   }
		   if ( ih.hkl_class().centric() ) {
			   Dano[ih].f() = 0.0;
			   Dano[ih].sigf() = 0.0;
		   }
	  }
  }
  else {
      truncate( isig, jsig, fsig, Sigma, scalef, spg1 );
  }


  
  // following code is for when truncate calc switched off - do not delete it
  // usually already have F's in this case; really just need to skip truncate calc

  /*
  for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
	  if ( !isig[ih].missing() ) {
		  float I = isig[ih].I();
		  float sigma = isig[ih].sigI();
		  HKL hkl = ih.hkl();
		  float weight = (float) CSym::ccp4spg_get_multiplicity( spg1, hkl.h(), hkl.k(), hkl.l() );
		  float sqwt = sqrt(weight);
		  if (I < 0.0) {
			  fsig[ih].f() = 0.0;
			  fsig[ih].sigf() = 0.0;
		  }
		  else {
			  fsig[ih].f() = sqrt(I)*scalef*sqwt;
			  fsig[ih].sigf() = 0.5*(sigma/sqrt(I))*scalef*sqwt; //check this
		  }
	  }
  }*/


  // create initial E
  clipper::HKL_data<clipper::data32::E_sigE> esig( hklinf );
  esig.compute( fsig, clipper::data32::Compute_EsigE_from_FsigF() );

  // calc E-scaling
  std::vector<double> params_init1( nprm, 1.0 );
  clipper::BasisFn_spline basis_fo1( esig, nprm, 2.0 );
  clipper::TargetFn_scaleEsq<clipper::data32::E_sigE> target_fo1( esig );
  clipper::ResolutionFn escale( hklinf, basis_fo1, target_fo1, params_init1 );

  // apply E-scaling
  for ( HRI ih = esig.first(); !ih.last(); ih.next() )
    if ( !esig[ih].missing() ) esig[ih].scale( sqrt( escale.f(ih) ) );

  //const int nprm2 = 60; 

  clipper::HKL_data<clipper::data32::E_sigE> esig1 = esig;
  clipper::HKL_data<clipper::data32::E_sigE> esig2 = esig;

  for ( HRI ih = esig.first(); !ih.last(); ih.next() ) {
	  if( ih.hkl_class().centric() ) esig1[ih].set_null();  // esig1 is acentrics
	  if( !ih.hkl_class().centric() ) esig2[ih].set_null(); // esig2 is centrics
  }


  // moments of E using clipper binning

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
    
  int nprmk = 60;
  std::vector<double> params_initk( nprmk, 1.0 );
  clipper::BasisFn_binner basis_fn1( fs1, nprmk, 1.5 );  // equal increments in invresolsq bins
  //clipper::BasisFn_binner basis_fn1( fs1, nprmk, 1.0 );  // equal volume bins


  clipper::TargetFn_meanFnth<clipper::data32::F_sigF> target_fn1( fs1, 1.0 );
  clipper::TargetFn_meanFnth<clipper::data32::F_sigF> target_fn2( fs1, 2.0 );
  clipper::TargetFn_meanFnth<clipper::data32::F_sigF> target_fn3( fs1, 3.0 );
  clipper::TargetFn_meanFnth<clipper::data32::F_sigF> target_fn4( fs1, 4.0 );
  clipper::TargetFn_meanFnth<clipper::data32::F_sigF> target_fn6( fs1, 6.0 );
  clipper::TargetFn_meanFnth<clipper::data32::F_sigF> target_fn8( fs1, 8.0 );

  clipper::ResolutionFn f1( hklinf, basis_fn1, target_fn1, params_initk );
  clipper::ResolutionFn f2( hklinf, basis_fn1, target_fn2, params_initk );
  clipper::ResolutionFn f3( hklinf, basis_fn1, target_fn3, params_initk );
  clipper::ResolutionFn f4( hklinf, basis_fn1, target_fn4, params_initk );
  clipper::ResolutionFn f6( hklinf, basis_fn1, target_fn6, params_initk );
  clipper::ResolutionFn f8( hklinf, basis_fn1, target_fn8, params_initk );
   
  //printf("$TABLE: Acentric moments of E for k=1,3,4,6,8:\n");
  printf("$TABLE: Acentric moments of E for k=1,3,4:\n");
  printf("$GRAPHS");
  printf(": 1st & 3rd moments of E (Expected values = 0.886, 1.329, Perfect twin = 0.94, 1.175):0|%5.3fx0|2:1,2,3:\n", maxres);
  printf(": 4th moment of E (Expected value = 2, Perfect Twin = 1.5):0|%5.3fx0|5:1,4:\n", maxres);
  //printf(": 6th & 8th moments of E (Expected value = 6, 24, Perfect Twin 3, 7.5):0|%5.3fx0|48:1,5,6:\n", maxres);

  //printf("$$ 1/resol^2 <E> <E**3> <E**4> <E**6> <E**8> $$\n$$\n");
  printf("$$ 1/resol^2 <E> <E**3> <E**4> $$\n$$\n");

  double mean1, mean3, mean4, mean6, mean8;
  mean1 = mean3 = mean4 = mean6 = mean8 = 0.0;
  for (int i=0; i<nbins; i++) {
	  double res = double(i+1) * maxres / double(nbins);   // equal increments in invresolsq bins
	  //double res = maxres * pow( double(i+1)/double(nbins), 0.666666 );  // equal volume bins
	  double i1 = basis_fn1.f_s( res, f2.params() );
	  printf("%10.6f %10.6f %10.6f %10.6f \n", res,
              basis_fn1.f_s( res, f1.params() )/pow(i1,0.5),
              basis_fn1.f_s( res, f3.params() )/pow(i1,1.5), 
	          basis_fn1.f_s( res, f4.params() )/pow(i1,2.0) ); 
              //basis_fn1.f_s( res, f6.params() )/pow(i1,3.0), 
              //basis_fn1.f_s( res, f8.params() )/pow(i1,4.0) ); 

	  mean1 += basis_fn1.f_s( res, f1.params() )/pow(i1,0.5);
      mean3 += basis_fn1.f_s( res, f3.params() )/pow(i1,1.5); 
	  mean4 += basis_fn1.f_s( res, f4.params() )/pow(i1,2.0); 
      mean6 += basis_fn1.f_s( res, f6.params() )/pow(i1,3.0); 
      mean8 += basis_fn1.f_s( res, f8.params() )/pow(i1,4.0);
  }
  mean1 /= double(nbins);
  mean3 /= double(nbins);
  mean4 /= double(nbins);
  mean6 /= double(nbins);
  mean8 /= double(nbins);

  printf("\nMEAN ACENTRIC MOMENTS OF E:\n\n");
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

  printf("$$\n\n");

  //centrics
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
  clipper::TargetFn_meanFnth<clipper::data32::F_sigF> target_fn6c( fs2, 6.0 );
  clipper::TargetFn_meanFnth<clipper::data32::F_sigF> target_fn8c( fs2, 8.0 );

  clipper::ResolutionFn f1c( hklinf, basis_fn2, target_fn1c, params_initk2 );
  clipper::ResolutionFn f2c( hklinf, basis_fn2, target_fn2c, params_initk2 );
  clipper::ResolutionFn f3c( hklinf, basis_fn2, target_fn3c, params_initk2 );
  clipper::ResolutionFn f4c( hklinf, basis_fn2, target_fn4c, params_initk2 );
  clipper::ResolutionFn f6c( hklinf, basis_fn2, target_fn6c, params_initk2 );
  clipper::ResolutionFn f8c( hklinf, basis_fn2, target_fn8c, params_initk2 );
   
  //printf("$TABLE: Centric moments of E for k=1,3,4,6,8:\n");
  printf("$TABLE: Centric moments of E for k=1,3,4:\n");
  printf("$GRAPHS");
  printf(": 1st & 3rd moments of E (Expected = 0.798, 1.596, Perfect Twin = 0.886, 1.329):0|%5.3fx0|4:1,2,3:\n", maxres);
  printf(": 4th moment of E (Expected = 3, Perfect Twin = 2):0|%5.3fx0|5:1,4:\n", maxres);
  //printf(": 6th & 8th moments of E (Expected = 15, 105, Perfect Twin = 6, 24):0|%5.3fx0|120:1,5,6:\n", maxres);

  //printf("$$ 1/resol^2 <E> <E**3> <E**4> <E**6> <E**8> $$\n$$\n");
  printf("$$ 1/resol^2 <E> <E**3> <E**4> $$\n$$\n");


  for (int i=0; i<nbins; i++) {
	  double res = double(i+1) * maxres / double(nbins);   // equal increments in invresolsq bins
	  //double res = maxres * pow( double(i+1)/double(nbins), 0.666666 );  // equal volume bins
	  double i1 = basis_fn2.f_s( res, f2c.params() );
	  printf("%10.6f %10.6f %10.6f %10.6f\n", res,
              basis_fn2.f_s( res, f1c.params() )/pow(i1,0.5),
              basis_fn2.f_s( res, f3c.params() )/pow(i1,1.5),
	  	      basis_fn1.f_s( res, f4c.params() )/pow(i1,2.0) ); 
              //basis_fn1.f_s( res, f6c.params() )/pow(i1,3.0), 
              //basis_fn1.f_s( res, f8c.params() )/pow(i1,4.0) ); 
  }

  printf("$$\n\n");


  // find the range of intensities
  clipper::Range<double> intensity_range_centric;
  clipper::Range<double> intensity_range_acentric;

  for ( HRI ih = esig.first(); !ih.last(); ih.next() ) {
	  if ( !esig[ih].missing() ) {
         if ( ih.hkl_class().centric() ) intensity_range_centric.include( pow(esig[ih].E(),2) );
		 else intensity_range_acentric.include( pow(esig[ih].E(),2) );
	  }
  }
  //printf("C: %20.5f %20.5f\n",intensity_range_centric.max(),intensity_range_centric.min());
  //printf("A: %20.5f %20.5f\n",intensity_range_acentric.max(),intensity_range_acentric.min());


  // construct cumulative distribution function for intensity
  clipper::Generic_ordinal intensity_ord_c;
  clipper::Generic_ordinal intensity_ord_a;
  intensity_ord_c.init( intensity_range_centric );
  intensity_ord_a.init( intensity_range_acentric );
  for ( HRI ih = esig.first(); !ih.last(); ih.next() ) {
  	  if ( !esig[ih].missing() ) {
          if ( ih.hkl_class().centric() ) intensity_ord_c.accumulate( esig[ih].E()*esig[ih].E() );
		  else intensity_ord_a.accumulate( esig[ih].E()*esig[ih].E() );
	  }
  }
  intensity_ord_c.prep_ordinal();
  intensity_ord_a.prep_ordinal();


  // construct another cumulative distribution function for intensity (using Z rather than E)
  clipper::Range<double> intensity_range_centric2;
  clipper::Range<double> intensity_range_acentric2;
  // changed from jsig to isig

  for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
	  if ( !isig[ih].missing() ) {
         if ( ih.hkl_class().centric() ) intensity_range_centric2.include( isig[ih].I()/Sigma.f(ih) );
		 else intensity_range_acentric2.include( isig[ih].I()/Sigma.f(ih) );
	  }
  }
  //printf("C2: %20.5f %20.5f\n",intensity_range_centric2.max(),intensity_range_centric2.min());
  //printf("A2: %20.5f %20.5f\n",intensity_range_acentric2.max(),intensity_range_acentric2.min());

  clipper::Generic_ordinal intensity_ord_c2;
  clipper::Generic_ordinal intensity_ord_a2;
  intensity_ord_c2.init( intensity_range_centric2 );
  intensity_ord_a2.init( intensity_range_acentric2 );
  for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
  	  if ( !isig[ih].missing() ) {
          if ( ih.hkl_class().centric() ) intensity_ord_c2.accumulate( isig[ih].I()/Sigma.f(ih)  );
		  else intensity_ord_a2.accumulate( isig[ih].I()/Sigma.f(ih) );
	  }
  }
  intensity_ord_c2.prep_ordinal();
  intensity_ord_a2.prep_ordinal();

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



  printf("$TABLE: Cumulative intensity distribution:\n");
  printf("$GRAPHS");
  printf(": Cumulative intensity distribution (Acentric and centric):N:1,2,3,4,5:\n$$");
  printf(" Z Acent_theor Acent_obser Cent_theor Cent_obser $$\n$$\n");
  double x = 0.0;
  double deltax=0.04;
  for (int i=0; i<=50; i++) {
	  printf("%10.5f %8.5f %8.5f %8.5f %8.5f\n", x, acen[i], intensity_ord_a.ordinal(x), cen[i], intensity_ord_c.ordinal(x));
	  x += deltax;
  }
  printf("$$\n\n");

// count entries where acentric distribution lower than expected - sign of twinning
  x = 0.08;
  deltax = 0.12;
  int ntw = 0;
  for (int i=0;i<4; i++) {
	  if ( (acen[3*i+2] - intensity_ord_a.ordinal(x))/acen[3*i+2] > 0.4 ) ntw ++;
  }
  if (ntw > 2) printf("\nWARNING: ****  Cumulative Distribution shows Possible Twinning ****\n");


  // falloff calculation (Yorgo Modis)

  Mat33<float> transf;
  //Cell cell = hklinf.cell();

  // calculate the matrix that transforms the cell coordinates h,k,l to Cartesian coordinates.
  tricart (cell, transf);

  nbins = 60;
  std::vector<float> somov(nbins,0.0);
  std::vector<float> somsdov(nbins,0.0);
  std::vector<int> numov(nbins,0);
  std::vector<float> enumov(nbins,0.0);

  /*std::vector<Vec3<> > somdir(nbins, Vec3<> (0.0,0.0,0.0));
  std::vector<Vec3<> > somsddir(nbins, Vec3<> (0.0,0.0,0.0));
  std::vector<Vec3<> > numdir(nbins, Vec3<> (0,0,0)); // will this work for integer?
  std::vector<Vec3<> > enumdir(nbins, Vec3<> (0.0,0.0,0.0)); */
  //std::vector<float> somdir[3](nbins,0.0);//,(nbins,0.0),(nbins,0.0)};
  //std::vector<float> somsddir(nbins,0.0);
  //std::vector<int> numdir(nbins,0);
  //std::vector<float> enumdir(nbins,0.0);

  float somdir[3][60];
  float somsddir[3][60];
  int numdir[3][60];
		  HKL hkl = ih.hkl();

  for (int i=0;i<3;i++){
	  for (int j=0;j<60;j++){
		  somdir[i][j] = somsddir[i][j] = enumdir[i][j] = 0.0;
		  numdir[i][j] = 0;
	  }
  }

  float cone = 30.0; //hardwired for now
  float ang[3];
  float cosang;
  Vec3<int> jhkl;
  Vec3<float> jhkl2;

  int nsymp = spg1->nsymop_prim;
  //printf("nsymp = %d  nsym = %d\n",nsymp,nsym);

  for ( HRI ih = fsig.first(); !ih.last(); ih.next() ) {
	  if ( !fsig[ih].missing() ) {
		  // bin number different in C because arrays start at zero
		  int bin = int( double(nbins) * ih.invresolsq() / maxres - 0.001);
		  if (bin >= nbins || bin < 0) printf("Warning: (Modis) illegal bin number %d\n", bin);
		  HKL hkl = ih.hkl();
		  float epsiln = (float) CSym::ccp4spg_get_multiplicity( spg1, hkl.h(), hkl.k(), hkl.l() ); 
		  epsiln /= spg1->nsymop;
		  if ( !ih.hkl_class().centric() ) epsiln *= 0.5;
		  //printf("%8.4f %8.4f\n", epsiln, ih.hkl_class().epsilonc());
          
		  for (int isym = 1; isym <= nsymp*2; isym++) {
			  CSym::ccp4spg_generate_indices(spg1, isym, hkl.h(), hkl.k(), hkl.l(), &jhkl[0], &jhkl[1], &jhkl[2]);

// convert h,k,l into Cartesian coordinates.
			  for (int i=0; i<3; i++) { jhkl2[i] = (float) jhkl[i];}
   
			  Vec3<float> hc = transf * jhkl2;
 
			  for (int j=0;j<3;j++) {
                  cosang = fabs(hc[j])/sqrt(ih.invresolsq());
// cosang can stray just past 1.0
                  if (cosang > 1.0) cosang = 1.0;
                  ang[j] = acos(cosang);
				  if ( ang[j] < Util::d2rad(cone) ) {
                      somdir[j][bin] += fsig[ih].f()*epsiln;
					  somsddir[j][bin] += epsiln*fsig[ih].f()/fsig[ih].sigf();
                      enumdir[j][bin] += epsiln;
                      numdir[j][bin]++;
			      }
			  }
              somov[bin] += fsig[ih].f()*epsiln;
              somsdov[bin] += epsiln*fsig[ih].f()/fsig[ih].sigf();
              enumov[bin] += epsiln;
              numov[bin]++;
		  }
	  }
  }

  for (int i=0;i<60;i++) {
	  for (int j=0;j<3;j++) {
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

  printf("$TABLE: Anisotropy analysis (Yorgo Modis):\n");
  printf("$GRAPHS");
  printf(": Mn(F) v resolution:N:1,2,3,4,5:\n");
  printf(": Mn(F/sd) v resolution:N:1,6,7,8,9:\n");
  printf(": No. reflections v resolution:N:1,10,11,12,13:\n");
  printf("$$ 1/resol^2 Mn(F(d1)) Mn(F(d2)) Mn(F(d3)) Mn(F(ov) Mn(F/sd(d1)) Mn(F/sd(d2)) Mn(F/sd(d3)) Mn(F/sd(ov))");
  printf(" N(d1) N(d2) N(d3) N(ov) $$\n$$\n");


  for(int i=0;i<nbins;i++){
	  double res = maxres*(double(i)+0.5)/double(nbins);
	  printf("%10.6f %12.4e %12.4e %12.4e %12.4e ",res,somdir[0][i],somdir[1][i],somdir[2][i],somov[i]);
	  printf("%12.4e %12.4e %12.4e %12.4e ",somsddir[0][i],somsddir[1][i],somsddir[2][i],somsdov[i]);
	  printf("%8d %8d %8d %8d\n",numdir[0][i],numdir[1][i],numdir[2][i],numov[i]);
  }
  printf("$$\n\n");

  fclose(floggraph);

  // output data

  mtzout.open_append( argv[mtzinarg], outfile );
  //mtzout.export_hkl_data( jsig, outcol );
  mtzout.export_hkl_data( fsig, outcol+"MEAN" );
  if (anomalous) {
	   mtzout.export_hkl_data( fsig_plus, outcol+"(+)" );
	   mtzout.export_hkl_data( fsig_minus, outcol+"(-)" );
	   mtzout.export_hkl_data( Dano, outcol+"_ANO" );
  }
  mtzout.close_append();

  return(0);
}


int truncate(  HKL_data<data32::I_sigI> isig,   HKL_data<data32::I_sigI>& jsig,   HKL_data<data32::F_sigF>& fsig,
			   clipper::ResolutionFn Sigma, float scalef, CSym::CCP4SPG *spg1)
{
  typedef clipper::HKL_data_base::HKL_reference_index HRI;
  //FILE *checkfile;
  //checkfile = fopen("checku.txt", "w");
  float J, sigJ, F, sigF;
  int iflag;

  for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
	  if ( !isig[ih].missing() ) {
		  float I = isig[ih].I();
		  float sigma = isig[ih].sigI();
		  float S = Sigma.f(ih);
		  HKL hkl = ih.hkl();
		  float weight = (float) CSym::ccp4spg_get_multiplicity( spg1, hkl.h(), hkl.k(), hkl.l() );
		  if( fabs( ih.hkl_class().epsilon() - weight ) > 0.001) printf("epsilon %f != weight %f", ih.hkl_class().epsilon(), weight);
		  float sqwt = sqrt(weight);

		  I /= weight;
		  sigma /= weight;

		  // handle acentric and centric reflections separately
		  if ( ih.hkl_class().centric() ) iflag = truncate_centric(I,sigma,S,J,sigJ,F,sigF);
		  else iflag = truncate_acentric(I,sigma,S,J,sigJ,F,sigF);	
		  //if ( !ih.hkl_class().centric()  && I < 0 ) 
			  //fprintf(checkfile,"%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %8.4f %8.4f %8.4f\n", I,sigma,S,J,sigJ,F,sigF,weight,
			  //ih.hkl_class().epsilon());
		  if (iflag) {
			  jsig[ih].I() = J;
			  jsig[ih].sigI() = sigJ;
			  fsig[ih].f() = F*scalef*sqwt;
			  fsig[ih].sigf() = sigF*scalef*sqwt;
			  //fprintf(checkfile,"%12.6f %12.6f %12.6f\n", I,fsig[ih].f(),fsig[ih].sigf());
		  }
		  HKL hkl = ih.hkl();
  }
  //fclose(checkfile);
  return(1);
}


int truncate_acentric(float I, float sigma, float S, float& J, float& sigJ, float& F, float& sigF)
{
  // look up tables taken from truncate.f 
  // tables give values from h = -4.0 to h = 3.0 in steps of 0.1
  // declare as doubles for now, to avoid compiler warnings  

  double ZJ[71] = 
               {0.226,0.230,0.235,0.240,0.246,0.251,0.257,0.263,0.270,
          0.276,0.283,0.290,0.298,0.306,0.314,0.323,0.332,0.341,0.351,
          0.362,0.373,0.385,0.397,0.410,0.424,0.439,0.454,0.470,0.487,
          0.505,0.525,0.545,0.567,0.590,0.615,0.641,0.668,0.698,0.729,
          0.762,0.798,0.835,0.875,0.917,0.962,1.009,1.059,1.112,1.167,
          1.226,1.287,1.352,1.419,1.490,1.563,1.639,1.717,1.798,1.882,
          1.967,2.055,2.145,2.236,2.329,2.422,2.518,2.614,2.710,2.808,
		  2.906,3.004};

  double ZJSD[71] =
               {0.217,0.221,0.226,0.230,0.235,0.240,0.245,0.250,0.255,
          0.261,0.267,0.273,0.279,0.286,0.292,0.299,0.307,0.314,0.322,
          0.330,0.339,0.348,0.357,0.367,0.377,0.387,0.398,0.409,0.421,
          0.433,0.446,0.459,0.473,0.488,0.503,0.518,0.535,0.551,0.568,
          0.586,0.604,0.622,0.641,0.660,0.679,0.698,0.718,0.737,0.757,
          0.776,0.795,0.813,0.831,0.848,0.865,0.881,0.895,0.909,0.921,
          0.933,0.943,0.953,0.961,0.968,0.974,0.980,0.984,0.988,0.991,
		  0.994,0.996};

  double ZF[71] =
               {0.423,0.428,0.432,0.437,0.442,0.447,0.453,0.458,0.464,
          0.469,0.475,0.482,0.488,0.495,0.502,0.509,0.516,0.524,0.532,
          0.540,0.549,0.557,0.567,0.576,0.586,0.597,0.608,0.619,0.631,
          0.643,0.656,0.670,0.684,0.699,0.714,0.730,0.747,0.765,0.783,
          0.802,0.822,0.843,0.865,0.887,0.911,0.935,0.960,0.987,1.014,
          1.042,1.070,1.100,1.130,1.161,1.192,1.224,1.257,1.289,1.322,
          1.355,1.388,1.421,1.454,1.487,1.519,1.551,1.583,1.615,1.646,
		  1.676,1.706};

  double ZFSD[71] = 
               {0.216,0.218,0.220,0.222,0.224,0.226,0.229,0.231,0.234,
          0.236,0.239,0.241,0.244,0.247,0.250,0.253,0.256,0.259,0.262,
          0.266,0.269,0.272,0.276,0.279,0.283,0.287,0.291,0.295,0.298,
          0.302,0.307,0.311,0.315,0.319,0.324,0.328,0.332,0.337,0.341,
          0.345,0.349,0.353,0.357,0.360,0.364,0.367,0.369,0.372,0.374,
          0.375,0.376,0.377,0.377,0.377,0.376,0.374,0.372,0.369,0.366,
          0.362,0.358,0.353,0.348,0.343,0.338,0.332,0.327,0.321,0.315,
		  0.310,0.304};

  float h,x,delta;
  int n;
  
  // Bayesian statistics tells us to modify I/sigma by subtracting off sigma/S
  // where S is the mean intensity in the resolution shell
  h = I/sigma - sigma/S;
  // reject as unphysical reflections for which I < -3.7 sigma, or h < -4.0
  if (I/sigma < -3.7 || h < -4.0 ) {
	  printf("unphys: %f %f %f %f\n",I,sigma,S,h);
	  return(0);
  }
  else {
	  if (h < 3.0) {
		  // use look up table if -4.0 < h < 3.0
          x = 10.0*(h+4.0);
          n = int(x);
          delta = x-n;
		  // linear interpolation
          J = (1.0-delta)*ZJ[n] + delta*ZJ[n+1];
          sigJ = (1.0-delta)*ZJSD[n] + delta*ZJSD[n+1];
          F = (1.0-delta)*ZF[n] + delta*ZF[n+1];
          sigF = (1.0-delta)*ZFSD[n] + delta*ZFSD[n+1];
		  // look up table gives J/sigma, so multiply by sigma to get output intensity
		  J *= sigma;
		  sigJ *= sigma;
		  F *= sqrt(sigma);
		  sigF *= sqrt(sigma);
	  }
	  else {
		  // if h > 4.0 intensities are unchanged by truncate
		  J = h*sigma;
		  sigJ = sigma;
		  F = sqrt(J);
		  sigF = 0.5*sigma/F;
	  }
	  return(1);
  }
}


int truncate_centric(float I, float sigma, float S, float& J, float& sigJ, float& F, float& sigF)
{
  // look up tables taken from truncate.f 
  // tables give values from h = -4.0 to h = 4.0 in steps of 0.1

	  double ZJ[81] = 
		       {0.114,0.116,0.119,0.122,0.124,0.127,0.130,0.134,0.137,
          0.141,0.145,0.148,0.153,0.157,0.162,0.166,0.172,0.177,0.183,
          0.189,0.195,0.202,0.209,0.217,0.225,0.234,0.243,0.253,0.263,
          0.275,0.287,0.300,0.314,0.329,0.345,0.363,0.382,0.402,0.425,
          0.449,0.475,0.503,0.534,0.567,0.603,0.642,0.684,0.730,0.779,
          0.833,0.890,0.952,1.018,1.089,1.164,1.244,1.327,1.416,1.508,
          1.603,1.703,1.805,1.909,2.015,2.123,2.233,2.343,2.453,2.564,
          2.674,2.784,2.894,3.003,3.112,3.220,3.328,3.435,3.541,3.647,
	      3.753,3.863};  // last value corrected (was 3.962 in truncate.f)

      double ZJSD[81] =
	           {0.158,0.161,0.165,0.168,0.172,0.176,0.179,0.184,0.188,
          0.192,0.197,0.202,0.207,0.212,0.218,0.224,0.230,0.236,0.243,
          0.250,0.257,0.265,0.273,0.282,0.291,0.300,0.310,0.321,0.332,
          0.343,0.355,0.368,0.382,0.397,0.412,0.428,0.445,0.463,0.481,
          0.501,0.521,0.543,0.565,0.589,0.613,0.638,0.664,0.691,0.718,
          0.745,0.773,0.801,0.828,0.855,0.881,0.906,0.929,0.951,0.971,
          0.989,1.004,1.018,1.029,1.038,1.044,1.049,1.052,1.054,1.054,
          1.053,1.051,1.049,1.047,1.044,1.041,1.039,1.036,1.034,1.031,
	      1.029,1.028};

      double ZF[81] =
	           {0.269,0.272,0.276,0.279,0.282,0.286,0.289,0.293,0.297,
          0.301,0.305,0.309,0.314,0.318,0.323,0.328,0.333,0.339,0.344,
          0.350,0.356,0.363,0.370,0.377,0.384,0.392,0.400,0.409,0.418,
          0.427,0.438,0.448,0.460,0.471,0.484,0.498,0.512,0.527,0.543,
          0.560,0.578,0.597,0.618,0.639,0.662,0.687,0.713,0.740,0.769,
          0.800,0.832,0.866,0.901,0.938,0.976,1.016,1.057,1.098,1.140,
          1.183,1.227,1.270,1.313,1.356,1.398,1.439,1.480,1.519,1.558,
          1.595,1.632,1.667,1.701,1.735,1.767,1.799,1.829,1.859,1.889,
	      1.917,1.945};

      double ZFSD[81] = 
	           {0.203,0.205,0.207,0.209,0.211,0.214,0.216,0.219,0.222,
          0.224,0.227,0.230,0.233,0.236,0.239,0.243,0.246,0.250,0.253,
          0.257,0.261,0.265,0.269,0.273,0.278,0.283,0.288,0.293,0.298,
          0.303,0.309,0.314,0.320,0.327,0.333,0.340,0.346,0.353,0.361,
          0.368,0.375,0.383,0.390,0.398,0.405,0.413,0.420,0.427,0.433,
          0.440,0.445,0.450,0.454,0.457,0.459,0.460,0.460,0.458,0.455,
          0.451,0.445,0.438,0.431,0.422,0.412,0.402,0.392,0.381,0.370,
          0.360,0.349,0.339,0.330,0.321,0.312,0.304,0.297,0.290,0.284,
	      0.278,0.272};

  float h,x,delta;
  float c1,c2,e2,e4;
  int n;
  // Bayesian statistics tells us to modify I/sigma by subtracting off sigma/2S
  // where S is the mean intensity in the resolution shell
  h = I/sigma - 0.5*sigma/S;
  // reject as unphysical reflections for which I < -3.7 sigma, or h < -4.0
  if (I/sigma < -3.7 || h < -4.0 ) {
	  //printf("unphys: %f %f %f %f\n",I,sigma,S,h);
	  return(0);
  }
  else {
	  if (h < 4.0) {
		  // use look up table if -4.0 < h < 4.0
          x = 10.0*(h+4.0);
          n = int(x);
          delta = x-n;
		  // linear interpolation
          J = (1.0-delta)*ZJ[n] + delta*ZJ[n+1];
          sigJ = (1.0-delta)*ZJSD[n] + delta*ZJSD[n+1];
          F = (1.0-delta)*ZF[n] + delta*ZF[n+1];
          sigF = (1.0-delta)*ZFSD[n] + delta*ZFSD[n+1];
		  // look up table gives J/sigma, so multiply by sigma to get output intensity
		  J *= sigma;
		  sigJ *= sigma;
		  F *= sqrt(sigma);
		  sigF *= sqrt(sigma);
	  }
	  else {
		  // if h > 4.0 use asymptotic formulas in French and Wilson Appendix
		  e2 = 1.0/(h*h);
		  e4 = e2*e2;
		  c1 = (1.0 - 0.375*e2 - 87.0/128.0*e4)*sqrt(h);
          c2 = sqrt((15.0/32.0*e4 + 0.25*e2)*h);

		  J = h*sigma*(1.0 - 0.5*e2 - 0.75*e4);
		  sigJ = 2.0*sigma*c1*c2;
		  F = c1*sqrt(sigma);
		  sigF = c2*sqrt(sigma);
	  }
	  return(1);
  }
}


void tricart(Cell cell, Mat33<float>& transf)
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
		    HKL hkl = ih.hkl();
    transf(0,0) = cell.a_star()*s3*sw;
    transf(0,1) = 0.0;
    transf(0,2) = 0.0;
		    HKL twin;
    transf(1,1) = cell.b_star();
    transf(1,2) = cell.c_star()*c1;
    transf(2,0) = cell.a_star()*s3*cw;
    transf(2,1) = 0.0;
    transf(2,2) = cell.c_star()*s1;
	return;
}

void Htest( HKL_data<data32::I_sigI> isig, Mat33<int> twinop, int &itwin, bool debug )
{
	typedef clipper::HKL_data_base::HKL_reference_index HRI;
    double HT=0.0;
    double HT2=0.0;
    double NT=0.0;
	std::vector<int> pdf(50,0);
    Vec3<int> jhkl, jhkl2;
    for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
	    if ( !isig[ih].missing() && !ih.hkl_class().centric() ) {
		    HKL hkl = ih.hkl();
		    jhkl[0] = hkl.h();
		    jhkl[1] = hkl.k();
		    jhkl[2] = hkl.l();
		    HKL twin;
		    jhkl2 = twinop*jhkl;
		    twin.h() = jhkl2[0];
		    twin.k() = jhkl2[1];
		    twin.l() = jhkl2[2];
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
	printf("Applying the H test for twinning: (Yeates Acta A44 142 (1980))\n");
	if (alpha > 0.05) {
		printf("H test suggests data is twinned\n");
		printf("Twinning fraction = %5.2f\n",alpha);
		itwin = 1;
	}
	else {
		printf("No twinning detected for this twinning operator\n");
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

