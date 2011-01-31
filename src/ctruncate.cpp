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
#include "ccp4_fortran.h"
#include "intensity_target.h"  // contains additions to resol_targetfn.h
#include "intensity_scale.h"   // contains additions to sfscale.cpp, sfscale.h, function_object_bases.h
#include "alt_hkl_datatypes.h"

#include "cpsf_utils.h"

#include <mmdb/mmdb_tables.h>


using namespace clipper;

// replacement for Wilson/Truncate

int truncate(  HKL_data<data32::I_sigI> isig,   HKL_data<data32::I_sigI>& jsig,   HKL_data<data32::F_sigF>& fsig,
			   clipper::ResolutionFn Sigma, float scalef, CSym::CCP4SPG *spg1, int& nrej, bool debug);
int truncate(  HKL_data<data32::J_sigJ_ano> isig,   HKL_data<data32::J_sigJ_ano>& jsig,   HKL_data<data32::G_sigG_ano>& fsig,
			   clipper::ResolutionFn Sigma, float scalef, CSym::CCP4SPG *spg1, int& nrej, bool debug);
int truncate_acentric(float I, float sigma, float S, float& J, float& sigJ, float& F, float& sigF, int& nrej, bool debug);
int truncate_centric(float I, float sigma, float S, float& J, float& sigJ, float& F, float& sigF, int& nrej, bool debug);
void straight_line_fit(std::vector<float> x, std::vector<float> y, std::vector<float> w, int n, float &a, float &b, float &siga, float &sigb);
void tricart(Cell cell, Mat33<float>& transf);
void Htest( HKL_data<data32::I_sigI> isig, Mat33<int> twinop, int &itwin, int scalefac, String s, CCP4Program& prog, bool debug );
void MatrixToString( Mat33<int> twinoper, String &s );

extern "C" void FORTRAN_CALL ( YYY_CELL2TG, yyy_cell2tg,
	   ( clipper::Cell& cell, double& sc_tol, int& ng, int *uu_g, int *u_g, 
		 int& lc, int& nc, int& nc2, int *uu_c, double *sc_c, int& ivb, int& ierr ),
	   ( clipper::Cell& cell, double& sc_tol, int& ng, int *uu_g, int *u_g, 
		 int& lc, int& nc, int& nc2, int *uu_c, double *sc_c, int& ivb, int& ierr ),
	   ( clipper::Cell& cell, double& sc_tol, int& ng, int *uu_g, int *u_g, 
		 int& lc, int& nc, int& nc2, int *uu_c, double *sc_c, int& ivb, int& ierr ));


int main(int argc, char **argv)
{
  CCP4Program prog( "ctruncate", "1.0.12", "$Date: 2009/12/11" );
  
  // defaults
  clipper::String outfile = "ctruncate_out.mtz";
  clipper::String outcol = "";
  clipper::String appendcol = "";
  clipper::String meancol = "/*/*/[IMEAN,SIGIMEAN]";
  //clipper::String meancol = "NONE";
  clipper::String pluscol = "/*/*/[I(+),SIGI(+)]";
  clipper::String minuscol = "/*/*/[I(-),SIGI(-)]";
  clipper::String anocols = "/*/*/[I(+),SIGI(+),I(-),SIGI(-)]";
  clipper::String ipfile = "NONE";
  clipper::String twintest = "first_principles";
  clipper::String ipseq = "NONE";

  bool aniso = true;
  bool debug = false;
  bool amplitudes = false;
  bool anomalous = false;

  int mtzinarg = 0;
  int mtzoutarg = 0;
  int nbins = 60;
  int ncbins = 60;
  int nresidues = 0;

  clipper::Resolution reso_Patt = clipper::Resolution( 4.0 );

  // clipper seems to use its own column labels, then append yours

  CCP4MTZfile mtzfile, mtzout;
  HKL_info hklinf, hklp;

  // command input
  prog.summary_beg();
  printf("\nUSER SUPPLIED INPUT:\n");
  CCP4CommandInput args( argc, argv, true ); 
  prog.summary_end();
  int arg = 0;
  while ( ++arg < args.size() ) {
    if ( args[arg] == "-mtzin" || args[arg] == "-hklin") {
		if ( ++arg < args.size() ) {
			ipfile = args[arg];
			mtzinarg = arg;
		}
    } else if ( args[arg] == "-mtzout" || args[arg] == "-hklout") {
		if ( ++arg < args.size() ) {
			outfile = args[arg];
			mtzoutarg = arg;
		}
    } else if ( args[arg] == "-colin" ) {
      if ( ++arg < args.size() ) meancol = args[arg];
    } else if ( args[arg] == "-colplus" ) {
      if ( ++arg < args.size() ) pluscol = args[arg];
	  anomalous = true;
	  printf("obsolete argument - use -colano instead\n");
	  printf("e.g. -colano /*/*/[I(+),SIGI(+),I(-),SIGI(-)]\n");
      return(0);
    } else if ( args[arg] == "-colminus" ) {
      if ( ++arg < args.size() ) minuscol = args[arg];
	  anomalous = true;
	  printf("obsolete argument - use -colano instead\n");
	  printf("e.g. -colano /*/*/[I(+),SIGI(+),I(-),SIGI(-)]\n");
      return(0);
    } else if ( args[arg] == "-colano" ) {
      if ( ++arg < args.size() ) anocols = args[arg];
	  anomalous = true;
    } else if ( args[arg] == "-colout" ) {
      if ( ++arg < args.size() ) appendcol = args[arg];
    } else if ( args[arg] == "-nbins" ) {
		if ( ++arg < args.size() ) nbins = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-nres" ) {
		if ( ++arg < args.size() ) nresidues = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-twintest" ) {
      if ( ++arg < args.size() ) twintest = args[arg];
    } else if ( args[arg] == "-seqin" ) {
      if ( ++arg < args.size() ) ipseq = args[arg];
    } else if ( args[arg] == "-tNCSres" ) {
      if ( ++arg < args.size() ) reso_Patt = clipper::Resolution( clipper::String(args[arg]).f() );
    } else if ( args[arg] == "-no-aniso" ) {
      aniso = false;
    } else if ( args[arg] == "-amplitudes" ) {
      amplitudes = true;
	} else if ( args[arg] == "-debug" ) {
      debug = true;
	} else {
	  printf("Unrecognised argument\n");
	  return(0);
	}

  }
  if (anomalous) {
    clipper::CCP4MTZ_type_registry::add_group( "G_sigG_ano", "FANO" );
    clipper::CCP4MTZ_type_registry::add_group( "J_sigJ_ano", "IANO" );
  }
  if ( args.size() <= 1 ) {
	  CCP4::ccperror(1,"Usage: ctruncate -mtzin <filename>  -mtzout <filename>  -colin <colpath> -colano <colpath> ");
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
  HKL_data<data32::J_sigJ_ano> isig_ano(hklinf);   // raw I(+) and sigma and I(-) and sigma
  HKL_data<data32::J_sigJ_ano> jsig_ano(hklinf);   // post-truncate anomalous I and sigma
  HKL_data<data32::G_sigG_ano> fsig_ano(hklinf);   // post-truncate anomalous F and sigma 
  HKL_data<data32::D_sigD> Dano(hklinf);   // anomalous difference and sigma 
  HKL_data<data32::I_sigI> ianiso(hklinf);   // anisotropy corrected I and sigma

//  clipper::MTZcrystal cxtl;
//  mtzfile.import_crystal( cxtl, meancol );
//  clipper::HKL_data<clipper::data32::F_sigF> faniso( hklinf, cxtl );  // don't seem to need crystal info
  clipper::HKL_data<clipper::data32::F_sigF> faniso( hklinf );

  if (amplitudes) {
      mtzfile.import_hkl_data( fsig, meancol );
  }
  else {
      //meancol = "/*/*/[" + meancol + ",SIG" + meancol + "]";
      mtzfile.import_hkl_data( isig, meancol );

      if (anomalous) {
          mtzfile.import_hkl_data( isig_ano, anocols );
          //pluscol = "/*/*/[" + pluscol + ",SIG" + pluscol + "]";
	      //mtzfile.import_hkl_data( isig_plus, pluscol );
          //minuscol = "/*/*/[" + minuscol + ",SIG" + minuscol + "]";
	      //mtzfile.import_hkl_data( isig_minus, minuscol );
      }
  }

  prog.summary_beg();
  printf("\n\nCRYSTAL INFO:\n\n");
  std::cout << "Crystal/dataset names: " << mtzfile.assigned_paths()[0].notail() << "\n"; 
  /*std::cout << mtzfile.assigned_paths()[0].notail().notail().tail() << "\n";  // crystal name
  std::cout << mtzfile.assigned_paths()[0].notail().tail() << "\n";  //dataset name
  String xtlname = mtzfile.assigned_paths()[0].notail().notail().tail();
  String setname = mtzfile.assigned_paths()[0].notail().tail();*/
  prog.summary_end();
  // need this mumbo-jumbo in order to write to output file

  if ( outcol[0] != '/' ) outcol = mtzfile.assigned_paths()[0].notail()+"/"+outcol;

  // hkl_list contains only those (h,k,l) for which at least data column is not NaN.
  // hkl_info contains all (h,k,l) out to the resolution limit regardless of whether there is any measured data.
  // will need hkl_list when writing to output file.
  clipper::Spacegroup spgr = mtzfile.spacegroup();
  clipper::Cell      cell1 = mtzfile.cell();
  clipper::Resolution reso = mtzfile.resolution();

  // limit resolution of Patterson calculation for tNCS (default 4 A)
  reso_Patt = clipper::Resolution( clipper::Util::max( reso.limit(), reso_Patt.limit() ) );

  HKL_info hkl_list;
  hkl_list.init( spgr, cell1, reso );
  mtzfile.import_hkl_list(hkl_list);

  MTZdataset cset; 
  MTZcrystal cxtl; 
  mtzfile.import_crystal ( cxtl, meancol );
  mtzfile.import_dataset ( cset, meancol );

  mtzfile.close_read();

  if (amplitudes) {
	  for ( HRI ih = fsig.first(); !ih.last(); ih.next() ) {
		  if ( !fsig[ih].missing() )
		  //isig[ih] = datatypes::I_sigI<float>( fsig[ih].f()*fsig[ih].f(), 2.0*fsig[ih].f()*fsig[ih].sigf() );
		  isig[ih] = clipper::data32::I_sigI( fsig[ih].f()*fsig[ih].f(), 2.0*fsig[ih].f()*fsig[ih].sigf() );
		  //printf("%f %f \n",fsig[ih].f(), isig[ih].I() );
      }
  }

  int Ncentric = 0;
  int Nreflections = 0;
  for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
    if ( !isig[ih].missing() ) Nreflections++;
	if ( ih.hkl_class().centric() && !isig[ih].missing()) Ncentric++;
  }
  printf("\nNcentric = %d\n", Ncentric);
  ncbins = std::min( Ncentric/10, nbins);
  printf("Number of centric bins = %d\n", ncbins);

  prog.summary_beg();

  printf("Cell parameters: %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n", cell1.a(), cell1.b(), cell1.c(), 
	      Util::rad2d( cell1.alpha() ), Util::rad2d( cell1.beta() ), Util::rad2d( cell1.gamma() ) );
  printf("\nNumber of reflections: %d\n", Nreflections);
  clipper::Grid_sampling grid;
  //clipper::String opfile = "patterson.map";

  // can't seem to get max resolution from clipper, so use CCP4 methods
  CMtz::MTZ *mtz1=NULL;
  int read_refs=1;  // not sure what read_refs actually does - reads reflections presumably
  float minres,maxres;
  mtz1 = CMtz::MtzGet(argv[mtzinarg], read_refs);

  // read title
  char title[72];
  CMtz::ccp4_lrtitl(mtz1, title);

  CMtz::MtzResLimits(mtz1,&minres,&maxres);
  float invopt = maxres;
  float resopt = 1.0/sqrt(invopt);
  printf("\nMinimum resolution = %7.3f A\nMaximum resolution = %7.3f A\n",1.0/sqrt(minres),1.0/sqrt(maxres));
  if (debug) printf("Minimum resolution = %f \nMaximum resolution = %f \n\n",minres,maxres);
  prog.summary_end();
  CSym::CCP4SPG *spg1 = CSym::ccp4spg_load_by_ccp4_num(CMtz::MtzSpacegroupNumber(mtz1));
  prog.summary_beg();

  // Clipper changes H3 to R3 so print out old spacegroup symbol instead
  //std::cout << "\nSpacegroup: " << spgr.symbol_hm() << " (number " << spgr.descr().spacegroup_number() << ")" << std::endl;
  
  char spacegroup[20];
  strcpy(spacegroup,spg1->symbol_old);
  printf("\nSpacegroup: %s (number %4d)\n", spg1->symbol_old, spg1->spg_ccp4_num);

  char pointgroup[20];
  strcpy(pointgroup,spg1->point_group);
  printf("Pointgroup: %s\n\n",pointgroup);
  prog.summary_end();

  // check for pseudo translation (taken from cpatterson)
  // get Patterson spacegroup
  clipper::Spacegroup
    pspgr( clipper::Spgr_descr( spgr.generator_ops().patterson_ops() ) );
  hklp.init( pspgr, cell1, reso_Patt, true );

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
  if ( grid.is_null() ) grid.init( pspgr, cell1, reso_Patt );

  // make xmap
  clipper::Xmap<float> patterson( pspgr, cell1, grid );
  patterson.fft_from( fphi );


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
  prog.summary_beg();
  printf("\n\nTRANSLATIONAL NCS:\n\n");
  if ( debug || (ratio > 0.2 && dist2 > 0.01) ) { 
	  printf("Translational NCS has been detected (with resolution limited to %5.2f A)\n", reso_Patt.limit() );
      printf("Ratio = %f\n",ratio);
      printf("Vector = (%6.3f, %6.3f, %6.3f)\n",c0[0],c0[1],c0[2]);
  }
  else {
	  printf("No translational NCS detected (with resolution limited to %5.2f A)\n", reso_Patt.limit() );
  }
  prog.summary_end();


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
	  double Itotal = 0.0;
      double FFtotal = 0.0;
	  prog.summary_beg();
	  printf("\n\nANISOTROPY CORRECTION (using intensities):\n");

/*
      clipper::HKL_data<clipper::data32::F_sigF> faniso( hklinf );
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
*/
	  for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {  
	  	  double I = isig[ih].I();
	      double sigI = isig[ih].sigI();
	      if ( I > 0.0 ) Itotal += I;
     	  ianiso[ih] = clipper::data32::I_sigI( I, sigI );
	  }

      // scale intensities 
	  clipper::Iscale_aniso<float> sfscl( 3.0 );     
	  sfscl( ianiso );                                        

      

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

	  // Eigenvalue calculation
	  Matrix<> mat( 3, 3, 0.0 );
	  for (int i=0; i<3; i++) {
		  for (int j=0; j<3; j++) {
              mat(i,j) = sfscl.u_aniso_orth()(i,j);
		  }
	  }
	  // add the isotropic part
	  // isotropic part is NOT diag(1,1,1) - would need to calculate properly
	  // should use eigenvalue differences rather then ratios

	  /*for (int i=0; i<3; i++) {
          mat(i,i) += 1.0;
	  }
      std::vector<ftype> v = mat.eigen( true );
	  printf("\nEigenvalues: %8.4f %8.4f %8.4f\n", v[0],v[1],v[2]);
	  printf("Eigenvalue ratios: %8.4f %8.4f %8.4f\n", v[0]/v[2], v[1]/v[2], v[2]/v[2]); 
	  if ( v[0] <= 0.0 ) CCP4::ccperror(1, "Anisotropy correction failed - negative eigenvalue.");
          invopt = maxres*v[0]/v[2];
	  resopt = 1.0/sqrt(invopt);
	  printf("Resolution limit in weakest direction = %7.3f A\n",resopt);
	  if ( v[0]/v[2] < 0.5 ) printf("\nWARNING! WARNING! WARNING! Your data is severely anisotropic\n");
	  */

	  prog.summary_end();
	  printf("\n");

/*
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
*/

      for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {    
		  if ( !isig[ih].missing() ) {
		      FFtotal += ianiso[ih].I();
		  }
	  }                                                         
	  
	  double scalefac = Itotal/FFtotal;
	  if (debug) printf("\nscalefactor = %6.3f %8.3f %8.3f\n\n",scalefac,Itotal,FFtotal);
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

  // truncate anisotropically corrected data at resolution limit in weakest direction
  for ( HRI ih = ianiso.first(); !ih.last(); ih.next() ) {
	  if (ih.invresolsq() > invopt) ianiso[ih].set_null();  
  }


  // calculate Sigma (mean intensity in resolution shell) 
  // use intensities uncorrected for anisotropy
  const int nprm = nbins;

  std::vector<double> params_init( nprm, 1.0 );
  clipper::BasisFn_spline basis_fo( isig, nprm );
  TargetFn_meanInth<clipper::data32::I_sigI> target_fo( isig, 1 );
  clipper::ResolutionFn Sigma( hklinf, basis_fo, target_fo, params_init );


  
  if (debug) {
      FILE *ftestfile;
      ftestfile = fopen("sigma.txt","w");
      for (int i=0; i!=nprm; ++i) {
          double res = maxres * pow( double(i+1)/nprm, 0.666666 );
	      fprintf(ftestfile,"%10.6f %10.6f \n", res, basis_fo.f_s( res, Sigma.params() ));
      }
      fclose(ftestfile);
  }

  // calculate moments of Z using truncate methods


  std::vector<double> Na(nbins,0.0),  Nc(ncbins,0.0);
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
		  else if (Ncentric) {
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

  printf("$TABLE: Acentric moments of E using Truncate method:\n");
  printf("$GRAPHS");
  printf(": 4th moment of E (Expected value = 2, Perfect Twin = 1.5):0|%5.3fx0|5:1,4:\n", maxres);
  printf(": 1st & 3rd moments of E (Expected values = 0.886, 1.329, Perfect twin = 0.94, 1.175):0|%5.3fx0|2:1,2,3:\n", maxres);
  //printf(": 6th & 8th moments of E (Expected value = 6, 24, Perfect Twin 3, 7.5):0|%5.3fx0|48:1,5,6:\n", maxres);
  printf("$$ 1/resol^2 <E> <E**3> <E**4> <E**6> <E**8> $$\n$$\n");

  for (int i=0; i<nbins; i++) {
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

  if (Ncentric) {
      printf("$TABLE: Centric moments of E using Truncate method:\n");
      printf("$GRAPHS");
      printf(": 4th moment of E (Expected = 3, Perfect Twin = 2):0|%5.3fx0|5:1,4:\n", maxres);
      printf(": 1st & 3rd moments of E (Expected = 0.798, 1.596, Perfect Twin = 0.886, 1.329):0|%5.3fx0|4:1,2,3:\n", maxres);
      //printf(": 6th & 8th moments of E (Expected = 15, 105, Perfect Twin = 6, 24):0|%5.3fx0|120:1,5,6:\n", maxres);
      printf("$$ 1/resol^2 <E> <E**3> <E**4> <E**6> <E**8> $$\n$$\n");

      for (int i=0; i<ncbins; i++) {
	      double res = maxres * pow((double(i) + 0.5)/double(ncbins), 0.666666);
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
  }

  prog.summary_beg();
  printf("\n\nTWINNING ANALYSIS:\n\n");
  int itwin = 0;
  int scalefac;
  Cell cell = hklinf.cell();

  // H test for twinning

  if (twintest != "table") {

      scalefac = 12;           // scale factor for integer symops
      double sc_tol = 0.05;    // tolerance for determining candidate twinops
      int lc = 48;             // maximum number of twinops that can be stored
      int nc = 0;
      int ivb = 0;
      int ierr = 0;
      int nc2;

      std::vector< Vec3<int> > trans;
      std::vector< Mat33<int> > rot;
      std::vector< Mat33<int> > twin(lc);   // NB need to dimension these o/w output from fortran corrupted
      std::vector<double> score(lc);

      int nsymops = spgr.num_symops();
      //printf("nsymops = %d\n",nsymops);
      Grid g( 12, 12, 12 );
      for (int i=0; i<nsymops; i++) {
	      Isymop isymop( spgr.symop(i), g );
	      /*for (int j=0; j<3; j++) {
		      printf("%6d %6d %6d   %6d\n", isymop.rot()(j,0), isymop.rot()(j,1), isymop.rot()(j,2), isymop.trn()[j] );
	      }
	      printf("\n");*/
	      trans.push_back( isymop.trn() );
	      // need to transpose matrix before passing to fortran
	      rot.push_back( isymop.rot().transpose() );
      }

      int *vv = &trans[0][0];
      int *ww = &rot[0](0,0);
      int *uu_c = &twin[0](0,0);
      double *sc_c = &score[0];


      FORTRAN_CALL ( YYY_CELL2TG, yyy_cell2tg,
                 (cell1, sc_tol, nsymops, ww, vv, lc, nc, nc2, uu_c, sc_c, ivb, ierr),
				 (cell1, sc_tol, nsymops, ww, vv, lc, nc, nc2, uu_c, sc_c, ivb, ierr),
				 (cell1, sc_tol, nsymops, ww, vv, lc, nc, nc2, uu_c, sc_c, ivb, ierr));

      //printf("nc = %d  ierr = %d\n",nc,ierr);
	  printf("First principles calculation of potential twinning operators using code by Andrey Lebedev:\n");
	  printf("First principles calculation has found %d potential twinning operators\n\n", nc-1);


      if (nc > 1) {
	      for (int k=1; k<nc; k++) {
			  // Transpose matrix once because passed from Fortran to C, then a second time because Andrey's
			  // convention is that (h,k,l) is postmultiplied by the twinning operator, whereas mine is that
			  // (h,k,l) is premultiplied by the the twinning operator. Net result: don't do anything!
              Mat33<int> twinoper = twin[k];
			  String s;
			  MatrixToString(twinoper,s);
			  std::cout << "Twinning operator: " << s << "\n";
	          /*for (int i=0; i<3; i++) {
		          // Divide by 12 (scale factor for integer syops)
		          printf("%7.4f %7.4f %7.4f\n",double(twinoper(i,0))/12.0, double(twinoper(i,1))/12.0, double(twinoper(i,2))/12.0 );
		      }
			  printf("\n");*/
	          Htest(ianiso, twinoper, itwin, scalefac, s, prog, debug);
	      }
      }
  }


  if (twintest != "first_principles") {
	  printf("\n   Potential twinning operators found from tables:\n\n");
      Mat33<int> twinop(0,0,0,0,0,0,0,0,0);
      int sg = CMtz::MtzSpacegroupNumber(mtz1);
      scalefac = 1;
	  String s;
      if ( (sg >= 75 && sg <= 80) || sg == 146 || (sg >= 168 && sg <= 173) || (sg >= 195 && sg <= 199) ) { 
	      printf("Twinning operator k, h, -l\n");
		  s = "k, h, -l";
	      twinop(0,1) = 1;
	      twinop(1,0) = 1;
	      twinop(2,2) = -1;
	      Htest ( ianiso, twinop, itwin, scalefac, s, prog, debug );
      }
      else if( sg >= 149 && sg <= 154 ) {
	      printf("Twinning operator -h, -k, l\n");
		  s = "-h, -k, l";
	      twinop(0,0) = -1;
	      twinop(1,1) = -1;
	      twinop(2,2) = 1;
	      Htest ( ianiso, twinop, itwin, scalefac, s, prog, debug );
      }
      else if( sg >= 143 && sg <= 145 ) {
	      printf("Twinning operator k, h, -l\n");
		  s = "k, h, -l";
	      twinop(0,1) = 1;
	      twinop(1,0) = 1;
	      twinop(2,2) = -1;
	      Htest( ianiso, twinop, itwin, scalefac, s, prog, debug );

	      printf("Twinning operator -k, -h, -l\n");
		  s = "-k, -h, -l";
	      twinop(0,1) = -1;
	      twinop(1,0) = -1;
	      Htest( ianiso, twinop, itwin, scalefac, s, prog, debug );

	      printf("Twinning operator -h, -k, l\n");
		  s = "-h, -k, l";
          twinop(0,1) = 0;
	      twinop(1,0) = 0;
	      twinop(0,0) = -1;
	      twinop(1,1) = -1;
	      twinop(2,2) = 1;
	      Htest( ianiso, twinop, itwin, scalefac, s, prog, debug );
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
	          Htest( ianiso, twinop, itwin, scalefac, s, prog, debug );
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
	          Htest( ianiso, twinop, itwin, scalefac, s, prog, debug );
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
	          Htest( ianiso, twinop, itwin, scalefac, s, prog, debug );
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
	          Htest( ianiso, twinop, itwin, scalefac, s, prog, debug );
	      }

	      if ( cell.beta() < Util::d2rad(93.0) ) {
		      printf("Twinning operator -h, -k, l\n");
			  s = "-h, -k, l";
		      twinop(0,0) = -1;
		      twinop(0,2) = 0;
		      twinop(1,1) = -1;
		      twinop(2,0) = 0;
		      twinop(2,2) = 1;
		      Htest( ianiso, twinop, itwin, scalefac, s, prog, debug );
	      }	 

	      if ( fabs( cos(cell.beta()) + 0.5*cell.a()/cell.c() ) < 0.02 ) { 
		      printf("Twinning operator -h, -k, h+l\n");
			  s = "-h, -k, h+l";
		      twinop(0,0) = -1;
		      twinop(0,2) = 0;
		      twinop(1,1) = -1;
		      twinop(2,0) = 1;
		      twinop(2,2) = 1;
		      Htest( ianiso, twinop, itwin, scalefac, s, prog, debug );
	      }	  

	      if ( fabs( cos(cell.beta()) + 0.5*cell.c()/cell.a() ) < 0.02 ) { 
		      printf("Twinning operator h+l, -k, -l\n");
			  s = "h+l, -k, -l";
		      twinop(0,0) = 1;
		      twinop(0,2) = 1;
		      twinop(1,1) = -1;
		      twinop(2,0) = 0;
		      twinop(2,2) = -1;
		      Htest( ianiso, twinop, itwin, scalefac, s, prog, debug );
	      }	  
	      if ( fabs( cos(cell.beta()) + cell.c()/cell.a() ) < 0.02 ) { 
		      printf("Twinning operator h+2l, -k, -l\n");
			  s = "h+2l, -k, -l";
		      twinop(0,0) = 1;
		      twinop(0,2) = 2;
		      twinop(1,1) = -1;
		      twinop(2,0) = 0;
		      twinop(2,2) = -1;
		      Htest( ianiso, twinop, itwin, scalefac, s, prog, debug );
	      }
      }
  }
  if (itwin) {
	  //printf("\nData has been truncated at %6.2f A resolution\n",resopt);
  	  if (aniso) printf("Anisotropy correction has been applied before calculating H\n\n");
	  else printf("Anisotropy correction has not been applied before calculating H\n\n");
  }

  // L test for twinning

  double LT=0.0;
  double LT2=0.0;
  double NLT=0.0;
  std::vector<int> cdf(20,0);

  for ( HRI ih = ianiso.first(); !ih.last(); ih.next() ) {
	  if ( !ianiso[ih].missing() && !ih.hkl_class().centric() ) {
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
  				          double I1 = ianiso[ih].I();
		                  double I2 = ianiso[hkl2].I();
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
	  printf("\nApplying the L test for twinning: (Padilla and Yeates Acta Cryst. D59 1124 (2003))\n");
	  printf("L test suggests data is twinned\n");
	  printf("L statistic = %6.3f  (untwinned 0.5 perfect twin 0.375)\n\n", Lav);
	  itwin = 1;
	  printf("All data regardless of I/sigma(I) has been included in the L test\n");
	  //printf("Data has been truncated at %6.2f A resolution\n",resopt);
	  if (aniso) printf("Anisotropy correction has been applied before calculating L\n\n");
	  else printf("Anisotropy correction has not been applied before calculating L\n\n");
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


  //Wilson pre
  
  prog.summary_beg();
  printf("\nWILSON SCALING:\n\n");
  int nsym = spg1->nsymop;
  clipper::MMoleculeSequence seq;

  std::string name[5] = { "C", "N", "O", "H", "S" };
  int numatoms[5] = { 0, 0, 0, 0, 0 }; 

// Use single letter residue names from mmdb - note that C appears twice

//                   A   R   N   D   C   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
  int Catoms[21] = { 3,  6,  4,  4,  3,  3,  5,  5,  2,  6,  6,  6,  6,  5,  9,  5,  3,  4, 11,  9,  5 };
  int Hatoms[21] = { 7, 14,  8,  7,  7,  7, 10,  9,  5,  9, 13, 13, 14, 11, 11,  9,  7,  9, 12, 11, 11 };
  int Natoms[21] = { 1,  4,  2,  1,  1,  1,  2,  1,  1,  3,  1,  1,  2,  1,  1,  1,  1,  1,  2,  1,  1 };
  int Oatoms[21] = { 2,  2,  3,  4,  2,  2,  3,  4,  2,  2,  2,  2,  2,  2,  2,  2,  3,  3,  2,  3,  2 };
  int Satoms[21] = { 0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0 };

  if ( ipseq != "NONE" ) {
    clipper::SEQfile seqf;
    seqf.read_file( ipseq );
    seqf.import_molecule_sequence( seq );
	MPolymerSequence poly = seq[0];
	String sequence = poly.sequence();
	//std::cout << poly.sequence() << "\n";
	for (int i=0; i<sequence.length(); i++) {
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
	printf("User supplied sequence contains %d C, %d N, %d O, %d H, %d S atoms\n", numatoms[0],numatoms[1],numatoms[2],numatoms[3],numatoms[4]);
  }

  else if (nresidues > 0) {
	  printf("User supplied number of residues = %d\n",nresidues);
  }
  else {
      nresidues = int(0.5*hklinf.cell().volume()/(nsym*157));
      printf("Estimated number of residues = %d\n",nresidues);
  }
  prog.summary_end();

  if ( ipseq == "NONE" ) {
      numatoms[0] = 5*nresidues;
      numatoms[1] = int(1.35*nresidues);
      numatoms[2] = int(1.5*nresidues);
      numatoms[3] = 8*nresidues;
  }


  // Wilson plot
  
  std::vector<float> xi, yi, wi, xxi, yyi, yy2i; 
  float totalscat; 
  float minres_scaling = 0.0625;   // 4 Angstroms

  for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
	  if ( !isig[ih].missing() && Sigma.f(ih) > 0.0) {
		  float lnS = -log(Sigma.f(ih));
		  float res = ih.invresolsq();

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
		  lnS += log(totalscat);
		  
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

  int nobs = xi.size();
  //printf("%d %d %d\n", xi.size(), yi.size(), wi.size());
  float a,b,siga,sigb,a1,b1;
  b = 0.0;

  if ( wi.size() > 200 && maxres > 0.0816) {               // 3.5 Angstroms
      straight_line_fit(xi,yi,wi,nobs,a,b,siga,sigb);
	  prog.summary_beg();
	  printf("\nResults from Clipper style Wilson plot:\n");
      a *= 2.0;
      printf ("B = %6.3f intercept = %6.3f siga = %6.3f sigb = %6.3f\n",a,b,siga,sigb);
      printf("scale factor on intensity = %10.4f\n\n",(exp(b)));
	  prog.summary_end();
  }
  else {
	  printf("Too few high resolution points to determine B factor and Wilson scale factor\n");
  }
  
  printf("$TABLE: Wilson plot:\n");
  printf("$GRAPHS");
  //printf(": Wilson plot:0|0.1111x-7|-5:1,2:\n$$");  // limits hardwired
  printf(": Wilson plot - estimated B factor = %5.1f :A:1,2:\n$$", a);  
  printf(" 1/resol^2 ln(I/I_th) $$\n$$\n");

  for ( int i=0; i!=nbins; ++i ) {
		float res = maxres*(float(i)+0.5)/float(nbins); 
		float totalscat = 0;
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
	  printf("%10.5f %10.5f\n", res,log(basis_fo.f_s( res, Sigma.params() ))-log(totalscat));
  }

  printf("$$\n\n");



  // Truncate style Wilson plot

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

 
  // apply the Truncate procedure, unless amplitudes have been input

  // if something went wrong with Wilson scaling, B could be negative, giving exponentially large scaled SF's
  // so only scale if B positive
  float scalef = 1.0;
  if ( b > 0 ) scalef = sqrt(exp(b));
  int nrej = 0; 

  if (!amplitudes) {
      if (anomalous) {
	      truncate( isig_ano, jsig_ano, fsig_ano, Sigma, scalef, spg1, nrej, debug );
	      int iwarn = 0;
	      for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			  if ( !Util::is_nan(fsig_ano[ih].f_pl() )  &&  !Util::is_nan(fsig_ano[ih].f_mi() ) ) {
			      fsig[ih].f() = 0.5 * ( fsig_ano[ih].f_pl() + fsig_ano[ih].f_mi() );
			      fsig[ih].sigf() = 0.5 * sqrt( pow( fsig_ano[ih].sigf_pl(), 2 ) + pow( fsig_ano[ih].sigf_mi(), 2 ) );
			      Dano[ih].d() = fsig_ano[ih].f_pl() - fsig_ano[ih].f_mi();
			      Dano[ih].sigd() = 2.0 * fsig[ih].sigf();
		      }
		      else if ( !Util::is_nan(fsig_ano[ih].f_pl() ) ) {
			      fsig[ih].f() = fsig_ano[ih].f_pl();
			      fsig[ih].sigf() = fsig_ano[ih].sigf_pl();
		      }
		      else if ( !Util::is_nan(fsig_ano[ih].f_mi() ) ) {
			      fsig[ih].f() = fsig_ano[ih].f_mi();
			      fsig[ih].sigf() = fsig_ano[ih].sigf_mi();
		      }
		      else if ( !isig[ih].missing() && iwarn != 1 ) {
			      printf("\nWARNING: Imean exists but I(+), I(-) do not\n\n");
			      iwarn = 1;
		      }
		      if ( ih.hkl_class().centric() ) {
			      Dano[ih].d() = 0.0;
		  HKL hkl = ih.hkl();
		      }
	      }
      }

      else {
          truncate( isig, jsig, fsig, Sigma, scalef, spg1, nrej, debug );
      }
  }
  prog.summary_beg();
  printf("\nINTENSITY TO AMPLITUDE CONVERSION:\n\n");
  printf("%d intensities have been rejected as unphysical\n", nrej);
  prog.summary_end();


  
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


  // moments of E using clipper binning

  //acentrics
  clipper::HKL_data<clipper::data32::F_sigF> fs1( hklinf );
  for ( HRI ih = ianiso.first(); !ih.last(); ih.next() ) {
	  if ( !ianiso[ih].missing() && !ih.hkl_class().centric() ) {
	      double I = ianiso[ih].I();
	      double sigI = ianiso[ih].sigI();
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
  printf(": 4th moment of E (Expected value = 2, Perfect Twin = 1.5):0|%5.3fx0|5:1,4:\n", invopt);
  printf(": 1st & 3rd moments of E (Expected values = 0.886, 1.329, Perfect twin = 0.94, 1.175):0|%5.3fx0|2:1,2,3:\n", invopt);
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
  printf("$$\n\n");

  mean1 /= double(nbins);
  mean3 /= double(nbins);
  mean4 /= double(nbins);
  mean6 /= double(nbins);
  mean8 /= double(nbins);
 
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

  nprmk = ncbins;
  //centrics
  if (Ncentric) {
      clipper::HKL_data<clipper::data32::F_sigF> fs2( hklinf );
      for ( HRI ih = ianiso.first(); !ih.last(); ih.next() ) {
	      if ( !ianiso[ih].missing() && ih.hkl_class().centric() ) {
	          double I = ianiso[ih].I();
	          double sigI = ianiso[ih].sigI();
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
      printf(": 4th moment of E (Expected = 3, Perfect Twin = 2):0|%5.3fx0|5:1,4:\n", invopt);
      printf(": 1st & 3rd moments of E (Expected = 0.798, 1.596, Perfect Twin = 0.886, 1.329):0|%5.3fx0|4:1,2,3:\n", invopt);
      //printf(": 6th & 8th moments of E (Expected = 15, 105, Perfect Twin = 6, 24):0|%5.3fx0|120:1,5,6:\n", maxres);

      //printf("$$ 1/resol^2 <E> <E**3> <E**4> <E**6> <E**8> $$\n$$\n");
      printf("$$ 1/resol^2 <E> <E**3> <E**4> $$\n$$\n");


      for (int i=0; i<ncbins; i++) {
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


  // construct cumulative distribution function for intensity (using Z rather than E)
  clipper::Range<double> intensity_range_centric;
  clipper::Range<double> intensity_range_acentric;
  // changed from jsig to isig

  for ( HRI ih = ianiso.first(); !ih.last(); ih.next() ) {
	  if ( !ianiso[ih].missing() ) {
         if ( ih.hkl_class().centric() ) intensity_range_centric.include( (ianiso[ih].I()/ih.hkl_class().epsilon() )/Sigma.f(ih) );
		 else intensity_range_acentric.include( (ianiso[ih].I()/ih.hkl_class().epsilon() ) /Sigma.f(ih) );
	  }
  }
  //printf("C2: %20.5f %20.5f\n",intensity_range_centric.max(),intensity_range_centric.min());
  //printf("A2: %20.5f %20.5f\n",intensity_range_acentric.max(),intensity_range_acentric.min());

  clipper::Generic_ordinal intensity_ord_c;
  clipper::Generic_ordinal intensity_ord_a;
  intensity_ord_c.init( intensity_range_centric );
  intensity_ord_a.init( intensity_range_acentric );
  for ( HRI ih = ianiso.first(); !ih.last(); ih.next() ) {
  	  if ( !ianiso[ih].missing() ) {
		  if ( ih.hkl_class().centric() ) intensity_ord_c.accumulate( ( ianiso[ih].I()/ih.hkl_class().epsilon() )/Sigma.f(ih)  );
		  else intensity_ord_a.accumulate( ( ianiso[ih].I()/ih.hkl_class().epsilon() )/Sigma.f(ih) );
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
	  if (Ncentric) printf("%10.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n", x, acen[i], acen_twin[i], intensity_ord_a.ordinal(x), cen[i], intensity_ord_c.ordinal(x));
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
  if (ntw > 2) {
	  prog.summary_beg();
	  printf("\nWARNING: ****  Cumulative Distribution shows Possible Twinning ****\n\n");
	  prog.summary_end();
  }


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

  int nzerosigma = 0;
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
					  if ( fsig[ih].sigf() > 0.0 ) somsddir[j][bin] += epsiln*fsig[ih].f()/fsig[ih].sigf();
                      enumdir[j][bin] += epsiln;
                      numdir[j][bin]++;
			      }
			  }
              somov[bin] += fsig[ih].f()*epsiln;
              if ( fsig[ih].sigf() > 0.0 ) somsdov[bin] += epsiln*fsig[ih].f()/fsig[ih].sigf();
			  else nzerosigma++;
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


  // output data
  if (!amplitudes) {
      //mtzout.open_append( argv[mtzinarg], outfile );
	  mtzout.open_write( outfile );
	  mtzout.export_crystal ( cxtl, outcol );
      mtzout.export_dataset ( cset, outcol );
      mtzout.export_hkl_info( hkl_list );
      //mtzout.export_hkl_data( jsig, outcol );
	  clipper::String labels;
	  if (appendcol == "") labels = outcol + "[F,SIGF]";
	  else labels = outcol + "[F_" + appendcol + ",SIGF_" + appendcol + "]";
	  mtzout.export_hkl_data( fsig, labels );
      if (anomalous) {
	      if (appendcol == "") labels = outcol + "[F(+),SIGF(+),F(-),SIGF(-)]";
	      else labels = outcol + "[F_" + appendcol + "(+),SIGF_" + appendcol + "(+),F_" + appendcol + "(-),SIGF_" + appendcol + "(-)]";
	      mtzout.export_hkl_data( fsig_ano, labels );
	      if (appendcol == "") labels = outcol + "[DANO,SIGDANO]";
	      else labels = outcol + "[DANO_" + appendcol + ",SIGDANO_" + appendcol + "]";
		  mtzout.export_hkl_data( Dano, labels );
      }
	  if (appendcol != "") {
		  String::size_type loc = meancol.find(",",0);
          meancol.insert(loc,"_"+appendcol);
		  loc = meancol.find("]",0);
		  meancol.insert(loc,"_"+appendcol);
	  }
	  mtzout.export_hkl_data( isig, outcol + meancol.tail() );

	  if (anomalous) {
		  if (appendcol != "") {
		      String::size_type loc = anocols.find("+",0);
              anocols.insert(loc-1,"_"+appendcol);
			  loc = anocols.find(",",0);
		      loc = anocols.find("+",loc+1);
		      anocols.insert(loc-1,"_"+appendcol);
			  loc = anocols.find("-",0);
              anocols.insert(loc-1,"_"+appendcol);
			  loc = anocols.find(",",loc);
		      loc = anocols.find("-",loc+1);
		      anocols.insert(loc-1,"_"+appendcol);
	      }
		  mtzout.export_hkl_data( isig_ano, outcol + anocols.tail() );
	  }

      //mtzout.close_append();
	  mtzout.close_write();

	  // Clipper will change H3 to R3, so change it back

	  CMtz::MTZ *mtz2=NULL;
      read_refs=1;  // need to read in reflections, otherwise they won't be written out
      mtz2 = CMtz::MtzGet(argv[mtzoutarg], read_refs);
	  // write title to output file
	  strncpy( mtz2->title, title, 71 );
	  if (spacegroup[0] == 'H') {
	      strcpy(mtz2->mtzsymm.spcgrpname,spacegroup);
	  }
		  HKL hkl = ih.hkl();
      CMtz::MtzFree( mtz2 );
  }
  CMtz::MtzFree( mtz1 );
  prog.set_termination_message( "Normal termination" );

  return(0);
}


int truncate(  HKL_data<data32::I_sigI> isig,   HKL_data<data32::I_sigI>& jsig,   HKL_data<data32::F_sigF>& fsig,
			   clipper::ResolutionFn Sigma, float scalef, CSym::CCP4SPG *spg1, int& nrej, bool debug)
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
		  if ( ih.hkl_class().centric() ) iflag = truncate_centric(I, sigma, S, J, sigJ, F, sigF, nrej, debug);
		  else iflag = truncate_acentric(I, sigma, S, J, sigJ, F, sigF, nrej, debug);	
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


int truncate(  HKL_data<data32::J_sigJ_ano> isig,   HKL_data<data32::J_sigJ_ano>& jsig,   HKL_data<data32::G_sigG_ano>& fsig,
			   clipper::ResolutionFn Sigma, float scalef, CSym::CCP4SPG *spg1, int& nrej, bool debug)

// takes anomalous I's as input. 

{
  typedef clipper::HKL_data_base::HKL_reference_index HRI;
  //FILE *checkfile;
  //checkfile = fopen("checku.txt", "w");
  float J, sigJ, F, sigF;
  int iflag;

  for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
	  if ( !Util::is_nan(isig[ih].I_pl() ) ) {
		  float I = isig[ih].I_pl();
		  float sigma = isig[ih].sigI_pl();
		  float S = Sigma.f(ih);
		  HKL hkl = ih.hkl();
		  float weight = (float) CSym::ccp4spg_get_multiplicity( spg1, hkl.h(), hkl.k(), hkl.l() );
		  if( fabs( ih.hkl_class().epsilon() - weight ) > 0.001) printf("epsilon %f != weight %f", ih.hkl_class().epsilon(), weight);
		  float sqwt = sqrt(weight);
		  HKL hkl = ih.hkl();
		  I /= weight;
		  sigma /= weight;

		  // handle acentric and centric reflections separately
		  if ( ih.hkl_class().centric() ) iflag = truncate_centric(I, sigma, S, J, sigJ, F, sigF, nrej, debug);
		  else iflag = truncate_acentric(I, sigma, S, J, sigJ, F, sigF, nrej, debug);	
		  //if ( !ih.hkl_class().centric()  && I < 0 ) 
			  //fprintf(checkfile,"%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %8.4f %8.4f %8.4f\n", I,sigma,S,J,sigJ,F,sigF,weight,
			  //ih.hkl_class().epsilon());
		  if (iflag) {
			  jsig[ih].I_pl() = J;
			  jsig[ih].sigI_pl() = sigJ;
			  //jsig[ih] = datatypes::I_sigI_ano<float>( J, sigJ );
			  fsig[ih].f_pl() = F*scalef*sqwt;
			  fsig[ih].sigf_pl() = sigF*scalef*sqwt;
			  //fprintf(checkfile,"%12.6f %12.6f %12.6f\n", I,fsig[ih].f(),fsig[ih].sigf());
		  }
	  }

	  if ( !Util::is_nan(isig[ih].I_mi() ) ) {
		  float I = isig[ih].I_mi();
		  float sigma = isig[ih].sigI_mi();
		  float S = Sigma.f(ih);
		  HKL hkl = ih.hkl();
		  float weight = (float) CSym::ccp4spg_get_multiplicity( spg1, hkl.h(), hkl.k(), hkl.l() );
		  if( fabs( ih.hkl_class().epsilon() - weight ) > 0.001) printf("epsilon %f != weight %f", ih.hkl_class().epsilon(), weight);
		  float sqwt = sqrt(weight);

		  I /= weight;
		  sigma /= weight;

		  // handle acentric and centric reflections separately
		  if ( ih.hkl_class().centric() ) iflag = truncate_centric(I, sigma, S, J, sigJ, F, sigF, nrej, debug);
		  else iflag = truncate_acentric(I, sigma, S, J, sigJ, F, sigF, nrej, debug);	
		  //if ( !ih.hkl_class().centric()  && I < 0 ) 
			  //fprintf(checkfile,"%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %8.4f %8.4f %8.4f\n", I,sigma,S,J,sigJ,F,sigF,weight,
			  //ih.hkl_class().epsilon());
		  if (iflag) {
			  jsig[ih].I_mi() = J;
			  jsig[ih].sigI_mi() = sigJ;
			  //jsig[ih] = datatypes::I_sigI_ano<float>( J, sigJ );
			  fsig[ih].f_mi() = F*scalef*sqwt;
			  fsig[ih].sigf_mi() = sigF*scalef*sqwt;
			  //fprintf(checkfile,"%12.6f %12.6f %12.6f\n", I,fsig[ih].f(),fsig[ih].sigf());
		  }
	  }
  }
  //fclose(checkfile);
  return(1);
}


int truncate_acentric(float I, float sigma, float S, float& J, float& sigJ, float& F, float& sigF, int& nrej, bool debug)
{
  // look up tables calculated using quadpack
  // tables give values from h = -4.0 to h = 3.0 in steps of 0.1
  // declare as doubles for now, to avoid compiler warnings  

  double ZJ[71] = {

  0.22561, 0.23037, 0.23531, 0.24046, 0.24581, 0.25139, 0.25720, 0.26327, 0.26959, 0.27620,
  0.28310, 0.29032, 0.29787, 0.30577, 0.31406, 0.32274, 0.33186, 0.34143, 0.35150, 0.36208,
  0.37322, 0.38495, 0.39731, 0.41036, 0.42413, 0.43868, 0.45406, 0.47033, 0.48755, 0.50580,
  0.52514, 0.54564, 0.56740, 0.59050, 0.61503, 0.64108, 0.66876, 0.69817, 0.72942, 0.76262,
  0.79788, 0.83533, 0.87507, 0.91722, 0.96188, 1.00916, 1.05915, 1.11192, 1.16756, 1.22611,
  1.28760, 1.35205, 1.41944, 1.48974, 1.56288, 1.63879, 1.71735, 1.79844, 1.88189, 1.96756,
  2.05525, 2.14478, 2.23597, 2.32863, 2.42258, 2.51764, 2.61365, 2.71046, 2.80794, 2.90596,
  3.00444  };

  double ZJSD[71] = {

  0.21604, 0.22024, 0.22459, 0.22910, 0.23377, 0.23861, 0.24362, 0.24882, 0.25422, 0.25982,
  0.26563, 0.27167, 0.27794, 0.28446, 0.29124, 0.29828, 0.30562, 0.31324, 0.32118, 0.32945,
  0.33805, 0.34701, 0.35634, 0.36606, 0.37618, 0.38671, 0.39768, 0.40910, 0.42099, 0.43335,
  0.44620, 0.45956, 0.47343, 0.48781, 0.50272, 0.51815, 0.53410, 0.55056, 0.56751, 0.58494,
  0.60281, 0.62109, 0.63974, 0.65869, 0.67789, 0.69726, 0.71673, 0.73619, 0.75555, 0.77470,
  0.79353, 0.81192, 0.82977, 0.84696, 0.86339, 0.87895, 0.89357, 0.90718, 0.91972, 0.93117,
  0.94152, 0.95076, 0.95894, 0.96609, 0.97226, 0.97755, 0.98200, 0.98573, 0.98880, 0.99130,
  0.99331  };

  double ZF[71] = {

  0.42331, 0.42784, 0.43250, 0.43731, 0.44226, 0.44737, 0.45264, 0.45808, 0.46369, 0.46949,
  0.47548, 0.48168, 0.48809, 0.49473, 0.50160, 0.50873, 0.51611, 0.52378, 0.53173, 0.53999,
  0.54857, 0.55749, 0.56677, 0.57643, 0.58649, 0.59696, 0.60788, 0.61927, 0.63115, 0.64354,
  0.65648, 0.66999, 0.68410, 0.69884, 0.71424, 0.73033, 0.74715, 0.76471, 0.78305, 0.80220,
  0.82218, 0.84301, 0.86472, 0.88732, 0.91082, 0.93522, 0.96052, 0.98672, 1.01379, 1.04171,
  1.07043, 1.09992, 1.13013, 1.16097, 1.19239, 1.22431, 1.25663, 1.28928, 1.32215, 1.35516,
  1.38822, 1.42125, 1.45416, 1.48689, 1.51936, 1.55152, 1.58334, 1.61476, 1.64576, 1.67632,
  1.70644  };

  double ZFSD[71] = {

  0.21545, 0.21753, 0.21966, 0.22185, 0.22409, 0.22639, 0.22874, 0.23116, 0.23363, 0.23618,
  0.23878, 0.24146, 0.24420, 0.24702, 0.24990, 0.25287, 0.25591, 0.25903, 0.26222, 0.26550,
  0.26886, 0.27230, 0.27583, 0.27944, 0.28313, 0.28690, 0.29075, 0.29468, 0.29868, 0.30274,
  0.30687, 0.31106, 0.31529, 0.31956, 0.32386, 0.32816, 0.33246, 0.33673, 0.34095, 0.34510,
  0.34915, 0.35307, 0.35683, 0.36039, 0.36372, 0.36678, 0.36951, 0.37190, 0.37389, 0.37544,
  0.37653, 0.37711, 0.37716, 0.37667, 0.37561, 0.37398, 0.37179, 0.36906, 0.36581, 0.36207,
  0.35789, 0.35332, 0.34841, 0.34323, 0.33783, 0.33228, 0.32663, 0.32095, 0.31529, 0.30968,
  0.30416  };

  float h,x,delta;
  int n;
  
  // Bayesian statistics tells us to modify I/sigma by subtracting off sigma/S
  // where S is the mean intensity in the resolution shell
  h = I/sigma - sigma/S;
  // reject as unphysical reflections for which I < -3.7 sigma, or h < -4.0
  if (I/sigma < -3.7 || h < -4.0 ) {
	  nrej++;
	  if (debug) printf("unphys: %f %f %f %f\n",I,sigma,S,h);
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


int truncate_centric(float I, float sigma, float S, float& J, float& sigJ, float& F, float& sigF, int& nrej, bool debug)
{
  // look up tables calculated using quadpack
  // tables give values from h = -4.0 to h = 4.0 in steps of 0.1

  double ZJ[91] = {

  0.11544, 0.11798, 0.12063, 0.12340, 0.12628, 0.12929, 0.13244, 0.13573, 0.13917, 0.14279,
  0.14657, 0.15055, 0.15473, 0.15912, 0.16374, 0.16862, 0.17376, 0.17918, 0.18492, 0.19099,
  0.19742, 0.20424, 0.21148, 0.21917, 0.22736, 0.23609, 0.24540, 0.25535, 0.26599, 0.27739,
  0.28960, 0.30272, 0.31681, 0.33198, 0.34833, 0.36596, 0.38499, 0.40556, 0.42781, 0.45190,
  0.47799, 0.50627, 0.53692, 0.57016, 0.60619, 0.64523, 0.68752, 0.73327, 0.78271, 0.83605,
  0.89346, 0.95513, 1.02117, 1.09167, 1.16666, 1.24610, 1.32989, 1.41787, 1.50981, 1.60539,
  1.70427, 1.80605, 1.91030, 2.01659, 2.12447, 2.23353, 2.34338, 2.45369, 2.56417, 2.67456,
  2.78470, 2.89443, 3.00367, 3.11236, 3.22047, 3.32800, 3.43496, 3.54139, 3.64732, 3.75278,
  3.85783, 3.96249, 4.06681, 4.17083, 4.27457, 4.37807, 4.48134, 4.58442, 4.68732, 4.79006,
  4.89265  };

  double ZJSD[91] = {

  0.15779, 0.16105, 0.16444, 0.16796, 0.17162, 0.17543, 0.17939, 0.18351, 0.18781, 0.19229,
  0.19697, 0.20185, 0.20694, 0.21227, 0.21784, 0.22367, 0.22977, 0.23617, 0.24287, 0.24990,
  0.25728, 0.26503, 0.27317, 0.28173, 0.29073, 0.30020, 0.31018, 0.32068, 0.33175, 0.34341,
  0.35571, 0.36867, 0.38233, 0.39673, 0.41191, 0.42790, 0.44473, 0.46244, 0.48106, 0.50060,
  0.52108, 0.54251, 0.56489, 0.58819, 0.61238, 0.63741, 0.66320, 0.68964, 0.71661, 0.74395,
  0.77148, 0.79898, 0.82620, 0.85289, 0.87877, 0.90354, 0.92694, 0.94869, 0.96857, 0.98639,
  1.00200, 1.01532, 1.02636, 1.03515, 1.04181, 1.04651, 1.04946, 1.05089, 1.05106, 1.05021,
  1.04859, 1.04642, 1.04389, 1.04115, 1.03835, 1.03558, 1.03291, 1.03039, 1.02805, 1.02591,
  1.02396, 1.02219, 1.02061, 1.01919, 1.01791, 1.01677, 1.01574, 1.01482, 1.01398, 1.01322,
  1.01253  };

  double ZF[91] = {

  0.27272, 0.27577, 0.27892, 0.28217, 0.28553, 0.28900, 0.29259, 0.29630, 0.30015, 0.30413,
  0.30826, 0.31255, 0.31700, 0.32163, 0.32644, 0.33145, 0.33666, 0.34210, 0.34777, 0.35369,
  0.35988, 0.36635, 0.37312, 0.38022, 0.38767, 0.39549, 0.40370, 0.41234, 0.42144, 0.43103,
  0.44115, 0.45183, 0.46311, 0.47505, 0.48769, 0.50108, 0.51528, 0.53035, 0.54634, 0.56332,
  0.58137, 0.60054, 0.62092, 0.64256, 0.66554, 0.68993, 0.71578, 0.74314, 0.77207, 0.80257,
  0.83468, 0.86836, 0.90359, 0.94031, 0.97842, 1.01780, 1.05830, 1.09974, 1.14193, 1.18465,
  1.22766, 1.27075, 1.31368, 1.35625, 1.39826, 1.43957, 1.48003, 1.51954, 1.55804, 1.59549,
  1.63188, 1.66722, 1.70153, 1.73485, 1.76725, 1.79876, 1.82946, 1.85939, 1.88862, 1.91719,
  1.94516, 1.97257, 1.99946, 2.02587, 2.05182, 2.07736, 2.10249, 2.12726, 2.15168, 2.17576,
  2.19952  };

  double ZFSD[91] = {

  0.20265, 0.20478, 0.20697, 0.20923, 0.21155, 0.21394, 0.21640, 0.21894, 0.22156, 0.22425,
  0.22704, 0.22991, 0.23288, 0.23595, 0.23912, 0.24240, 0.24579, 0.24930, 0.25293, 0.25669,
  0.26059, 0.26462, 0.26880, 0.27314, 0.27763, 0.28228, 0.28710, 0.29210, 0.29728, 0.30265,
  0.30821, 0.31396, 0.31991, 0.32605, 0.33240, 0.33893, 0.34565, 0.35255, 0.35962, 0.36683,
  0.37417, 0.38159, 0.38908, 0.39658, 0.40403, 0.41138, 0.41855, 0.42546, 0.43200, 0.43809,
  0.44360, 0.44842, 0.45243, 0.45551, 0.45756, 0.45846, 0.45814, 0.45655, 0.45364, 0.44944,
  0.44397, 0.43732, 0.42959, 0.42092, 0.41149, 0.40146, 0.39103, 0.38038, 0.36969, 0.35912,
  0.34880, 0.33886, 0.32936, 0.32037, 0.31193, 0.30405, 0.29671, 0.28990, 0.28360, 0.27776,
  0.27234, 0.26731, 0.26263, 0.25826, 0.25417, 0.25032, 0.24670, 0.24327, 0.24002, 0.23693,
  0.23398  };


  float h,x,delta;
  float c1,c2,e2,e4,e6;
  int n;
  // Bayesian statistics tells us to modify I/sigma by subtracting off sigma/2S
  // where S is the mean intensity in the resolution shell
  h = I/sigma - 0.5*sigma/S;
  // reject as unphysical reflections for which I < -3.7 sigma, or h < -4.0
  if (I/sigma < -3.7 || h < -4.0 ) {
	  nrej++;
	  if (debug) printf("unphys: %f %f %f %f\n",I,sigma,S,h);
	  return(0);
  }
  else {
	  if (h < 5.0) {
		  // use look up table if -4.0 < h < 5.0
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
		  // if h > 5.0 use asymptotic formulas in French and Wilson Appendix, with extra term
		  e2 = 1.0/(h*h);
		  e4 = e2*e2;
		  e6 = e2*e4;
		  c1 = (1.0 - 0.375*e2 - 87.0/128.0*e4 - 2889.0/1024.0*e6)*sqrt(h);
          c2 = sqrt((273.0/128.0*e6 + 15.0/32.0*e4 + 0.25*e2)*h);

		  J = h*sigma*(1.0 - 0.5*e2 - 0.75*e4 - 3.0*e6);
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

void Htest( HKL_data<data32::I_sigI> isig, Mat33<int> twinop, int &itwin, int scalefac, String s, CCP4Program& prog, bool debug )
{
	typedef clipper::HKL_data_base::HKL_reference_index HRI;
    double HT=0.0;
    double HT2=0.0;
    double NT=0.0;
	std::vector<int> pdf(50,0);
    std::vector<int> cdf(20,0);
    Vec3<int> jhkl, jhkl2;
    for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
	    if ( !isig[ih].missing() && !ih.hkl_class().centric() ) {
		    HKL hkl = ih.hkl();
		    jhkl[0] = hkl.h();
		    jhkl[1] = hkl.k();
		    jhkl[2] = hkl.l();
		    HKL twin;
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


// convert twinning operator from a matrix to a string
// code modified from Symop::format

void MatrixToString( Mat33<int> op, String &s )
{
  String t, hkl="hkl";
  for ( int i = 0; i < 3; i++ ) {
    t = "";
	for ( int j = 0; j < 3; j++ ) {
      if ( op(i,j) != 0 ) {
	    t += ( op(i,j) > 0 ) ? "+" : "-";
	    if ( abs( op(i,j) ) != 12 )
	      t += String::rational( fabs( float( op(i,j) )/12.0 ), 24 );
	    t += hkl[j];
      }
    }
    s += t.substr( ( t[0] == '+' ) ? 1 : 0 );
    if ( i < 2 ) s+= ", ";
  }
}

