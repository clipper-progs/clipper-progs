//
//     CTRUNCATE
//     Copyright (C) 2006-2011 Norman Stein, Charles Ballard
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
#include "ctruncate_truncate.h"
#include "ctruncate_utils.h"
#include "ctruncate_twin.h"
#include "ctruncate_parity.h"
#include "ctruncate_moments.h"
#include "ctruncate_analyse.h"
#include "ctruncate_matthews.h"

#include <mmdb/mmdb_tables.h>


using namespace clipper;
using namespace ctruncate;

// replacement for Wilson/Truncate


extern "C" void FORTRAN_CALL ( YYY_CELL2TG, yyy_cell2tg,
	   ( clipper::Cell& cell, double& sc_tol, int& ng, int *uu_g, int *u_g, 
		 int& lc, int& nc, int& nc2, int *uu_c, double *sc_c, int& ivb, int& ierr ),
	   ( clipper::Cell& cell, double& sc_tol, int& ng, int *uu_g, int *u_g, 
		 int& lc, int& nc, int& nc2, int *uu_c, double *sc_c, int& ivb, int& ierr ),
	   ( clipper::Cell& cell, double& sc_tol, int& ng, int *uu_g, int *u_g, 
		 int& lc, int& nc, int& nc2, int *uu_c, double *sc_c, int& ivb, int& ierr ));


int main(int argc, char **argv)
{
  CCP4Program prog( "ctruncate", "1.0.13", "$Date: 2011/02/07" );
  
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
	int nprm = 60;

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
	clipper::CCP4MTZ_type_registry::add_type( "ISym", "Y", 1.0);
	clipper::CCP4MTZ_type_registry::add_group( "ISym", "ISYM" );
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
  HKL_data<data32::ISym> freidal_sym(hklinf);

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

	if ( ipseq != "NONE" ) {
		clipper::SEQfile seqf;
		seqf.read_file( ipseq );
		
		ctruncate::Matthews cmath(true,false);
		int nmol = cmath(cell1, spgr, seqf, resopt);
		std::cout << "Expected number of molecules in ASU : " << nmol << std::endl;
		cmath.summary();
	} else if (nresidues > 0) {		
		ctruncate::Matthews cmath(true,false);
		int nmol = cmath(cell1, spgr, nresidues, resopt);
		std::cout << "Expected number of molecules in ASU : " << nmol << std::endl;
		cmath.summary();
		
	}
	
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


  // calculate moments of Z using truncate methods

	moments_Z(isig, maxres, nbins);

  prog.summary_beg();
  printf("\n\nTWINNING ANALYSIS:\n\n");
  bool itwin = false;

	{
		//printf("\nData has been truncated at %6.2f A resolution\n",resopt);
		if (aniso) printf("Anisotropy correction has been applied before calculating H-test\n\n");
		else printf("Anisotropy correction has not been applied before calculating H-test\n\n");
	}
	
  // H test for twinning

  if (twintest != "table") {
	  itwin = Htest_driver_fp(ianiso, prog, debug);
  }

  if (twintest != "first_principles") {
	  itwin = Htest_driver_fp(ianiso, prog, debug);
  }

	 // L test for twinning
	{
		//printf("\nData has been truncated at %6.2f A resolution\n",resopt);
		if (aniso) printf("Anisotropy correction has been applied before calculating L-test\n\n");
		else printf("Anisotropy correction has not been applied before calculating L-test\n\n");
	}
 
	itwin = Ltest_driver(ianiso, prog, debug);
   

  //printf("Starting parity group analysis:\n");

  //Parity group analysis

	ctruncate::parity(isig, maxres, nbins);

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
	  numatoms[4] = int(0.05*nresidues);
  }


  // Wilson plot
	
  	nprm = std::max(int(sqrt(float(Nreflections))),nbins );
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
  b = 0.0;

  if ( wi.size() > 200 && maxres > maxres_scaling) {               // 3.5 Angstroms
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
  
	// Sigma or Normalisation curve
	// calculate Sigma (mean intensity in resolution shell) 
	// use intensities uncorrected for anisotropy
	
	int nprm2 = std::floor(nprm/3.0);
	
	HKL_data<data32::I_sigI> xsig(hklinf);  // knock out ice rings and centric
	for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
		double reso = ih.hkl().invresolsq(hklinf.cell());
		xsig[ih] = clipper::data32::I_sigI( isig[ih].I(), isig[ih].sigI() );
		if ( ih.hkl_class().centric() ) xsig[ih].I() = xsig[ih].sigI() = clipper::Util::nan(); // loose centrics
		if ( icerings.InRing(reso) != -1 )
			xsig[ih].I() = xsig[ih].sigI() = clipper::Util::nan(); // loose ice rings
	}		
	std::vector<double> params_ice( nprm2, 1.0 );
	clipper::BasisFn_spline basis_fo( xsig, nprm2, 2.0 );
	TargetFn_meanInth<clipper::data32::I_sigI> target_fo( xsig, 1 );
	clipper::ResolutionFn Sigma( hklinf, basis_fo, target_fo, params_ice );
	
	// end of Norm calc

	
	// wilson plot plus Norm curve
  printf("$TABLE: Wilson plot:\n");
  printf("$GRAPHS");
  //printf(": Wilson plot:0|0.1111x-7|-5:1,2:\n$$");  // limits hardwired
  printf(": Wilson plot - estimated B factor = %5.1f :A:1,2,3,4:\n$$", a);  
  printf(" 1/resol^2 ln(I/I_th) Sigma Overall-B $$\n$$\n");

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
	  printf("%10.5f %10.5f %10.5f %10.5f \n", res,log(basis_fo_wilson.f_s( res, wilsonplot.params() ))-log(totalscat),
			 log(basis_fo.f_s( res, Sigma.params() ))-log(totalscat),-0.5*a*res-b);
  }

  printf("$$\n\n");

	if (debug) {
		FILE *ftestfile;
		ftestfile = fopen("sigma.txt","w");
		for (int i=0; i!=nprm; ++i) {
			double res = maxres * pow( double(i+1)/nprm, 0.666666 );
			fprintf(ftestfile,"%10.6f %10.6f %10.6f \n", res, basis_fo_wilson.f_s( res, wilsonplot.params() ),basis_fo.f_s( res, Sigma.params() ));
		}
		fclose(ftestfile);
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
			  freidal_sym[ih].isym() = 0; //mimic old truncate
			  if ( !Util::is_nan(fsig_ano[ih].f_pl() )  &&  !Util::is_nan(fsig_ano[ih].f_mi() ) ) {
			      fsig[ih].f() = 0.5 * ( fsig_ano[ih].f_pl() + fsig_ano[ih].f_mi() );
			      fsig[ih].sigf() = 0.5 * sqrt( pow( fsig_ano[ih].sigf_pl(), 2 ) + pow( fsig_ano[ih].sigf_mi(), 2 ) );
			      Dano[ih].d() = fsig_ano[ih].f_pl() - fsig_ano[ih].f_mi();
			      Dano[ih].sigd() = 2.0 * fsig[ih].sigf();
				  freidal_sym[ih].isym() = 0;
		      }
		      else if ( !Util::is_nan(fsig_ano[ih].f_pl() ) ) {
			      fsig[ih].f() = fsig_ano[ih].f_pl();
			      fsig[ih].sigf() = fsig_ano[ih].sigf_pl();
				  freidal_sym[ih].isym() = 1;
		      }
		      else if ( !Util::is_nan(fsig_ano[ih].f_mi() ) ) {
			      fsig[ih].f() = fsig_ano[ih].f_mi();
			      fsig[ih].sigf() = fsig_ano[ih].sigf_mi();
				  freidal_sym[ih].isym() = 2;
		      }
		      else if ( !isig[ih].missing() && iwarn != 1 ) {
			      printf("\nWARNING: Imean exists but I(+), I(-) do not\n\n");
			      iwarn = 1;
		      }
		      if ( ih.hkl_class().centric() ) {
			      Dano[ih].d() = 0.0;
			      Dano[ih].sigd() = 0.0;
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
  // moments_Z(ianiso,invopt,nbins,prog);


  // construct cumulative distribution function for intensity (using Z rather than E)
  int ntw = cumulative_plot(isig, Sigma);
	
  if (ntw > 2) {
	  prog.summary_beg();
	  printf("\nWARNING: ****  Cumulative Distribution shows Possible Twinning ****\n\n");
	  prog.summary_end();
  }


  // falloff calculation (Yorgo Modis)

  Mat33<float> transf;
  //Cell cell = hklinf.cell();

  // calculate the matrix that transforms the cell coordinates h,k,l to Cartesian coordinates.
  tricart (cell1, transf);

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
  float enumdir[3][60];

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

	
	{
		ctruncate::PattPeak patt_peak(std::sqrt(maxres));
		
		float opt_res = patt_peak(basis_fo, Sigma);
		
		float width_patt = 2.0f*patt_peak.sigma();
		
		float b_patt = 2.0f*clipper::Util::twopi()*clipper::Util::twopi()*std::pow(width_patt/2.0f,2.0f);
		
		prog.summary_beg();
		std::cout << "Estimated Optical Resolution: " << opt_res << std::endl;
		prog.summary_end();
		
		
		
	}
	

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
	      if (appendcol == "") labels = outcol + "[ISYM]";
	      else labels = outcol + "[ISYM_" + appendcol + "]";
		  mtzout.export_hkl_data( freidal_sym, labels );
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
	  CMtz::MtzPut( mtz2, outfile.c_str() );
      CMtz::MtzFree( mtz2 );
  }
  CMtz::MtzFree( mtz1 );
  prog.set_termination_message( "Normal termination" );

  return(0);
}





