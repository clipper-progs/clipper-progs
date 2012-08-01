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
#include "ccp4/ccp4_general.h"
#include "ccp4/ccp4_program.h"
#include "clipper/ccp4/ccp4_mtz_io.h"
#include "ccp4/cmtzlib.h"
#include "ccp4/csymlib.h"
#include <iostream>
#include <math.h>
#include "ccp4/ccp4_fortran.h"
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
#include "ctruncate_wilson.h"

#include <mmdb/mmdb_tables.h>


using namespace clipper;
using namespace ctruncate;

// replacement for Wilson/Truncate


int main(int argc, char **argv)
{
  CCP4Program prog( "ctruncate", "1.6.1", "$Date: 2012/08/01" );
  
  // defaults
  clipper::String outfile = "ctruncate_out.mtz";
  clipper::String outcol = "";
  clipper::String freecol = "/*/*/[FreeR_flag]";
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
  bool freein = false;
  bool amplitudes = false;
  bool anomalous = false;

  int mtzinarg = 0;
  int mtzoutarg = 0;
  int nbins = 60;
  int ncbins = 60;
  int nresidues = 0;
	int nprm = 60;

  clipper::Resolution reso_Patt = clipper::Resolution( 4.0 );
  clipper::Resolution reso_Twin = clipper::Resolution( 0.1 );
  clipper::Resolution reso_trunc = clipper::Resolution( 0.1 );

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
    } else if ( args[arg] == "-freein" ) {
      freein = true;
      if ( ++arg < args.size() ) freecol = args[arg];
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
    } else if ( args[arg] == "-twinres" ) {
        if ( ++arg < args.size() ) reso_Twin = clipper::Resolution( clipper::String(args[arg]).f() );
    } else if ( args[arg] == "-reso" ) {
        if ( ++arg < args.size() ) reso_trunc = clipper::Resolution( clipper::String(args[arg]).f() );
    } else if ( args[arg] == "-no-aniso" ) {
      aniso = false;
    } else if ( args[arg] == "-amplitudes" ) {
      amplitudes = true;
	} else if ( args[arg] == "-debug" ) {
      debug = true;
    } else if ( args[arg] == "-i" ) {
      CCP4::ccp4_prog_info();
      return(0);
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

  //mtzfile.open_read( args[1] );
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
  HKL_data<data32::Flag> free(hklinf);

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
  if (freein) mtzfile.import_hkl_data( free, freecol );

  prog.summary_beg();
  printf("\nCRYSTAL INFO:\n\n");
  std::cout << "Crystal/dataset names: " << mtzfile.assigned_paths()[0].notail() << "\n"; 
  /*std::cout << mtzfile.assigned_paths()[0].notail().notail().tail() << "\n";  // crystal name
  std::cout << mtzfile.assigned_paths()[0].notail().tail() << "\n";  //dataset name
  String xtlname = mtzfile.assigned_paths()[0].notail().notail().tail();
  String setname = mtzfile.assigned_paths()[0].notail().tail();*/
  prog.summary_end();
  printf("\n");

  // need this mumbo-jumbo in order to write to output file

  if ( outcol[0] != '/' ) outcol = mtzfile.assigned_paths()[0].notail()+"/"+outcol;

  // hkl_list contains only those (h,k,l) for which at least data column is not NaN.
  // hkl_info contains all (h,k,l) out to the resolution limit regardless of whether there is any measured data.
  // will need hkl_list when writing to output file.
  clipper::Spacegroup spgr = mtzfile.spacegroup();
  clipper::Cell      cell1 = mtzfile.cell();
  clipper::Resolution reso = mtzfile.resolution();

	// limit resolution for truncation, and hence output
  reso_trunc = clipper::Resolution( clipper::Util::max( reso.limit(), reso_trunc.limit() ) );

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
  mtz1 = CMtz::MtzGet(args[mtzinarg].c_str(), read_refs);

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
		prog.summary_beg(); printf("CELL CONTENTS:\n\n");
		
		clipper::SEQfile seqf;
		seqf.read_file( ipseq );
		
		ctruncate::Matthews cmath(true,false);
		int nmol = cmath(cell1, spgr, seqf, resopt);
		std::cout << "Expected number of molecules in ASU : " << nmol << std::endl;
		prog.summary_end();
		cmath.summary();
	} else if (nresidues > 0) {		
		prog.summary_beg(); printf("CELL CONTENTS:\n\n");
		ctruncate::Matthews cmath(true,false);
		int nmol = cmath(cell1, spgr, nresidues, resopt);
		std::cout << "Expected number of molecules in ASU : " << nmol << std::endl;
		prog.summary_end();
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

  int npeak = 5;

  std::vector<int> ppks = pksch( patterson );

  float top_peak = patterson.get_data( ppks[0] );
  float next_peak = patterson.get_data( ppks[1] );
  clipper::Coord_frac c0 = patterson.coord_of( ppks[1] ).coord_frac(grid);
  float ratio = next_peak/top_peak;
  float dist2 = pow(c0[0], 2.0) + pow(c0[1], 2.0) + pow(c0[2], 2.0);
  // look for peaks > 20% of origin peak and at least 0.1 distant from origin
  // precentage estimate is Zwartz CCP4 Newsletter 42
  const double aval = 0.0679;
  const double bval = 3.56;
  double pval = (1.0 - std::exp(-std::pow(ratio/(aval*(1.0-ratio)),-bval)) )*100.0;;
  prog.summary_beg();
  printf("\n\nTRANSLATIONAL NCS:\n");
  if ( debug || (ratio > 0.2 && dist2 > 0.01) ) { 
	  printf("Translational NCS has been detected at (%6.3f, %6.3f, %6.3f).\n  The probability, based on peak ratio, is %5.2f that this is\n by chance (with resolution limited to %5.2f A). \n", c0[0],c0[1],c0[2],pval,reso_Patt.limit() );
      printf("This will have a major impact on the twinning estimates and effectiveness of the truncate procedure\n");
      printf("Peak Ratio = %5.2f \n",ratio);
      printf("Peak Vector = (%6.3f, %6.3f, %6.3f)\n",c0[0],c0[1],c0[2]);
  }
  else {
	  printf("No translational NCS detected (with resolution limited to %5.2f A)\n", reso_Patt.limit() );
      if ( dist2 > 0.01 ) printf("Top off origin peak at (%6.3f, %6.3f, %6.3f) with a probability of %5.2f\n",c0[0],c0[1],c0[2],100.0-pval); 
  }
  prog.summary_end();
  printf("\n");

  // falloff and completeness for input intensities
    //if (!amplitudes) yorgo_modis_plot(isig,maxres,60,prog);

    // anisotropy estimation
    clipper::U_aniso_orth uao;
	  double Itotal = 0.0;
    
	  prog.summary_beg();
    printf("\nANISOTROPY ESTIMATION (using intensities):\n");

	  for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {  
	  	  double I = isig[ih].I();
	      double sigI = isig[ih].sigI();
	      if ( I > 0.0 ) Itotal += I;
     	  ianiso[ih] = clipper::data32::I_sigI( I, sigI );
	  }
    
    AnisoCorr<Iscale_logLikeAniso<float>, clipper::datatypes::I_sigI<float>, float > llscl(ianiso);
	  
    uao = llscl.u_aniso_orth(Scaling::I);
         
    // Eigenvalue calculation
    AnisoDirection<float> direction(uao);

    std::vector<float> v = direction.eigenValues();
    float max = direction.max();

    printf("\nEigenvalues: %8.4f %8.4f %8.4f\n", v[0],v[1],v[2]);
    printf("Eigenvalue ratios: %8.4f %8.4f %8.4f\n", v[0]/max, v[1]/max, v[2]/max);
    if ( v[0] <= 0.0 ) CCP4::ccperror(1, "Anisotropy correction failed - negative eigenvalue.");
    invopt = maxres*v[0]/v[2];
    resopt = 1.0/sqrt(invopt);
    printf("Resolution limit in weakest direction = %7.3f A\n",resopt);
    if ( v[0]/max < 0.5 ) printf("\nWARNING! WARNING! WARNING! Your data is severely anisotropic\n");
    prog.summary_end();
    printf("\n");
 
    //std::cout << "estimated Biso: " << clipper::Util::u2b(0.5*uao.u_iso()) << std::endl;
          
    //if (!amplitudes) yorgo_modis_plot(isig,maxres,60,prog, uao);
	  
    printf("\nAnisotropic U (orthogonal coords):\n\n");
    printf("| %8.4f %8.4f %8.4f |\n", uao(0,0) ,  uao(0,1) ,  uao(0,2)  );
    printf("| %8.4f %8.4f %8.4f |\n", uao(1,0) ,  uao(1,1) ,  uao(1,2)  );
    printf("| %8.4f %8.4f %8.4f |\n", uao(2,0) ,  uao(2,1) ,  uao(2,2)  );
	  
    clipper::U_aniso_frac uaf = uao.u_aniso_frac( cell1 );

    printf("\nAnisotropic U scaling (fractional coords):\n\n"); 

    printf("| %11.3e %11.3e %11.3e |\n", uaf(0,0) ,  uaf(0,1) ,  uaf(0,2)  );
    printf("| %11.3e %11.3e %11.3e |\n", uaf(1,0) ,  uaf(1,1) ,  uaf(1,2)  );
    printf("| %11.3e %11.3e %11.3e |\n", uaf(2,0) ,  uaf(2,1) ,  uaf(2,2)  );

    printf("\nAnisotropic B scaling (fractional coords):\n\n"); 

    printf("| %11.3e %11.3e %11.3e |\n",clipper::Util::u2b( uaf(0,0) ), clipper::Util::u2b( uaf(0,1) ), clipper::Util::u2b( uaf(0,2) ) );
    printf("| %11.3e %11.3e %11.3e |\n",clipper::Util::u2b( uaf(1,0) ), clipper::Util::u2b( uaf(1,1) ), clipper::Util::u2b( uaf(1,2) ) );
    printf("| %11.3e %11.3e %11.3e |\n",clipper::Util::u2b( uaf(2,0) ), clipper::Util::u2b( uaf(2,1) ), clipper::Util::u2b( uaf(2,2) ) );

    if (aniso)  {
        double FFtotal = 0.0;
        printf("\nANISOTROPY CORRECTION (using intensities):\n");

        Iscale_aniso<float> sfscl( 3.0 );
        sfscl( ianiso);
        
        printf("\nAnisotropic scaling (orthogonal coords):\n\n");
        
        printf("|%8.4f %8.4f %8.4f |\n", sfscl.u_aniso_orth(Scaling::I)(0,0), sfscl.u_aniso_orth(Scaling::I)(0,1), sfscl.u_aniso_orth(Scaling::I)(0,2) );
        printf("|%8.4f %8.4f %8.4f |\n", sfscl.u_aniso_orth(Scaling::I)(1,0), sfscl.u_aniso_orth(Scaling::I)(1,1), sfscl.u_aniso_orth(Scaling::I)(1,2) );
        printf("|%8.4f %8.4f %8.4f |\n", sfscl.u_aniso_orth(Scaling::I)(2,0), sfscl.u_aniso_orth(Scaling::I)(2,1), sfscl.u_aniso_orth(Scaling::I)(2,2) );

	printf("\n");

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

  // truncate anisotropically corrected data at resolution limit in weakest direction
  for ( HRI ih = ianiso.first(); !ih.last(); ih.next() ) {
	  if (ih.invresolsq() > invopt) ianiso[ih].set_null();  
  }

  // calculate moments of Z using truncate methods

	moments_Z(isig, maxres, nbins);

	
	HKL_data<data32::I_sigI>& iptr = (aniso) ? ianiso : isig;
	
  printf("\nTWINNING ANALYSIS:\n\n");
    float hval(0.0), lval(0.0);
    //reduced resolution range for twinning tests
    reso_Twin = clipper::Resolution( clipper::Util::max( clipper::ftype(resopt), reso_Twin.limit() ) );
  HKL_info hklt(spgr, cell1, reso_Twin, true);
    HKL_data<data32::I_sigI> itwin(hklt);
    for ( HRI ih = itwin.first() ; !ih.last() ; ih.next() ) {
        itwin[ih] = iptr[ih.hkl()];
    }

	{
		printf("\nData has been truncated at %6.2f A resolution\n",reso_Twin.limit());
		if (aniso) printf("Anisotropy correction has been applied before calculating twinning tests\n\n");
		else printf("Anisotropy correction has not been applied before calculating twinning tests\n\n");
	}
	
  // H test for twinning

  if (twintest != "table") {
	  hval = Htest_driver_fp(itwin, debug);
  }

  if (twintest != "first_principles") {
	  hval = Htest_driver_table(itwin, debug);
  }

	 // L test for twinning
	lval = Ltest_driver(iptr, debug);
 
    prog.summary_beg();
    twin_summary(hval,lval);
        prog.summary_end(); 
    printf("\n");
  //printf("Starting parity group analysis:\n");

  //Parity group analysis

	ctruncate::parity(isig, maxres, nbins);	
	
	// Ice rings
	ctruncate::Rings icerings;
	icerings.DefaultIceRings();
	icerings.ClearSums();
    float iceTolerance = 4.0f;
	
	for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
		if ( !isig[ih].missing() )
		if (! ih.hkl_class().centric() ) {
			double reso = ih.hkl().invresolsq(hklinf.cell());
			int ice = icerings.InRing(reso);
			if ( ice != -1 ) icerings.AddObs(ice,isig[ih],reso); //symmetry?
		}
	}
	
	for (int i = 0; i != icerings.Nrings(); ++i) icerings.SetReject(i, true);
	
  //Wilson pre
	//std::vector<float> wilson(2,0.0f);
	std::vector<double> param_gauss( 2, 0.0 );
	//int nprm2 = std::floor(nprm/3.0);
	int nprm2 = 12;
	
	HKL_data<data32::I_sigI> xsig(hklinf);  // knock out ice rings and centric
	
	{
		int nbins = 60;
		std::vector<float> sumov(nbins,0.0);
		std::vector<float> summeas(nbins,0.0);
		std::vector<float> iceov(icerings.Nrings(),0.0);
		std::vector<float> icemeas(icerings.Nrings(),0.0);
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			double reso = ih.hkl().invresolsq(hklinf.cell());
			int ring=icerings.InRing(reso);
			if ( ring == -1 ) {
				int bin = int( double(nbins) * ih.invresolsq() / maxres - 0.001);
				//if (bin >= nbins || bin < 0) printf("Warning: (completeness) illegal bin number %d\n", bin);
				if ( bin < nbins && bin >= 0 ) {
					sumov[bin] += 1.0;
					if ( !isig[ih].missing() ) summeas[bin] += 1.0;
				}
			} else {
				if ( ring <  icerings.Nrings() && ring >= 0 ) {
					iceov[ring] += 1.0;
					if ( !isig[ih].missing() ) icemeas[ring] += 1.0;
				}
			}
		}
		// smooth the mean values, also means values are over 3 bins so hopefully get fewer empty bins
		for (int i = 0 ; i != nbins ; ++i ) {
			if ( i == 0 ) {
				sumov[i] += sumov[i+1];
				summeas[i] += summeas[i+1];
			} else if ( i == (nbins-1) ) {
				sumov[i] += sumov[i-1];
				summeas[i] += summeas[i-1];
			} else {
				sumov[i] += sumov[i+1] + sumov[i-1];
				summeas[i] += summeas[i+1] + summeas[i-1];
			}
		}
		
		HKL_data<data32::I_sigI> tr1(hklinf);
		tr1 = data32::I_sigI(1.0f,1.0f);
		for ( HRI ih = tr1.first(); !ih.last(); ih.next() ) tr1[ih].scale( ih.hkl_class().epsilon() ); //tage into account  epsilon
		
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			double reso = ih.hkl().invresolsq(hklinf.cell());
			xsig[ih] = clipper::data32::I_sigI( (isig[ih].I()), isig[ih].sigI() );
			if ( ih.hkl_class().centric() ) xsig[ih].I() = xsig[ih].sigI() = clipper::Util::nan(); // loose centrics
			if ( icerings.InRing(reso) != -1 ) 
				if ( icerings.Reject( icerings.InRing(reso) ) ) xsig[ih].I() = xsig[ih].sigI() = clipper::Util::nan(); // loose ice rings
		}		
		
		// calc scale
		//std::vector<double> param_gauss( 2, 0.0 );
		clipper::BasisFn_log_gaussian gauss;
		clipper::TargetFn_scaleLogI1I2<clipper::data32::I_sigI,clipper::data32::I_sigI> tfn_gauss( xsig, tr1);
		ResolutionFn rfn_gauss( hklinf, gauss, tfn_gauss, param_gauss );
		param_gauss = rfn_gauss.params();
				
		//datatypes::Compute_scale_u_iso<float> compute_s(1.0,param_gauss[1]);
		//xsig.compute( isig, compute_s ); scale xsig for gaussian decay
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			if ( !isig[ih].missing() ) {
				double reso = ih.hkl().invresolsq(hklinf.cell());
				xsig[ih] = clipper::data32::I_sigI( (isig[ih].I()*exp(-param_gauss[1]*reso) ), isig[ih].sigI()*exp(-param_gauss[1]*reso) );
			}
		}
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			double reso = ih.hkl().invresolsq(hklinf.cell());
			if ( ih.hkl_class().centric() ) xsig[ih].I() = xsig[ih].sigI() = clipper::Util::nan(); // loose centrics
			if ( icerings.InRing(reso) != -1 ) 
				if ( icerings.Reject( icerings.InRing(reso) ) ) xsig[ih].I() = xsig[ih].sigI() = clipper::Util::nan(); // loose ice rings
		}		
		
		std::vector<double> params_ice( nprm2, 1.0 );
		clipper::BasisFn_spline basis_ice( xsig, nprm2, 2.0 );
		TargetFn_meanInth<clipper::data32::I_sigI> target_ice( xsig, 1 );
		clipper::ResolutionFn Sigma( hklinf, basis_ice, target_ice, params_ice );
		
		for (int i = 0; i != icerings.Nrings(); ++i) {
			bool reject = false;
			float reso = icerings.MeanSSqr(i);
			if ( reso <= maxres && icerings.MeanSigI(i) > 0.0f ) {
				float imean = icerings.MeanI(i);
				float sigImean = icerings.MeanSigI(i);
				float expectedI = exp( log(basis_ice.f_s( reso, Sigma.params()) ) + param_gauss[1]*reso);
				if ((imean-expectedI)/sigImean > iceTolerance) reject = true;
			}
			icerings.SetReject(i, reject);
		}
		
		bool icer = false;
		for ( int i = 0; i != icerings.Nrings(); ++i) if (icerings.Reject(i) ) icer = true;
		printf("\n");
		prog.summary_beg();
		printf("\nICE RINGS:\n\n");
		if (icer) printf("Possible Ice Rings\n\n");
		else printf("No Ice Rings detected\n\n");
		prog.summary_end();
                printf("\n"); 
		
		if ( icerings.MeanSSqr(0) <= maxres && icerings.MeanSigI(0) > 0.0f ) {
			printf("Ice Ring Summary:\n");
			printf(" reso mean_I mean_Sigma Estimated_I Zscore Completeness Ave_Completeness\n");
			for ( int i = 0; i != icerings.Nrings(); ++i) {
				float reso = icerings.MeanSSqr(i);
				if ( reso <= maxres && icerings.MeanSigI(i) > 0.0f ) {
					int bin = int( double(nbins) * reso / maxres - 0.001);
					float imean = icerings.MeanI(i);
					float sigImean = icerings.MeanSigI(i);
					float expectedI = exp( log(basis_ice.f_s( reso, Sigma.params()) ) + param_gauss[1]*reso);
					printf("%6.2f %-10.2f %-10.2f %-10.2f %-6.2f %-6.2f %-6.2f\n",1.0f/std::sqrt(reso),imean,sigImean,expectedI,(imean-expectedI)/sigImean,icemeas[i]/iceov[i],summeas[bin]/sumov[bin]);
				}
			}
		}
	}
	
	// Sigma or Normalisation curve
	// calculate Sigma (mean intensity in resolution shell) 
	// use intensities uncorrected for anisotropy
	for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
		double reso = ih.hkl().invresolsq(hklinf.cell());
		xsig[ih] = clipper::data32::I_sigI( (isig[ih].I()*exp(-param_gauss[1]*reso) ), isig[ih].sigI()*exp(-param_gauss[1]*reso) );
		if ( ih.hkl_class().centric() ) xsig[ih].I() = xsig[ih].sigI() = clipper::Util::nan(); // loose centrics
		if ( icerings.InRing(reso) != -1 ) 
			if ( icerings.Reject( icerings.InRing(reso) ) ) xsig[ih].I() = xsig[ih].sigI() = clipper::Util::nan(); // loose ice rings
	}

		
	std::vector<double> params( nprm2, 1.0 );
	clipper::BasisFn_spline basis_fo( xsig, nprm2, 2.0 );
	TargetFn_meanInth<clipper::data32::I_sigI> target_fo( xsig, 1 );
	clipper::ResolutionFn Sigma( hklinf, basis_fo, target_fo, params );
	
	for ( HRI ih = xsig.first(); !ih.last(); ih.next() ) {
		double reso = ih.hkl().invresolsq(hklinf.cell());
		xsig[ih].I() = exp( log(Sigma.f(ih)) +param_gauss[1]*reso );  //not using multiplicity  THIS DOES NOT CONFORM!!!!!!
	}
	
	//Wilson plot
	//std::vector<float> wilson(2,0.0f);
	//nprm = 60*std::max(int(sqrt(float(Nreflections))),nbins );
	clipper::MMoleculeSequence seq;
	
	WilsonB wilson( WilsonB::BEST);
	if ( ipseq != "NONE" ) {
		clipper::SEQfile seqf;
		seqf.read_file( ipseq );
		seqf.import_molecule_sequence( seq );
		MPolymerSequence poly = seq[0];
		//wilson = wilson_plot(isig,poly,maxres,nprm, prog, xsig);
		wilson(isig,poly,&icerings);
	} else if (nresidues > 0) {
		//wilson = wilson_plot(isig,nresidues,maxres,nprm, prog, xsig);
		wilson(isig,nresidues,&icerings);
	} else {
		//wilson = wilson_plot(isig,maxres,nprm, prog, xsig);
		wilson(isig,&icerings);
	}		
	// end of Norm calc
	clipper::String comment("Smooth");
        printf("\n");
	prog.summary_beg();
	wilson.summary();
	prog.summary_end();
        printf("\n");
	wilson.plot(xsig,comment);
    
  // apply the Truncate procedure, unless amplitudes have been input

  // if something went wrong with Wilson scaling, B could be negative, giving exponentially large scaled SF's
  // so only scale if B positive
  float scalef = 1.0;

  if ( wilson.intercept() > 0 ) scalef = sqrt(wilson.intercept() );
  int nrej = 0; 

  if (!amplitudes) {
      if (anomalous) {
	      truncate( isig_ano, jsig_ano, fsig_ano, xsig, scalef, spg1, reso_trunc, nrej, debug );
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
          truncate( isig, jsig, fsig, xsig, scalef, spg1, reso_trunc, nrej, debug );
      }
  }
  printf("\n");
  prog.summary_beg();
  printf("\nINTENSITY TO AMPLITUDE CONVERSION:\n\n");
  printf("%d intensities have been rejected as unphysical\n", nrej);
  prog.summary_end();
  printf("\n");


  
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
  int ntw = cumulative_plot(isig, xsig);
	
  if (ntw > 2) {
	  prog.summary_beg();
	  printf("\nWARNING: ****  Cumulative Distribution shows Possible Twinning ****\n\n");
	  prog.summary_end();
  }


  // falloff calculation (Yorgo Modis)
	yorgo_modis_plot(fsig,maxres,60,prog,uao);
	
    // anomalous signal
    if (anomalous ) {
        if (amplitudes) AnomStats<float> anomstats(fsig_ano);
        else AnomStats<float> anomstats(isig_ano);
    }

	
	{
		ctruncate::PattPeak patt_peak(std::sqrt(maxres));
		
		float opt_res = patt_peak(xsig);
		
		float width_patt = 2.0f*patt_peak.sigma();
		
		float b_patt = 4.0f*clipper::Util::twopi2()*std::pow(width_patt/2.0f,2.0f);
		
		prog.summary_beg();
		std::cout << "Estimated Optical Resolution: " << opt_res << std::endl;
		prog.summary_end();
		
		
		
	}
	
	// I/Sigma and F/sqrt(Sigma) plots
	  if (!amplitudes) {
		//could use a clipper::Histogram
		clipper::Range<double> range(-5.0,10.0);
		/*for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			float eps = ih.hkl_class().epsilon();
			if (!isig[ih].missing() ) range.include(isig[ih].I()/(xsig[ih].I()));
						//if (!fsig[ih].missing() ) range.include(std::sqrt(eps)*fsig[ih].f()/std::sqrt(xsig[ih].I()) );
		}*/
		clipper::Histogram Icount(range,200);
		clipper::Histogram Jcount(range,200);
		clipper::Histogram Fcount(range,200);
		clipper::Histogram IScount(range,200);
		clipper::Histogram JScount(range,200);
		clipper::Histogram FScount(range,200);
		
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			float eps = ih.hkl_class().epsilon();
			if (!isig[ih].missing() ) {
				Icount.accumulate(isig[ih].I()/xsig[ih].I());
			}
		}
			
		for ( HRI ih = jsig.first(); !ih.last(); ih.next() ) {
			float eps = ih.hkl_class().epsilon();
			if (!jsig[ih].missing() ) {
				Jcount.accumulate(jsig[ih].I()/xsig[ih].I());
			}
		}	
		
		for ( HRI ih = fsig.first(); !ih.last(); ih.next() ) {
			float eps = ih.hkl_class().epsilon();
			if (!fsig[ih].missing() ) {
				Fcount.accumulate(fsig[ih].f()/(scalef*std::sqrt(xsig[ih].I()) ) );
			}
		}	
		
		  
		  for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			  float eps = ih.hkl_class().epsilon();
			  if (!isig[ih].missing() ) {
				  IScount.accumulate(isig[ih].I()/isig[ih].sigI() );
			  }
		  }
		  
		  for ( HRI ih = jsig.first(); !ih.last(); ih.next() ) {
			  float eps = ih.hkl_class().epsilon();
			  if (!jsig[ih].missing() ) {
				  JScount.accumulate(jsig[ih].I()/jsig[ih].sigI() );
			  }
		  }	
		  
		  for ( HRI ih = fsig.first(); !ih.last(); ih.next() ) {
			  float eps = ih.hkl_class().epsilon();
			  if (!fsig[ih].missing() ) {
				  FScount.accumulate(fsig[ih].f()/fsig[ih].sigf() );
			  }
		  }	
		  
		printf("$TABLE: Phil plot:\n");
		printf("$GRAPHS");
		printf(": Phil plot - normalised values:A:1,2,3,4:\n");  
		printf(": Phil plot - vs sigma:A:1,5,6,7:\n$$");  

		printf(" Value Io/Sigma I/Sigma F/Sigma**0.5 Io/sigIo I/sigI F/sigF$$\n$$\n");

		for ( int i=0; i!=200; ++i ) {
			float res = range.min()+float(i)*(range.max()-range.min())/200.0; 
			printf("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", 
					res,
					Icount.y(res),
				    Jcount.y(res),
					Fcount.y(res),
				   IScount.y(res),
				   JScount.y(res),
				   FScount.y(res)
				   );
		}
		
		printf("$$\n\n");
	}		

  // output data
  if (!amplitudes) {
      //mtzout.open_append( args[mtzinarg], outfile );
	  mtzout.open_write( outfile );
	  mtzout.export_crystal ( cxtl, outcol );
      mtzout.export_dataset ( cset, outcol );
      mtzout.export_hkl_info( hkl_list );
      //mtzout.export_hkl_data( jsig, outcol );
	  clipper::String labels;
	  if (appendcol == "") labels = outcol + "[F,SIGF]";
	  else labels = outcol + "[F_" + appendcol + ",SIGF_" + appendcol + "]";
	  mtzout.export_hkl_data( fsig, labels );

      if (freein) {
	  if (appendcol != "") {
		  String::size_type loc = freecol.find("]",0);
		  freecol.insert(loc,"_"+appendcol);
	  }
	  mtzout.export_hkl_data( free, outcol + freecol.tail() );
      }

      if (anomalous) {
		  if (appendcol == "") labels = outcol + "[DANO,SIGDANO]";
	      else labels = outcol + "[DANO_" + appendcol + ",SIGDANO_" + appendcol + "]";
		  mtzout.export_hkl_data( Dano, labels );
	      if (appendcol == "") labels = outcol + "[F(+),SIGF(+),F(-),SIGF(-)]";
	      else labels = outcol + "[F_" + appendcol + "(+),SIGF_" + appendcol + "(+),F_" + appendcol + "(-),SIGF_" + appendcol + "(-)]";
	      mtzout.export_hkl_data( fsig_ano, labels );
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
      mtz2 = CMtz::MtzGet(args[mtzoutarg].c_str(), read_refs);
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





