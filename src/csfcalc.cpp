// Clipper app to perform structure factor calculation
/* Copyright 2003-2004 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-mmdb.h>

extern "C" {
  #include <stdio.h>
  #include <stdlib.h>
}


int main( int argc, char** argv )
{
  CCP4Program prog( "csfcalc", "0.1", "$Date: 2004/06/01" );

  // defaults
  enum ANISO { NONE, FOBS, FCAL };
  clipper::String title;
  clipper::String ippdb = "NONE";
  clipper::String ipfile = "NONE";
  clipper::String ipcolfo = "NONE";
  clipper::String ipcolfree = "NONE";
  clipper::String opfile = "sfcalc.mtz";
  clipper::String opcol = "sfcalc";
  bool bulk  = true;
  ANISO aniso = NONE;
  int freeflag = 0;
  int n_refln = 1000;
  int n_param = 20;
  int verbose = 0;

  // command input
  CCP4CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if ( args[arg] == "-title" ) {
      if ( ++arg < args.size() ) title = args[arg];
    } else if ( args[arg] == "-pdbin" ) {
      if ( ++arg < args.size() ) ippdb = args[arg];
    } else if ( args[arg] == "-mtzin" ) {
      if ( ++arg < args.size() ) ipfile = args[arg];
    } else if ( args[arg] == "-colin-fo" ) {
      if ( ++arg < args.size() ) ipcolfo = args[arg];
    } else if ( args[arg] == "-colin-free" ) {
      if ( ++arg < args.size() ) ipcolfree = args[arg];
    } else if ( args[arg] == "-mtzout" ) {
      if ( ++arg < args.size() ) opfile = args[arg];
    } else if ( args[arg] == "-colout" ) {
      if ( ++arg < args.size() ) opcol = args[arg];
    } else if ( args[arg] == "-free-flag" ) {
      if ( ++arg < args.size() ) freeflag = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-num-reflns" ) {
      if ( ++arg < args.size() ) n_refln = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-num-params" ) {
      if ( ++arg < args.size() ) n_param = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-no-bulk" ) {
      bulk = false;
    } else if ( args[arg] == "-aniso-obs" ) {
      aniso = FOBS;
    } else if ( args[arg] == "-aniso-cal" ) {
      aniso = FCAL;
    } else if ( args[arg] == "-verbose" ) {
      if ( ++arg < args.size() ) verbose = clipper::String(args[arg]).i();
    } else {
      std::cout << "Unrecognized:\t" << args[arg] << "\n";
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "Usage: csfcalc\n\t-pdbin <filename>\n\t-mtzin <filename>\n\t-colin-fo <colpath>\n\t-colin-free <colpath>\n\t-mtzout <filename>\n\t-colout <colpath>\n\t-free-flag <free set>\n\t-num-reflns <reflns per spline param>\n\t-num-params <spline params>\n\t-no-bulk\n\t-aniso-obs\n\t-aniso-cal\nStructure factor calculation with bulk solvent correction.\n";
    exit(1);
  }

  // make data objects
  clipper::CCP4MTZfile mtzin, mtzout;
  clipper::MTZcrystal cxtl;
  clipper::HKL_info hkls;
  double bulkfrc, bulkscl;
  typedef clipper::HKL_data_base::HKL_reference_index HRI;
  using clipper::data32::F_sigF;  using clipper::data32::F_phi;
  using clipper::data32::Phi_fom; using clipper::data32::Flag;
  mtzin.set_column_label_mode( clipper::CCP4MTZfile::Legacy );

  // open file
  mtzin.open_read( ipfile );
  mtzin.import_hkl_info( hkls );
  mtzin.import_crystal( cxtl, ipcolfo+".F_sigF.F" );
  clipper::HKL_data<F_sigF> fo( hkls, cxtl );
  clipper::HKL_data<Flag> free( hkls, cxtl );
  mtzin.import_hkl_data( fo, ipcolfo );
  if ( ipcolfree != "NONE" ) mtzin.import_hkl_data( free, ipcolfree );
  if ( opcol[0] != '/' ) opcol = mtzin.assigned_paths()[0].notail()+"/"+opcol;
  mtzin.close_read();

  // atomic model
  clipper::MMDBManager mmdb;
  const int mmdbflags = ::mmdb::MMDBF_IgnoreBlankLines | ::mmdb::MMDBF_IgnoreDuplSeqNum | ::mmdb::MMDBF_IgnoreNonCoorPDBErrors | ::mmdb::MMDBF_IgnoreRemarks;
  mmdb.SetFlag( mmdbflags );
  mmdb.ReadPDBASCII( (char*)ippdb.c_str() );

  // get a list of all the atoms
  clipper::mmdb::PPCAtom psel;
  int hndl, nsel;
  hndl = mmdb.NewSelection();
  mmdb.SelectAtoms( hndl, 0, 0, ::mmdb::SKEY_NEW );
  mmdb.GetSelIndex( hndl, psel, nsel );
  clipper::MMDBAtom_list atoms( psel, nsel );
  mmdb.DeleteSelection( hndl );

  // calculate structure factors
  clipper::HKL_data<F_phi> fc( hkls, cxtl );
  if ( bulk ) {
    clipper::SFcalc_obs_bulk<float> sfcb;
    sfcb( fc, fo, atoms );
    bulkfrc = sfcb.bulk_frac();
    bulkscl = sfcb.bulk_scale();
  } else {
    clipper::SFcalc_aniso_fft<float> sfc;
    sfc( fc, atoms );
    bulkfrc = bulkscl = 0.0;
  }

  // do anisotropic scaling
  if ( aniso != NONE )  {
    clipper::SFscale_aniso<float>::TYPE F = clipper::SFscale_aniso<float>::F;
    clipper::SFscale_aniso<float> sfscl;
    if ( aniso == FOBS ) sfscl( fo, fc );  // scale Fobs
    if ( aniso == FCAL ) sfscl( fc, fo );  // scale Fcal
    std::cout << "\nAnisotropic scaling:\n"
	      << sfscl.u_aniso_orth(F).format() << "\n";
  }

  // now do sigmaa calc
  clipper::HKL_data<F_phi>   fb( hkls, cxtl ), fd( hkls, cxtl );
  clipper::HKL_data<Phi_fom> phiw( hkls, cxtl );
  clipper::HKL_data<Flag>    flag( hkls, cxtl );
  for ( HRI ih = flag.first(); !ih.last(); ih.next() )
    if ( !fo[ih].missing() && (free[ih].missing()||free[ih].flag()==freeflag) )
      flag[ih].flag() = clipper::SFweight_spline<float>::BOTH;
    else
      flag[ih].flag() = clipper::SFweight_spline<float>::NONE;

  // do sigmaa calc
  clipper::SFweight_spline<float> sfw( n_refln, n_param );
  sfw( fb, fd, phiw, fo, fc, flag );

  // calc abcd
  clipper::HKL_data<clipper::data32::ABCD> abcd( hkls );
  abcd.compute( phiw, clipper::data32::Compute_abcd_from_phifom() );

  // output data
  mtzout.open_append( ipfile, opfile );
  mtzout.export_hkl_data( abcd, opcol );
  mtzout.export_hkl_data( fc, opcol );
  mtzout.export_hkl_data( fb, opcol+"_BEST" );
  mtzout.export_hkl_data( fd, opcol+"_DIFF" );
  mtzout.close_append();

  // now calc R and R-free
  std::vector<double> params( n_param, 1.0 );
  clipper::BasisFn_spline basisfn( fo, n_param, 1.0 );
  clipper::TargetFn_scaleF1F2<F_phi,F_sigF> targetfn( fc, fo );
  clipper::ResolutionFn rfn( hkls, basisfn, targetfn, params );
  double r1w, f1w, r1f, f1f, Fo, Fc;
  r1w = f1w = r1f = f1f = 0.0;
  for ( HRI ih = fo.first(); !ih.last(); ih.next() )
    if ( !fo[ih].missing() ) {
      Fo = fo[ih].f();
      Fc = sqrt( rfn.f(ih) ) * fc[ih].f();
      if ( free[ih].flag() == freeflag ) {
	r1f += fabs( Fo - Fc );
	f1f += Fo;
      } else {
	r1w += fabs( Fo - Fc );
	f1w += Fo;
      }
    }
  r1f /= clipper::Util::max( f1f, 0.1 );
  r1w /= clipper::Util::max( f1w, 0.1 );
  std::cout << "\n R-factor      : " << r1w
	    << "\n Free R-factor : " << r1f << "\n";

  // DIAGNOSTIC OUTPUT
  if ( verbose > 1 ) {
    std::cout << "\n Bulk Correction Volume: " << bulkfrc;
    std::cout << "\n Bulk Correction Factor: " << bulkscl << "\n";
    std::cout << "\nNumber of spline params: " << sfw.params_scale().size() << "\n";
    clipper::BasisFn_spline basisfn( hkls, sfw.params_scale().size(), 1.0 );
    printf("\n $TABLE: Sigmaa statistics :\n $GRAPHS:scale vs resolution:N:1,2:\n        :lack of closure vs resolution:N:1,3:\n $$\n 1/resol^2   scale   lack_of_closure $$\n $$\n");
    for ( int i = 0; i <= 20.0; i++ ) {
      double s = hkls.resolution().invresolsq_limit()*double(i)/20.0;
      printf( "%6.3f %12.3f %12.3f\n",
	      s, basisfn.f_s(s,sfw.params_scale()), basisfn.f_s(s,sfw.params_error()) );
    }
    printf(" $$\n");
  }

}
