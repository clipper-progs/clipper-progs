// Clipper app to perform ffts and map stats
/* Copyright 2003-2004 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-ccp4.h>

extern "C" {
  #include <stdlib.h>
}



int main( int argc, char** argv )
{
  CCP4Program prog( "cmapcoeff", "0.1", "$Date: 2014/12/05" );

  // defaults
  clipper::String title;
  clipper::String ipfile = "NONE";
  clipper::String ipcola = "NONE";
  clipper::String ipcolf1 = "NONE";
  clipper::String ipcolf2 = "NONE";
  clipper::String ipcolh1 = "NONE";
  clipper::String ipcolh2 = "NONE";
  clipper::String ipcolp1 = "NONE";
  clipper::String ipcolp2 = "NONE";
  clipper::String opfile = "NONE";
  clipper::String opcol  = "NONE";
  bool scale1 = false;
  bool scale2 = false;
  double uvalue = 0.0;
  clipper::Resolution reso;

  // command input
  CCP4CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if ( args[arg] == "-title" ) {
      if ( ++arg < args.size() ) title = args[arg];
    } else if ( args[arg] == "-mtzin" ) {
      if ( ++arg < args.size() ) ipfile = args[arg];
    } else if ( args[arg] == "-mtzout" ) {
      if ( ++arg < args.size() ) opfile = args[arg];
    } else if ( args[arg] == "-colin-fano" ) {
      if ( ++arg < args.size() ) ipcola = args[arg];
    } else if ( args[arg] == "-colin-fo-1" ) {
      if ( ++arg < args.size() ) ipcolf1 = args[arg];
    } else if ( args[arg] == "-colin-fo-2" ) {
      if ( ++arg < args.size() ) ipcolf2 = args[arg];
    } else if ( args[arg] == "-colin-hl-1" ) {
      if ( ++arg < args.size() ) ipcolh1 = args[arg];
    } else if ( args[arg] == "-colin-hl-2" ) {
      if ( ++arg < args.size() ) ipcolh2 = args[arg];
    } else if ( args[arg] == "-colin-phifom-1" ) {
      if ( ++arg < args.size() ) ipcolp1 = args[arg];
    } else if ( args[arg] == "-colin-phifom-2" ) {
      if ( ++arg < args.size() ) ipcolp2 = args[arg];
    } else if ( args[arg] == "-colout" ) {
      if ( ++arg < args.size() ) opcol = args[arg];
    } else if ( args[arg] == "-scale-fo-1" ) {
      scale1 = true;
    } else if ( args[arg] == "-scale-fo-2" ) {
      scale2 = true;
    } else if ( args[arg] == "-u-value" ) {
      if ( ++arg < args.size() ) uvalue = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-b-value" ) {
      if ( ++arg < args.size() ) uvalue = clipper::Util::b2u(clipper::String(args[arg]).f());
   } else if ( args[arg] == "-resolution" ) {
      if ( ++arg < args.size() ) {
        reso = clipper::Resolution( clipper::String(args[arg]).f() );
      }
    } else {
      std::cout << "Unrecognized:\t" << args[arg] << "\n";
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "Usage: cmapcoeff\n\t-mtzin <filename>\n\t-mtzout <filename>\n\t-colin-fo-1 <colpath>\n\t-colin-fo-2 <colpath>\n\t-colin-fano <colpath>\n\t-colin-hl-1 <colpath>\n\t-colin-hl-2 <colpath>\n\t-colin-phi-fom-1 <colpath>\n\t-colin-phi-fom-2 <colpath>\n\t-scale-fo-1\n\t-scale-fo-2\n\t-u-value <U>\n\t-b-value <B>\n\t-resolution <reso>\nIf -colin-fano is supplied, an anomalous difference map is calculated. Otherwise a weighted, difference or double-difference map is calculated. TODO: SCALING AND OUTLIER REJECTION\n";
    exit(1);
  }

  // make data objects
  clipper::CCP4MTZfile mtzin, mtzout;
  clipper::HKL_info hkls;
  typedef clipper::HKL_data_base::HKL_reference_index HRI;
  mtzin.set_column_label_mode( clipper::CCP4MTZfile::Legacy );

  // open file
  mtzin.open_read( ipfile );
  if ( reso.is_null() ) reso = mtzin.resolution();
  mtzin.import_hkl_info( hkls );
  clipper::HKL_data<clipper::data32::F_sigF_ano> fano( hkls );
  clipper::HKL_data<clipper::data32::F_sigF>     fsig1( hkls ), fsig2( hkls );
  clipper::HKL_data<clipper::data32::ABCD>       abcd1( hkls ), abcd2( hkls );
  clipper::HKL_data<clipper::data32::Phi_fom>    phiw1( hkls ), phiw2( hkls );
  if ( ipcola != "NONE" ) mtzin.import_hkl_data( fano, ipcola );
  if ( ipcolf1 != "NONE" ) mtzin.import_hkl_data( fsig1, ipcolf1 );
  if ( ipcolf2 != "NONE" ) mtzin.import_hkl_data( fsig2, ipcolf2 );
  if ( ipcolh1 != "NONE" ) mtzin.import_hkl_data( abcd1, ipcolh1 );
  if ( ipcolh2 != "NONE" ) mtzin.import_hkl_data( abcd2, ipcolh2 );
  if ( ipcolp1 != "NONE" ) mtzin.import_hkl_data( phiw1, ipcolp1 );
  if ( ipcolp2 != "NONE" ) mtzin.import_hkl_data( phiw2, ipcolp2 );
  if ( opcol[0] != '/' ) opcol = mtzin.assigned_paths()[0].notail()+"/"+opcol;
  mtzin.close_read();

  // make phases if necessary
  if ( ipcolp1 == "NONE" && ipcolh1 != "NONE" )
    phiw1.compute( abcd1, clipper::data32::Compute_phifom_from_abcd() );
  if ( ipcolp2 == "NONE" && ipcolh2 != "NONE" )
    phiw2.compute( abcd2, clipper::data32::Compute_phifom_from_abcd() );

  // scaling
  if ( ipcolf1 != "NONE" && ipcolf2 != "NONE" && (scale1||scale2) ) {
    const int n_param = 10;
    std::vector<double> params( n_param, 1.0 );
    if ( scale1 ) {
      clipper::BasisFn_spline basisfn( fsig1, n_param, 1.0 );
      clipper::TargetFn_scaleF1F2<clipper::data32::F_sigF,clipper::data32::F_sigF> targetfn( fsig1, fsig2 );
      clipper::ResolutionFn rfn( hkls, basisfn, targetfn, params );
      for ( HRI ih = fsig1.first(); !ih.last(); ih.next() )
        if ( !fsig1[ih].missing() ) fsig1[ih].scale(sqrt(rfn.f(ih)));
    } else {
      clipper::BasisFn_spline basisfn( fsig2, n_param, 1.0 );
      clipper::TargetFn_scaleF1F2<clipper::data32::F_sigF,clipper::data32::F_sigF> targetfn( fsig2, fsig1 );
      clipper::ResolutionFn rfn( hkls, basisfn, targetfn, params );
      for ( HRI ih = fsig1.first(); !ih.last(); ih.next() )
        if ( !fsig2[ih].missing() ) fsig2[ih].scale(sqrt(rfn.f(ih)));      
    }
  }

  // map coefficients
  clipper::HKL_data<clipper::data32::F_phi> fphi1( hkls ), fphi2( hkls );
  if ( ipcola != "NONE" ) {
    // anomalous difference F
    fsig1.compute( fano, clipper::data32::Compute_diff_fsigf_from_fsigfano() );
    for ( HRI ih = fphi1.first(); !ih.last(); ih.next() )
      fphi1[ih].shift_phase( -0.5*clipper::Util::pi() );
    fphi1.compute( fsig1, phiw1, clipper::data32::Compute_fphi_from_fsigf_phifom() );
  } else {
    fphi1.compute( fsig1, phiw1, clipper::data32::Compute_fphi_from_fsigf_phifom() );
    if ( ipcolf2 != "NONE" ) {
      if ( ipcolp2 != "NONE" || ipcolh2 != "NONE" )
        fphi2.compute( fsig2, phiw2, clipper::data32::Compute_fphi_from_fsigf_phifom() );
      else
        fphi2.compute( fsig2, phiw1, clipper::data32::Compute_fphi_from_fsigf_phifom() );
      fphi1.compute( fphi1, fphi2, clipper::data32::Compute_sub_fphi() );
    }
  }

  // apply U value
  fphi1.compute( fphi1, clipper::data32::Compute_scale_u_iso_fphi(1.0,-uvalue) );

  // write to file
  mtzout.open_append( ipfile, opfile );
  mtzout.export_hkl_data( fphi1, opcol );
  mtzout.close_append();
}
