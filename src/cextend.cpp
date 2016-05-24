#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-ccp4.h>
#include <iostream>
//gcc needs this 
#include <stdlib.h>

int main( int argc, char** argv )
{
  //defaults
  clipper::String ipfile = "NONE";
  clipper::String ipcolfo = "";
  clipper::String ipcolfiso = "";
  clipper::String opfile = "NONE";
  clipper::String opcolfo = "*/*/[F,SIGF]";
  clipper::String opcolfiso = "*/*/[F_ISO,SIGF_ISO]";
  clipper::String opcole = "*/*/[E,SIGE]";
  clipper::String opcoleiso = "*/*/[E_ISO,SIGE_ISO]";
  clipper::Resolution reso;
  bool aniso = false;
  double reslim = 1.0;
  double resfilt = -1.0;
  const int nparm = 12; // as in cecalc
  
  // command input
  CCP4CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if ( args[arg] == "-mtzin" ) {
      if ( ++arg < args.size() ) ipfile = args[arg];
    } else if ( args[arg] == "-mtzout" ) {
      if ( ++arg < args.size() ) opfile = args[arg];
    } else if ( args[arg] == "-colin-fo" ) {
      if ( ++arg < args.size() ) ipcolfo = args[arg];
    } else if ( args[arg] == "-colin-fiso" ) {
      if ( ++arg < args.size() ) ipcolfiso = args[arg];
    } else if ( args[arg] == "-resolution" ) {
      if ( ++arg < args.size() ) reslim = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-aniso" ) {
      aniso = true;
    }  else {
      std::cout << "Unrecognized:\t" << args[arg] << "\n";
      args.clear();
    }
  }
    if ( args.size() <= 1 ) {
    std::cout << "Usage: cextend\n\t-mtzin <filename>\n\t-mtzout <filename>\n\t-colin-fo <colpath>\n\t[ optional -colin-fiso <colpath> ]\n\t [ optional -resolution <resolution> default 1.0 ] \
      \n\t-aniso anisotopy correct \nExtend hkls to higher resolution, optionally correct for aniostropy and calcuate normalised structure factors\n";
    exit(1);
  }
  
  // make data objects
  clipper::Spacegroup sg;
  clipper::Cell cell;
  clipper::CCP4MTZfile mtzin, mtzout;
  clipper::MTZcrystal xtal;
  clipper::MTZdataset dset;
  typedef clipper::HKL_data_base::HKL_reference_index HRI;

  // read  file
  mtzin.open_read( ipfile );
  
  // spacegroup and cell from first file
  sg = mtzin.spacegroup();
  cell = mtzin.cell();
  
  // extend resolution to requested resolution 
  reso = clipper::Resolution( reslim );
  clipper::HKL_info hkls( sg, cell, reso, true );
  
  // import F,SIGF and Fiso,SIGFiso
  clipper::HKL_data<clipper::data32::F_sigF> fsig( hkls );
  clipper::HKL_data<clipper::data32::F_sigF> fsigiso( hkls );
  mtzin.import_hkl_data( fsig, ipcolfo );
  
  // if anisotropy corrected data in input file read them in
  if ( ipcolfiso != "" ) {
    mtzin.import_hkl_data( fsigiso, ipcolfiso ); 
  }
  
  // otherwise read uncorrected Fs and correct them if requested
  else {
    mtzin.import_hkl_data( fsigiso, ipcolfo ); 
  }
  mtzin.close_read();

  // now correct for anisotropy 
  if (aniso) {
    clipper::SFscale_aniso<float>::TYPE F = clipper::SFscale_aniso<float>::F;
    clipper::SFscale_aniso<float>::MODE M = clipper::SFscale_aniso<float>::NORMAL;    
    clipper::SFscale_aniso<float> sfscale ( 3.0, M ); // rejection criterion for F/sigF from caniso
    sfscale( fsigiso, resfilt, 12 );
    std::cout << "\n Correcting for anisotropy. Scale factors: \n" << sfscale.u_aniso_orth(F).format() << "\n";
  }
  
  // Fs to Es
  clipper::HKL_data<clipper::data32::E_sigE> esig ( hkls );
  if (aniso) {
    esig.compute( fsigiso, clipper::data32::Compute_EsigE_from_FsigF() );
  } else {
    esig.compute( fsig, clipper::data32::Compute_EsigE_from_FsigF() );
  }
  
  // now calculate scaling
  std::vector<double> initial_params( nparm, 1.0 );
  clipper::BasisFn_spline basis_f( esig, nparm, 2.0 );
  clipper::TargetFn_scaleEsq<clipper::data32::E_sigE> target_f ( esig );
  clipper::ResolutionFn escale ( hkls, basis_f, target_f, initial_params );
  
  // apply scaling
  for ( HRI ih = esig.first(); !ih.last(); ih.next() )
    if ( !esig[ih].missing() ) esig[ih].scale( sqrt( escale.f(ih) ) );
  
  // write output
  mtzout.open_write( opfile );
  mtzout.export_crystal( xtal, opcolfo );
  mtzout.export_dataset( dset, opcolfo );
  mtzout.export_hkl_info( fsig.hkl_info() );
  mtzout.export_hkl_data( fsig, opcolfo );
  if ( ( ipcolfiso != "" ) || (aniso) ) {
    mtzout.export_hkl_data( fsigiso, opcolfiso );
  }
  if (aniso) {
    mtzout.export_hkl_data ( esig, opcoleiso );
  } else {
    mtzout.export_hkl_data ( esig, opcole );
  }
  mtzout.close_write();
}