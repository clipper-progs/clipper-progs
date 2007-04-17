// Clipper app to do E calc
/* Copyright 2003-2004 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>


int main( int argc, char** argv )
{
  CCP4Program prog( "cecalc", "0.1", "$Date: 2004/07/01" );

  // defaults
  clipper::String title;
  clipper::String ipfile = "NONE";
  clipper::String ipcolfo = "NONE";
  clipper::String opfile = "ecalc.mtz";
  clipper::String opcol = "ecalc";
  const int nprm = 12;

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
    } else if ( args[arg] == "-colin-fo" ) {
      if ( ++arg < args.size() ) ipcolfo = args[arg];
    } else if ( args[arg] == "-colout" ) {
      if ( ++arg < args.size() ) opcol = args[arg];
    } else {
      std::cout << "Unrecognized:\t" << args[arg] << "\n";
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "Usage: cecalc\n\t-mtzin <filename>\n\t-mtzout <filename>\n\t-colin-fo <colpath>\n\t-colout <colpath>\nCalculate E's from F's\n";
    exit(1);
  }

  // make data objects
  clipper::CCP4MTZfile mtzin, mtzout;
  clipper::HKL_info hkls;
  typedef clipper::HKL_data_base::HKL_reference_index HRI;

  // open file
  mtzin.open_read( ipfile );
  mtzin.import_hkl_info( hkls );
  clipper::HKL_data<clipper::data32::F_sigF> fsig( hkls );
  mtzin.import_hkl_data( fsig, ipcolfo );
  if ( opcol[0] != '/' ) opcol = mtzin.assigned_paths()[0].notail()+"/"+opcol;
  mtzin.close_read();

  // create initial E
  clipper::HKL_data<clipper::data32::E_sigE> esig( hkls );
  esig.compute( fsig, clipper::data32::Compute_EsigE_from_FsigF() );

  // calc E-scaling
  std::vector<double> params_init( nprm, 1.0 );
  clipper::BasisFn_spline basis_fo( esig, nprm, 2.0 );
  clipper::TargetFn_scaleEsq<clipper::data32::E_sigE> target_fo( esig );
  clipper::ResolutionFn escale( hkls, basis_fo, target_fo, params_init );

  // apply E-scaling
  for ( HRI ih = esig.first(); !ih.last(); ih.next() )
    if ( !esig[ih].missing() ) esig[ih].scale( sqrt( escale.f(ih) ) );

  // output data
  mtzout.open_append( ipfile, opfile );
  mtzout.export_hkl_data( esig, opcol );
  mtzout.close_append();


  // generate stats:
  double na, nc, sa, sc;
  na = nc = sa = sc = 0.0;
  for ( HRI ih = esig.first(); !ih.last(); ih.next() )
    if ( !esig[ih].missing() ) 
      if ( ih.hkl_class().centric() ) {
	nc += 1.0;
	sc += esig[ih].E()*esig[ih].E();
      } else {
	na += 1.0;
	sa += esig[ih].E()*esig[ih].E();
      }
  std::cout << "Number of reflections: " << clipper::Util::intr(na+nc) << "  Mean E^2: " << (sa+sc)/clipper::Util::max(na+nc,1.0) << std::endl;
  std::cout << "Number of acentrics  : " << clipper::Util::intr(na) << "  Mean E^2: " << (sa)/clipper::Util::max(na,1.0) << std::endl;
  std::cout << "Number of centrics   : " << clipper::Util::intr(nc) << "  Mean E^2: " << (sc)/clipper::Util::max(nc,1.0) << std::endl;
}
