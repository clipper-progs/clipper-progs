// Clipper app to perform structure factor calculation
/* Copyright 2003-2004 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-ccp4.h>

extern "C" {
#include <stdlib.h>
}
 

int main( int argc, char** argv )
{
  CCP4Program prog( "caniso", "0.1", "$Date: 2007/03/01" );

  // defaults
  enum ANISO { NONE, FOBS, FCAL };
  clipper::String title;
  clipper::String ipfile = "NONE";
  clipper::String ipcolfo = "NONE";
  clipper::String opfile = "aniso.mtz";
  clipper::String opcol = "aniso";

  // command input
  CCP4CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if ( args[arg] == "-title" ) {
      if ( ++arg < args.size() ) title = args[arg];
    } else if ( args[arg] == "-mtzin" ) {
      if ( ++arg < args.size() ) ipfile = args[arg];
    } else if ( args[arg] == "-colin-fo" ) {
      if ( ++arg < args.size() ) ipcolfo = args[arg];
    } else if ( args[arg] == "-mtzout" ) {
      if ( ++arg < args.size() ) opfile = args[arg];
    } else if ( args[arg] == "-colout" ) {
      if ( ++arg < args.size() ) opcol = args[arg];
    } else {
      std::cout << "Unrecognized:\t" << args[arg] << "\n";
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "Usage: caniso\n\t-mtzin <filename>\n\t-colin-fo <colpath>\n\t-mtzout <filename>\n\t-colout <colpath>\nAnisotropy correction.\n";
    exit(1);
  }

  // make data objects
  clipper::CCP4MTZfile mtzin, mtzout;
  clipper::MTZcrystal cxtl;
  clipper::HKL_info hkls;
  typedef clipper::HKL_data_base::HKL_reference_index HRI;
  using clipper::data32::F_sigF;

  // open file
  mtzin.open_read( ipfile );
  mtzin.import_hkl_info( hkls );
  mtzin.import_crystal( cxtl, ipcolfo );
  clipper::HKL_data<F_sigF> fo( hkls, cxtl );
  mtzin.import_hkl_data( fo, ipcolfo );
  if ( opcol[0] != '/' ) opcol = mtzin.assigned_paths()[0].notail()+"/"+opcol;
  mtzin.close_read();

  // scale structure factors
  clipper::SFscale_aniso<float>::TYPE F = clipper::SFscale_aniso<float>::F;
  clipper::SFscale_aniso<float> sfscl( 3.0 );
  sfscl( fo );
  std::cout << "\nAnisotropic scaling:\n"
	    << sfscl.u_aniso_orth(F).format() << "\n";

  // output data
  mtzout.open_append( ipfile, opfile );
  mtzout.export_hkl_data( fo, opcol );
  mtzout.close_append();
}
