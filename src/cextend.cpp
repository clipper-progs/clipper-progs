//C (C) 2015-2016 Huw Jenkins, The University of York
//L
//L  This library is free software and is distributed under the terms
//L  and conditions of version 2.1 of the GNU Lesser General Public
//L  Licence (LGPL) with the following additional clause:
//L
//L     `You may also combine or link a "work that uses the Library" to
//L     produce a work containing portions of the Library, and distribute
//L     that work under terms of your choice, provided that you give
//L     prominent notice with each copy of the work that the specified
//L     version of the Library is used in it, and that you include or
//L     provide public access to the complete corresponding
//L     machine-readable source code for the Library including whatever
//L     changes were used in the work. (i.e. If you make changes to the
//L     Library you must distribute those, but you do not need to
//L     distribute source or object code to those portions of the work
//L     not covered by this licence.)'
//L
//L  Note that this clause grants an additional right and does not impose
//L  any additional restriction, and so does not affect compatibility
//L  with the GNU General Public Licence (GPL). If you wish to negotiate
//L  other terms, please contact the maintainer.
//L
//L  You can redistribute it and/or modify the library under the terms of
//L  the GNU Lesser General Public License as published by the Free Software
//L  Foundation; either version 2.1 of the License, or (at your option) any
//L  later version.
//L
//L  This library is distributed in the hope that it will be useful, but
//L  WITHOUT ANY WARRANTY; without even the implied warranty of
//L  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//L  Lesser General Public License for more details.
//L
//L  You should have received a copy of the CCP4 licence and/or GNU
//L  Lesser General Public License along with this library; if not, write
//L  to the CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//L  The GNU Lesser General Public can also be obtained by writing to the
//L  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//L  MA 02111-1307 USA

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
  clipper::String ipcolphifom = "";
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
    } else if ( args[arg] == "-colin-phifom" ) {
      if ( ++arg < args.size() ) ipcolphifom = args[arg];
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
    std::cout << "Usage: cextend\n\t-mtzin <filename>\n\t-mtzout <filename>\n\t-colin-fo <colpath>\n\t[ optional -colin-fiso <colpath> ]\
      \n\t[ optional -colin-phifom <colpath> ]\n\t[ optional -resolution <resolution> default 1.0 ] \
      \n\t-aniso anisotropy correct \nExtend hkls to higher resolution, optionally correct for anisotropy and calculate normalised structure factors\n";
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
  clipper::HKL_data<clipper::data32::Phi_fom> phifom( hkls );
  mtzin.import_hkl_data( fsig, ipcolfo );
  
  // if anisotropy corrected data in input file read them in
  if ( ipcolfiso != "" ) {
    mtzin.import_hkl_data( fsigiso, ipcolfiso ); 
  }
  // otherwise read uncorrected Fs and correct them if requested
  else {
    mtzin.import_hkl_data( fsigiso, ipcolfo ); 
  }
  
  // if phases are supplied read them in too
  if ( ipcolphifom != "" ) {
    mtzin.import_hkl_data( phifom, ipcolphifom );
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
  if ( ipcolphifom != "") {
    mtzout.export_hkl_data ( phifom, ipcolphifom );
  }
  mtzout.close_write();
}