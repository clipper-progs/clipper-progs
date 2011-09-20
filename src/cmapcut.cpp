// Clipper mapcut
/* Copyright 2011 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-minimol.h>

#include <algorithm>

extern "C" {
#include <string.h>
#include <stdlib.h>
}


int main( int argc, char** argv )
{
  CCP4Program prog( "cmapcut", "0.1", "$Date: 2011/07/04" );

  std::cout << "\nCopyright 2011 Kevin Cowtan and University of York\n";
  std::cout << "All rights reserved. Please reference:\n";
  std::cout << " Cowtan K. (2011) 'mapcut software'.\n";

  // defaults
  clipper::String title;
  clipper::String ippdb_wrk = "NONE";
  clipper::String ipmap_wrk = "NONE";
  clipper::String ipmsk_wrk = "NONE";
  clipper::String ipmtz_wrk = "NONE";
  clipper::String opmsk = "NONE";
  clipper::String opmap = "NONE";
  clipper::String opmtz = "NONE";
  clipper::String ipcol_wrk_fc = "NONE";
  double radius = 5.0;         // mask radius
  double cellmu = 2.5;         // cell multiplier
  double resol = -1.0;         // resolution
  bool offset = false;         // shift density to origin
  bool respad = false;         // pad resolution by one reciprocal cell
  int verbose = 0;

  // command input
  CCP4CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if        ( args[arg] == "-title" ) {
      if ( ++arg < args.size() ) title = args[arg];
    } else if ( args[arg] == "-pdbin" ) {
      if ( ++arg < args.size() ) ippdb_wrk = args[arg];
    } else if ( args[arg] == "-mapin" ) {
      if ( ++arg < args.size() ) ipmap_wrk = args[arg];
    } else if ( args[arg] == "-mskin" ) {
      if ( ++arg < args.size() ) ipmsk_wrk = args[arg];
    } else if ( args[arg] == "-mtzin" ) {
      if ( ++arg < args.size() ) ipmtz_wrk = args[arg];
    } else if ( args[arg] == "-colin-fc" ) {
      if ( ++arg < args.size() ) ipcol_wrk_fc = args[arg];
    } else if ( args[arg] == "-mapout" ) {
      if ( ++arg < args.size() ) opmap = args[arg];
    } else if ( args[arg] == "-mskout" ) {
      if ( ++arg < args.size() ) opmsk = args[arg];
    } else if ( args[arg] == "-mtzout" ) {
      if ( ++arg < args.size() ) opmtz = args[arg];
    } else if ( args[arg] == "-mask-radius" ) {
      if ( ++arg < args.size() ) radius = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-cell-multiplier" ) {
      if ( ++arg < args.size() ) cellmu = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-resolution" ) {
      if ( ++arg < args.size() ) resol  = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-offset-to-origin" ) {
      offset = true;
    } else if ( args[arg] == "-pad-resolution" ) {
      respad = true;
    } else if ( args[arg] == "-verbose" ) {
      if ( ++arg < args.size() ) verbose = clipper::String(args[arg]).i();
    } else {
      std::cout << "\nUnrecognized:\t" << args[arg] << std::endl;
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "\nUsage: cmapcut\n\t-pdbin <filename>\n\t-mapin <filename>\n\t-mskin <filename>\n\t-mtzin <filename>\n\t-colin-fc <colpath>\n\t-mtzout <filename>\n\t-mapout <filename>\n\t-mskout <filename>\n\t-mask-radius <radius/A>       (default 5A)\n\t-cell-multiplier <multiplier> (default 2.5x)\n\t-resolution <resolution> (default as input)\n\t-offset-to-origin\n\t-pad-resolution\nCut electron density from a map and prepare a map or structure factors.\nNormal usage is to specify an MTZ (or map) to cut and a model (or mask)\nto cut with. The result is output as an MTZ or map." << std::endl;
    clipper::Message::message(clipper::Message_fatal("Invalid option."));
  }

  std::cout << std::endl;
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  typedef clipper::NXmap<float>::Map_reference_index NRI;
  using clipper::data32::F_phi;

  // numbers to output
  clipper::CCP4MTZfile mtzfile;
  mtzfile.set_column_label_mode( clipper::CCP4MTZfile::Legacy );
 
  // get work map
  clipper::Xmap<float> xmap;
  if ( ipmtz_wrk != "NONE" ) {  // from mtz
    // Get work reflection data
    clipper::HKL_data<F_phi> wrk_fp;
    mtzfile.open_read( ipmtz_wrk );
    mtzfile.import_hkl_data( wrk_fp, ipcol_wrk_fc );
    mtzfile.close_read();
    clipper::Grid_sampling grid( wrk_fp.hkl_info().spacegroup(),
				 wrk_fp.hkl_info().cell(),
				 wrk_fp.hkl_info().resolution() );
    xmap.init( wrk_fp.hkl_info().spacegroup(), wrk_fp.hkl_info().cell(), grid );
    xmap.fft_from( wrk_fp );
    if ( resol <= 0.0 ) resol = 1.0/sqrt(wrk_fp.invresolsq_range().max());
  } else {  // from map file
    clipper::CCP4MAPfile mapfile;
    mapfile.open_read( ipmap_wrk );
    mapfile.import_xmap( xmap );
    mapfile.close_read();
  }

  // map params
  clipper::Spacegroup    spgr = xmap.spacegroup();
  clipper::Cell          cell = xmap.cell();
  clipper::Grid_sampling grid = xmap.grid_sampling();

  // Get work xmsk
  clipper::NXmap<float> mask;
  if ( ippdb_wrk != "NONE" ) {  // from model
    clipper::MiniMol mol_tmp;
    clipper::MMDBfile mmdb;
    mmdb.read_file( ippdb_wrk );
    mmdb.import_minimol( mol_tmp );
    clipper::Atom_list atoms = mol_tmp.atom_list();
    clipper::Coord_map cm0(  1.0e20,  1.0e20,  1.0e20 );
    clipper::Coord_map cm1( -1.0e20, -1.0e20, -1.0e20 );
    for ( int i = 0; i < atoms.size(); i++ ) {
      const clipper::Coord_map cm = xmap.coord_map( atoms[i].coord_orth() );
      for ( int j = 0; j < 3; j++ ) {
	cm0[j] = std::min( cm0[j], cm[j] );
	cm1[j] = std::max( cm1[j], cm[j] );
      }
    }
    clipper::Grid_range gr0( cell, grid, radius+1.0 );
    clipper::Grid_range gr1( cm0.floor()+gr0.min(), cm1.ceil()+gr0.max() );
    mask.init( cell, grid, gr1 );
    clipper::EDcalc_mask<float> edcalc( radius );
    edcalc( mask, atoms );
  } else {  // from map file
    clipper::CCP4MAPfile mapfile;
    mapfile.open_read( ipmsk_wrk );
    mapfile.import_nxmap( mask );
    mapfile.close_read();
  }

  // find the bounds of the masked region
  clipper::Coord_orth co0(  1.0e20,  1.0e20,  1.0e20 );
  clipper::Coord_orth co1( -1.0e20, -1.0e20, -1.0e20 );
  clipper::Coord_orth com( 0.0, 0.0, 0.0 );
  double count = 0.0;
  for ( NRI ix = mask.first(); !ix.last(); ix.next() )
    if ( mask[ix] > 0.0 ) {
      const clipper::Coord_orth co = ix.coord_orth();
      for ( int j = 0; j < 3; j++ ) {
	co0[j] = std::min( co0[j], co[j] );
	co1[j] = std::max( co1[j], co[j] );
      }
      com += co;
      count += 1.0;
    }
  com = (1.0/count) * com;
  std::cout << "Centre of mass:" << std::endl
            << "  " << com.format() << std::endl;
  std::cout << "Masked region bounds:" << std::endl
            << "  " << co0.format() << std::endl
	    << "  " << co1.format() << std::endl;
  clipper::Coord_orth coff( 0.0, 0.0, 0.0 );
  if ( offset ) coff = com;

  // set cell
  const clipper::Coord_orth cd = cellmu * ( co1 - co0 );
  const clipper::Cell cellcut( clipper::Cell_descr( cd[0], cd[1], cd[2] ) );
  const clipper::Spacegroup spgrcut( clipper::Spacegroup::P1 );
  if ( resol <= 0.0 ) {  // get cell from map
    resol = 1.0e20;
    const clipper::Cell&          c = xmap.cell();
    const clipper::Grid_sampling& g = xmap.grid_sampling();
    resol = std::min( resol, 2.0/(c.a_star()*double(g.nu())) );
    resol = std::min( resol, 2.0/(c.b_star()*double(g.nv())) );
    resol = std::min( resol, 2.0/(c.c_star()*double(g.nw())) );
    resol *= 1.5;
  }
  if ( respad ) {  // pad resolution by one grid cell
    resol = 1.0 / ( 1.0/resol + sqrt( pow(cellcut.a(),-2.0) +
				      pow(cellcut.b(),-2.0) +
				      pow(cellcut.c(),-2.0) ) );
  }
  clipper::Resolution rescut( resol );
  clipper::Grid_sampling gridcut( spgrcut, cellcut, rescut );
  // reduce resolution if grid too large
  while ( double(gridcut.size()) > 0.5e9 ) {
    std::cout << "Warning: reducing resolution " << resol << " A" << std::endl;
    resol *= pow( double(gridcut.size())/0.4e9, 0.333 );
    rescut = clipper::Resolution( resol );
    gridcut = clipper::Grid_sampling( spgrcut, cellcut, rescut );
  }
  // output
  std::cout << "Cut map spacegroup: " << spgrcut.symbol_hm() << std::endl;
  std::cout << "Cut map cell      : " << cellcut.format() << std::endl;
  std::cout << "Cut map sampling  : " << gridcut.format() << std::endl;
  std::cout << "Cut map resolution: " << rescut.limit() << " A" << std::endl;

  // make the cut map
  clipper::Xmap<float> xcut( spgrcut, cellcut, gridcut );
  xcut = 0.0;

  // now find the bounds of the masked region on the cut map coords
  clipper::Coord_grid cg0(  999999,  999999,  999999 );
  clipper::Coord_grid cg1( -999999, -999999, -999999 );
  for ( NRI ix = mask.first(); !ix.last(); ix.next() )
    if ( mask[ix] > 0.0 ) {
      clipper::Coord_grid cg =
	xcut.coord_map(ix.coord_orth()-coff).coord_grid();
      for ( int j = 0; j < 3; j++ ) {
	cg0[j] = std::min( cg0[j], cg[j]-1 );
	cg1[j] = std::max( cg1[j], cg[j]+1 );
      }
    }

  // fill the cut map
  clipper::Xmap_base::Map_reference_coord i0( xcut, cg0 ), iu, iv, iw;
  for ( iu = i0; iu.coord().u() <= cg1.u(); iu.next_u() )
    for ( iv = iu; iv.coord().v() <= cg1.v(); iv.next_v() )
      for ( iw = iv; iw.coord().w() <= cg1.w(); iw.next_w() ) {
	const clipper::Coord_orth co = iw.coord_orth()+coff;
	clipper::Coord_grid cg = mask.coord_map(co).coord_grid();
	if ( mask.in_map( cg ) )
	  if ( mask.get_data( cg ) > 0.0 )
	    xcut[iw] = xmap.interp<clipper::Interp_cubic>(xmap.coord_map(co));
      }

  // make maps
  if ( opmap != "NONE" ) {
    clipper::CCP4MAPfile mapfile;
    mapfile.open_write( opmap );
    mapfile.export_xmap( xcut );
    mapfile.close_write();
  }
  if ( opmsk != "NONE" ) {
    clipper::CCP4MAPfile mapfile;
    mapfile.open_write( opmsk );
    mapfile.set_cell( cell );
    mapfile.export_nxmap( mask );
    mapfile.close_write();
  }

  // make structure factors
  if ( opmtz != "NONE" ) {
    clipper::HKL_sampling hklcut( cellcut, rescut );
    clipper::HKL_data<F_phi> fpcut( spgrcut, cellcut, hklcut );
    xcut.fft_to( fpcut );
    std::cout << "Number of reflections: " << fpcut.num_obs() << std::endl;
    clipper::CCP4MTZfile mtzfile;
    mtzfile.open_write( opmtz );
    mtzfile.export_hkl_info( fpcut.hkl_info() );
    mtzfile.export_hkl_data( fpcut, "/*/*/mapcut" );
    mtzfile.close_write();
  }
}
