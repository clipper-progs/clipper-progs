// Clipper app to move one model to match another using xtal symmetry
/* Copyright 2007 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-minimol.h>


class MapFilterFn_g5 : public clipper::MapFilterFn_base {
public: clipper::ftype operator() ( const clipper::ftype& radius ) const { return exp(-radius*radius/50.0); }
};


int main( int argc, char** argv )
{
  CCP4Program prog( "csymmatch", "0.2", "$Date: 2007/11/06" );

  // defaults
  clipper::String title;
  clipper::String ippdbref = "NONE";
  clipper::String ippdb    = "NONE";
  clipper::String oppdb    = "symmatch.pdb";

  // command input
  CCP4CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if ( args[arg] == "-title" ) {
      if ( ++arg < args.size() ) title = args[arg];
    } else if ( args[arg] == "-pdbin-ref" ) {
      if ( ++arg < args.size() ) ippdbref = args[arg];
    } else if ( args[arg] == "-pdbin" ) {
      if ( ++arg < args.size() ) ippdb = args[arg];
    } else if ( args[arg] == "-pdbout" ) {
      if ( ++arg < args.size() ) oppdb = args[arg];
    } else {
      std::cout << "Unrecognized:\t" << args[arg] << "\n";
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "Usage: csymmatch\n\t-pdbin-ref <filename>\n\t-pdbin <filename>\n\t-pdbout <filename>\nApply symmetry and cell shifts to each chain in 'pdbin' to obtain the best match to 'pdbin-ref'.\n";
    exit(1);
  }

  // atomic models
  const int mmdbflags = MMDBF_IgnoreBlankLines | MMDBF_IgnoreDuplSeqNum | MMDBF_IgnoreNonCoorPDBErrors | MMDBF_IgnoreRemarks;
  clipper::MMDBfile mmdbref, mmdbwrk;
  clipper::MiniMol molref, molwrk;
  mmdbref.SetFlag( mmdbflags );
  mmdbwrk.SetFlag( mmdbflags );
  mmdbref.read_file( ippdbref );
  mmdbwrk.read_file( ippdb    );
  mmdbref.import_minimol( molref );
  mmdbwrk.import_minimol( molwrk );

  clipper::Spacegroup spg1 = clipper::Spacegroup(clipper::Spacegroup::P1);
  clipper::Spacegroup spgr = mmdbwrk.spacegroup();
  clipper::Cell       cell = mmdbwrk.cell();
  if ( spgr.is_null() ) {
    std::cerr << "Cannot get spacegroup, check PDB file and CCP4 setup"
	      << std::endl;
    exit(1);
  }

  // user output
  std::cout << std::endl << "Reference molecule:" << std::endl;
  if ( !mmdbref.spacegroup().is_null() )
    std::cout << "  Spacegroup: " 
	      << mmdbref.spacegroup().symbol_hm() << std::endl;
  if ( !mmdbref.cell().is_null() )
    std::cout << "  Unit cell: " 
	      << mmdbref.cell().format() << std::endl;
  std::cout << "Moving molecule:" << std::endl;
  if ( !mmdbwrk.spacegroup().is_null() )
    std::cout << "  Spacegroup: " 
	      << mmdbwrk.spacegroup().symbol_hm() << std::endl;
  if ( !mmdbwrk.cell().is_null() )
    std::cout << "  Unit cell: " 
	      << mmdbwrk.cell().format() << std::endl;
  std::cout << std::endl;

  // calculate extent of model
  clipper::Atom_list atomr = molref.atom_list();
  clipper::Range<clipper::ftype> urange, vrange, wrange;
  clipper::Coord_frac cfr( 0.0, 0.0, 0.0 );
  for ( int i = 0; i < atomr.size(); i++ ) {
    clipper::Coord_frac cf = atomr[i].coord_orth().coord_frac( cell );
    cfr += cf;
    urange.include( cf.u() );
    vrange.include( cf.v() );
    wrange.include( cf.w() );
  }
  clipper::Coord_frac cf0( urange.min(), vrange.min(), wrange.min() );
  clipper::Coord_frac cf1( urange.max(), vrange.max(), wrange.max() );
  cfr = (1.0/double(atomr.size())) * cfr;

  // calculate mask using wrk cell and ref atoms
  clipper::Resolution reso( 5.0 );
  clipper::Grid_sampling grid( spg1, cell, reso );
  clipper::Grid_range    grng( grid,  cf0,  cf1 );
  grng.add_border(4);
  clipper::NXmap<float> nxmap( cell, grid, grng ), nxflt( cell, grid, grng );
  clipper::EDcalc_mask<float> maskcalc( 2.0 );
  nxmap = 0.0;
  maskcalc( nxmap, atomr );
  MapFilterFn_g5 fn;
  clipper::MapFilter_fft<float>
    fltr( fn, 1.0, clipper::MapFilter_fft<float>::Relative );
  fltr( nxflt, nxmap );

  // now score each chain, symmetry and offset in turn
  for ( int c = 0; c < molwrk.size(); c++ ) {
    double              bestscr = 0.0;
    int                 bestsym = 0;
    clipper::Coord_frac bestoff( 0.0, 0.0, 0.0 );
    const clipper::Coord_frac cfh( 0.5, 0.5, 0.5 );
    for ( int sym = 0; sym < spgr.num_symops(); sym++ ) {
      clipper::Atom_list atomw = molwrk[c].atom_list();
      clipper::RTop_orth rtop = spgr.symop(sym).rtop_orth( cell );
      clipper::Coord_orth cow( 0.0, 0.0, 0.0 );
      for ( int a = 0; a < atomw.size(); a++ ) {
	atomw[a].transform( rtop );
	cow += atomw[a].coord_orth();
      }
      if ( atomw.size() > 0 ) cow = (1.0/double(atomw.size())) * cow;
      clipper::Coord_frac cfw = cow.coord_frac( cell );
      clipper::Coord_frac cfwt = cfw.lattice_copy_near( cfr - cfh );
      clipper::Coord_frac off0 = cfwt - cfw;

      // try offsets
      for ( double du = 0.0; du <= 1.01; du += 1.0 )
	for ( double dv = 0.0; dv < 1.01; dv += 1.0 )
	  for ( double dw = 0.0; dw < 1.01; dw += 1.0 ) {
	    clipper::Coord_frac off( round( off0.u() ) + du,
				     round( off0.v() ) + dv,
				     round( off0.w() ) + dw );
	    clipper::Coord_orth ofo = off.coord_orth( cell );
	    double scr = 0.0;
	    for ( int a = 0; a < atomw.size(); a++ ) {
	      clipper::Coord_orth coa = atomw[a].coord_orth() + ofo;
	      clipper::Coord_grid cga = nxflt.coord_map( coa ).coord_grid();
	      if ( nxflt.in_map( cga ) ) scr += nxflt.get_data( cga );
	    }
	    if ( scr > bestscr ) {
	      bestscr = scr;
	      bestsym = sym;
	      bestoff = off;
	    }
	  }
    }
    // now transform using the best operator
    clipper::Coord_orth cot = bestoff.coord_orth( cell );
    clipper::RTop_orth rtop = spgr.symop(bestsym).rtop_orth( cell );
    rtop = clipper::RTop_orth( rtop.rot(), rtop.trn()+cot );
    molwrk[c].transform( rtop );

    // user output
    std::cout << "Chain: " << molwrk[c].id()
	      << " will be trasformed as follows:" << std::endl
	      << "  Symmetry operator: "
	      << spgr.symop(bestsym).format() << std::endl
	      << "  Lattice shift:     "
	      << bestoff.format() << std::endl;
  }

  // write file
  mmdbwrk.export_minimol( molwrk );
  mmdbwrk.write_file( oppdb );
}
