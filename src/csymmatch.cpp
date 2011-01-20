// Clipper app to move one model to match another using xtal symmetry
/* Copyright 2007 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-minimol.h>

extern "C" {
  #include <stdlib.h>
}


class MapFilterFn_g5 : public clipper::MapFilterFn_base {
public: clipper::ftype operator() ( const clipper::ftype& radius ) const { return exp(-radius*radius/50.0); }
};


int main( int argc, char** argv )
{
  CCP4Program prog( "csymmatch", "0.3", "$Date: 2009/07/13" );

  // defaults
  clipper::String title;
  clipper::String ippdbref = "NONE";
  clipper::String ippdb    = "NONE";
  clipper::String oppdb    = "symmatch.pdb";
  double crad = 2.0;
  bool omatch = false;

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
    } else if ( args[arg] == "-connectivity-radius" ) {
      if ( ++arg < args.size() ) crad = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-origin-hand" ) {
      omatch = true;
    } else {
      std::cout << "Unrecognized:\t" << args[arg] << "\n";
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "Usage: csymmatch\n\t-pdbin-ref <filename>\n\t-pdbin <filename>\n\t-pdbout <filename>\n\t-connectivity-radius <radius/A>\n\t-origin-hand\nApply symmetry and cell shifts to each chain in 'pdbin' to obtain the best match to 'pdbin-ref'.\n";
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

  // do origin matching if required
  if ( omatch ) {
    // calculate structure factors
    clipper::Spacegroup rspgr = molref.spacegroup();
    clipper::Cell       rcell = molref.cell();
    clipper::Resolution rreso( 3.0 );
    clipper::HKL_sampling hkls( rcell, rreso );
    clipper::HKL_data<clipper::data32::F_phi> fphi1( rspgr, rcell, hkls );
    clipper::HKL_data<clipper::data32::F_phi> fphi2( rspgr, rcell, hkls );
    clipper::SFcalc_iso_fft<float> sfc;
    clipper::Atom_list atoms1 = molref.atom_list();
    clipper::Atom_list atoms2 = molwrk.atom_list();
    sfc( fphi1, atoms1 );
    sfc( fphi2, atoms2 );
    // get origin shift
    bool invert = false;
    clipper::Coord_frac x( 0.0, 0.0, 0.0 );
    clipper::OriginMatch<float>( invert, x, fphi1, fphi2 );
    if ( invert ) std::cout << std::endl << " Change of hand  : YES" << std::endl;
    else          std::cout << std::endl << " Change of hand  : NO" << std::endl;
    std::cout << " Change of origin: " << x.format() << std::endl << std::endl;
    // now move the atoms
    clipper::Mat33<> rot = clipper::Mat33<>::identity();
    clipper::Coord_orth trn = x.coord_orth( rcell );
    if ( invert )
      rot = clipper::Mat33<>(-1.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,-1.0);
    clipper::RTop_orth rt( rot, trn );
    molwrk.transform( rt );
  }

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

  // make a list of work atom groups
  std::vector<std::vector<clipper::MAtomIndex> > groups;
  for ( int c = 0; c < molwrk.size(); c++ ) {
    std::vector<clipper::MAtomIndex> group;
    for ( int r0 = 0; r0 < molwrk[c].size(); r0++ ) {
      // add this monomer
      for ( int a0 = 0; a0 < molwrk[c][r0].size(); a0++ )
	group.push_back( clipper::MAtomIndex(c,r0,a0) );
      // check if next monomer is connected to this
      bool iscon = false;
      int r1 = r0 + 1;
      if ( r1 < molwrk[c].size() ) {
	for ( int a0 = 0; a0 < molwrk[c][r0].size(); a0++ )
	  for ( int a1 = 0; a1 < molwrk[c][r1].size(); a1++ )
	    if ( ( molwrk[c][r0][a0].coord_orth() -
		   molwrk[c][r1][a1].coord_orth() ).lengthsq() < crad*crad )
	      iscon = true;
      }
      // if not, add this group
      if ( !iscon ) {
	groups.push_back( group );
	group.clear();
      }
    }
  }

  // now score each chain, symmetry and offset in turn
  for ( int g = 0; g < groups.size(); g++ ) {
    const std::vector<clipper::MAtomIndex>& group = groups[g];
    clipper::Atom_list atoms;
    for ( int i = 0; i < group.size(); i++ )
      atoms.push_back( molwrk[group[i].polymer()][group[i].monomer()]
		             [group[i].atom()] );
    double              bestscr = 0.0;
    int                 bestsym = 0;
    clipper::Coord_frac bestoff( 0.0, 0.0, 0.0 );
    const clipper::Coord_frac cfh( 0.5, 0.5, 0.5 );
    for ( int sym = 0; sym < spgr.num_symops(); sym++ ) {
      clipper::Atom_list atomw = atoms;
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
	    clipper::Coord_frac off( rint( off0.u() ) + du,
				     rint( off0.v() ) + dv,
				     rint( off0.w() ) + dw );
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
    for ( int i = 0; i < group.size(); i++ )
      molwrk[group[i].polymer()][group[i].monomer()]
            [group[i].atom()].transform( rtop );

    // user output
    clipper::MAtomIndex g0 = group.front();
    clipper::MAtomIndex g1 = group.back();
    std::cout << "Chain: " << molwrk[g0.polymer()].id()
	      << " " << molwrk[g0.polymer()][g0.monomer()].id()
	      << "-" << molwrk[g1.polymer()][g1.monomer()].id()
	      << " will be transformed as follows:" << std::endl
	      << "  Symmetry operator: "
	      << spgr.symop(bestsym).format() << std::endl
	      << "  Lattice shift:     "
	      << bestoff.format() << std::endl
	      << "  with normalised score:  "
	      << bestscr/double(atoms.size()) << std::endl;
  }

  // write file
  mmdbwrk.export_minimol( molwrk );
  mmdbwrk.write_file( oppdb );
}
