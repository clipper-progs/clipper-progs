// Clipper composite omit calculation
/* Copyright 2011 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-minimol.h>



class OmitCoordinates {
 public:
  OmitCoordinates() {}
  OmitCoordinates( clipper::Spacegroup spgr, clipper::Cell cell, int nomit );
  const std::vector<clipper::Coord_orth>& coordinates() { return coords; }
  double radius() const { return rad; }
 private:
  std::vector<clipper::Coord_orth> coords;
  double rad;
};


class MapFilterFn_smooth : public clipper::MapFilterFn_base {
public:
  MapFilterFn_smooth( const clipper::ftype& r0, const clipper::ftype& r1 ) : r0_(r0), r1_(r1) {}
    clipper::ftype operator() ( const clipper::ftype& radius ) const 
    {
      if ( radius > r1_ ) return 0.0;
      if ( radius < r0_ ) return 1.0;
      clipper::ftype x = ( radius - r0_ ) / ( r1_ - r0_ );
      return (2.0*x*x*x-3.0*x*x+1.0);
    }
private:
    clipper::ftype r0_, r1_;
};


OmitCoordinates::OmitCoordinates( clipper::Spacegroup spgr, clipper::Cell cell, int nomit )
{
  typedef clipper::Xmap<char>::Map_reference_index MRI;

  // calculate sphere radius
  double vasu = cell.volume() / spgr.num_symops();
  double r0 = 1.10*pow( vasu/double(nomit), 0.333 );

  // make a map
  clipper::Resolution reso( 1.0 );
  clipper::Grid_sampling grid( spgr, cell, reso );
  clipper::Xmap<float> xmap( spgr, cell, grid );
  xmap = 0.0;

  clipper::Grid_range gm( cell, grid, r0 );
  clipper::Xmap_base::Map_reference_coord i0, iu, iv, iw ;

  // now start packing spheres in ASU
  double cut2 = r0*r0;
  for ( MRI ix = xmap.first(); !ix.last(); ix.next() ) {
    if ( xmap[ix] == 0.0 ) {
      clipper::Coord_orth x0 = ix.coord_orth();
      coords.push_back( x0 );
      clipper::Coord_grid cent = x0.coord_frac( cell ).coord_grid( grid );
      clipper::Coord_grid g0 = cent + gm.min();
      clipper::Coord_grid g1 = cent + gm.max();
      i0 = clipper::Xmap_base::Map_reference_coord( xmap, g0 );
      for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() )
	for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
	  for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() ) {
	    double r2 = ( iw.coord_orth() - x0 ).lengthsq();
	    if ( r2 < cut2 ) {
	      float r = 1.0-r2/cut2;
	      xmap[iw] = clipper::Util::max( xmap[iw], r );
	    }
	  }
    }
  }

  // find the lowest points remaining
  float fcut = 0.40;
  for ( int i = 0; i < 2*nomit; i++ ) {
    MRI iz = xmap.first();
    for ( MRI iy = xmap.first(); !iy.last(); iy.next() )
      if ( xmap[iy] < xmap[iz] ) iz = iy;
    if ( xmap[iz] < fcut ) {
      clipper::Coord_orth x0 = iz.coord_orth();
      coords.push_back( x0 );
      clipper::Coord_grid cent = x0.coord_frac( cell ).coord_grid( grid );
      clipper::Coord_grid g0 = cent + gm.min();
      clipper::Coord_grid g1 = cent + gm.max();
      i0 = clipper::Xmap_base::Map_reference_coord( xmap, g0 );
      for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() )
	for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
	  for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() ) {
	    double r2 = ( iw.coord_orth() - x0 ).lengthsq();
	    if ( r2 < cut2 ) {
	      float r = 1.0-r2/cut2;
	      xmap[iw] = clipper::Util::max( xmap[iw], r );
	    }
	  }
    } else {
      break;
    }
  }


  // check the map for unfilled points
  int n0, n1, n2;
  n0 = n1 = n2 = 0;
  for ( MRI ix = xmap.first(); !ix.last(); ix.next() ) {
    if      ( xmap[ix] < 0.01 ) n0++;
    else if ( xmap[ix] < fcut ) n1++;
    else                        n2++;
  }
  //std::cout << coords.size() << "\n";
  //std::cout << n0 << "\t" << n1 << "\t" << n2 << "\n";
  //for ( int i = 0; i < coords.size(); i++ ) std::cout << i << " " << coords[i].format() << "\n";

  rad = r0*sqrt(1.0-fcut);
}



int main( int argc, char** argv )
{
  CCP4Program prog( "comit", "0.1.0", "$Date: 2011/12/14" );

  std::cout << "\nCopyright 2011 Kevin Cowtan and University of York. All rights reserved.\n\n";
  std::cout << " \n\n";

  // defaults
  clipper::String title;
  clipper::String ippdb = "NONE";
  clipper::String opmap = "NONE";
  clipper::String ipmtz = "NONE";
  clipper::String opmtz = "NONE";
  clipper::String ipcol_fo = "NONE";
  clipper::String ipcol_fc = "NONE";
  clipper::String opcol = "omit";
  clipper::String prefix = "NONE";
  int nomit = 30;
  double rpad = 3.0;
  int n_refln = 1000;
  int n_param = 20;
  int verbose = 0;

  // command input
  CCP4CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if        ( args[arg] == "-title" ) {
      if ( ++arg < args.size() ) title = args[arg];
    } else if ( args[arg] == "-pdbin" ) {
      if ( ++arg < args.size() ) ippdb = args[arg];
    } else if ( args[arg] == "-mapout" ) {
      if ( ++arg < args.size() ) opmap = args[arg];
    } else if ( args[arg] == "-mtzin" ) {
      if ( ++arg < args.size() ) ipmtz = args[arg];
    } else if ( args[arg] == "-mtzout" ) {
      if ( ++arg < args.size() ) opmtz = args[arg];
    } else if ( args[arg] == "-prefix" ) {
      if ( ++arg < args.size() ) prefix = args[arg];
    } else if ( args[arg] == "-colin-fo" ) {
      if ( ++arg < args.size() ) ipcol_fo = args[arg];
    } else if ( args[arg] == "-colin-fc" ) {
      if ( ++arg < args.size() ) ipcol_fc = args[arg];
    } else if ( args[arg] == "-colout" ) {
      if ( ++arg < args.size() ) opcol = args[arg];
    } else if ( args[arg] == "-nomit" ) {
      if ( ++arg < args.size() ) nomit = clipper::String(args[arg]).i();;
    } else if ( args[arg] == "-pad-radius" ) {
      if ( ++arg < args.size() ) rpad = clipper::String(args[arg]).f();;
    } else if ( args[arg] == "-verbose" ) {
      if ( ++arg < args.size() ) verbose = clipper::String(args[arg]).i();
    } else {
      std::cout << "\nUnrecognized:\t" << args[arg] << std::endl;
      args.clear();
    }
  }
  if ( args.size() <= 1 )
    clipper::Message::message( clipper::Message_fatal( "\nUsage: comit\n\t-pdbin <filename>\n\t-mtzin <filename>\n\t-mtzout <filename>\n\t-mapout <filename>\n\t-prefix <filename>\n\t-colin-fo <colpath>\n\t-colin-fc <colpath>\n\t-colout <colpath>\n\t-nomit <number of spheres>\n\t-pad-radius <radius/A>\n.\n" ) );

  const int mmdbflags = MMDBF_IgnoreBlankLines | MMDBF_IgnoreDuplSeqNum | MMDBF_IgnoreNonCoorPDBErrors | MMDBF_IgnoreRemarks;
  typedef clipper::HKL_info::HKL_reference_index HRI;
  typedef clipper::Xmap<float>::Map_reference_index MRI;

  /* THERE ARE TWO METHODS IMPLEMENTED:

     1. Work froma single refmac dataset and do the omit internally by
     blanking out density. This is fast, and gives reasonable results
     as long as you start from FC_ALL,PHIC_ALL, not FWT,PHWT.

     2. Run comit once to prepare a list of models for input to
     refmac, and then again to accumulate the results from the
     multiple refmac runs into an omit map. This is slow, but may be
     better.
  */

  if ( prefix == "NONE" ) {  // WORK FROM A SINGLE REFMAC RUN

    /*
      Work froma single refmac dataset and do the omit internally by
      blanking out density. This is fast, and gives reasonable results
      as long as you start from FC_ALL,PHIC_ALL, not FWT,PHWT.
    */

    // read MTZ
    clipper::HKL_data<clipper::data32::F_sigF> fobs;
    clipper::HKL_data<clipper::data32::F_phi>  fphi;
    clipper::CCP4MTZfile mtz;
    mtz.set_column_label_mode( clipper::CCP4MTZfile::Legacy );
    mtz.open_read( ipmtz );
    mtz.import_hkl_data( fobs, ipcol_fo );
    mtz.import_hkl_data( fphi, ipcol_fc );
    if ( opcol[0] != '/' ) opcol = mtz.assigned_paths()[0].notail()+"/"+opcol;
    mtz.close_read();

    // prepare omit coords
    clipper::Spacegroup spgr = fobs.spacegroup();
    clipper::Cell       cell = fobs.cell();
    OmitCoordinates oc( spgr, cell, nomit );
    const std::vector<clipper::Coord_orth>& omits = oc.coordinates();

    // calculate map
    clipper::Grid_sampling grid( spgr, cell, fobs.resolution() );
    clipper::Xmap<float> xmap( spgr, cell, grid );
    xmap.fft_from( fphi );

    // set up sigmaa calc
    clipper::HKL_data<clipper::data32::Flag> modeflag( fobs );
    for ( HRI ih = modeflag.first(); !ih.last(); ih.next() )
      if ( !fobs[ih].missing() )
	modeflag[ih].flag() = clipper::SFweight_spline<float>::BOTH;
      else
	modeflag[ih].flag() = clipper::SFweight_spline<float>::NONE;
    clipper::HKL_data<clipper::data32::F_phi> fb( fphi ), fd( fphi );
    clipper::HKL_data<clipper::data32::Phi_fom> phiw( fobs );
    clipper::HKL_data<clipper::data32::ABCD> abcd( fphi );

    // accumulate omit map
    clipper::Xmap<float> xrslt( spgr, cell, grid ), xwght( spgr, cell, grid );
    xrslt = xwght = 0.0;
    float r0 = oc.radius() - rpad;
    float r1 = oc.radius();
    float r2 = oc.radius() + rpad;
    MapFilterFn_smooth flt_sml( r0, r1 );
    MapFilterFn_smooth flt_lrg( r1, r2 );
    for ( int z = 0; z < omits.size(); z++ ) {
      clipper::Xmap<float> xomit = xmap;
      clipper::Coord_orth x0 = omits[z];
      // make the mask parameters
      clipper::Grid_range gm( cell, grid, r2 );
      clipper::Xmap_base::Map_reference_coord i0, iu, iv, iw ;
      clipper::Coord_grid cent = x0.coord_frac( cell ).coord_grid( grid );
      clipper::Coord_grid g0 = cent + gm.min();
      clipper::Coord_grid g1 = cent + gm.max();
      i0 = clipper::Xmap_base::Map_reference_coord( xmap, g0 );
      for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() )
	for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
	  for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() ) {
	    const double r2 = ( iw.coord_orth() - x0 ).lengthsq();
	    const float wgt = ( 1.0 - flt_lrg( sqrt( r2 ) ) );
	    xomit[iw] *= wgt;
	  }
      // calculate structure factors from omitted map
      clipper::HKL_data<clipper::data32::F_phi> fphic( fphi );
      xomit.fft_to( fphic );

      // sigmaa phase compination
      clipper::SFweight_spline<float> sfw( n_refln, n_param );
      sfw( fb, fd, phiw, fobs, fphic, modeflag );
      fphic = fb;

      // calc new map
      xomit.fft_from( fphic );

      // now accumulate 
      for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() )
	for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
	  for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() ) {
	    const double r2 = ( iw.coord_orth() - x0 ).lengthsq();
	    const float wgt = flt_sml( sqrt( r2 ) );
	    xrslt[iw] += wgt*xomit[iw];
	    xwght[iw] += wgt;
	  }
    }

    // accumulate the omit map
    for ( MRI ix = xmap.first(); !ix.last(); ix.next() )
      xrslt[ix] /= xwght[ix];

    // write map
    if ( opmap != "NONE" ) {
      clipper::CCP4MAPfile mapout;
      mapout.open_write( opmap );
      mapout.export_xmap( xrslt );
      mapout.close_write();
    }

    // write mtz
    if ( opmtz != "NONE" ) {
      // calculate
      xrslt.fft_to( fphi );
      clipper::SFweight_spline<float> sfw( n_refln, n_param );
      sfw( fb, fd, phiw, fobs, fphi, modeflag );
      abcd.compute( phiw, clipper::data32::Compute_abcd_from_phifom() );

      // write results
      mtz.open_append( ipmtz, opmtz );
      mtz.export_hkl_data( fphi, opcol );
      mtz.export_hkl_data( abcd, opcol );
      mtz.close_append();
    }

  } else {                   // PREPARE/ACCUMULATE REPEATED REFMAC RUNS

    /*
      Run comit once to prepare a list of models for input to refmac,
      and then again to accumulate the results from the multiple
      refmac runs into an omit map. This is slow, but may be better.
    */

    if ( ippdb != "NONE" ) {  // FIRST RUN - CREATE OMIT MODELS

      // prepare list of omit coordinates
      clipper::MMDBfile mmdbwrk;
      clipper::MiniMol molwrk;
      mmdbwrk.SetFlag( mmdbflags );
      mmdbwrk.read_file( ippdb );
      mmdbwrk.import_minimol( molwrk );

      // prepare omit coords
      clipper::Spacegroup spgr = molwrk.spacegroup();
      clipper::Cell       cell = molwrk.cell();
      OmitCoordinates oc( spgr, cell, nomit );
      const std::vector<clipper::Coord_orth>& omits = oc.coordinates();

      // write mask models
      clipper::MiniMol molmsk( spgr, cell );
      clipper::MPolymer mp;
      for ( int i = 0; i < omits.size(); i++ ) {
	clipper::MMonomer mm;
	mm.set_type( "OMI" );
	mm.set_seqnum( i+1 );
	clipper::MAtom ma = clipper::MAtom::null();
	ma.set_element( "O" ); ma.set_id( "OMIT" ); ma.set_occupancy( 0.0 );
	ma.set_u_iso( clipper::Util::b2u( oc.radius() ) );
	ma.set_coord_orth( omits[i] );
	mm.insert( ma );
	mp.insert( mm );
      }
      molmsk.insert( mp );
      clipper::MMDBfile mmdbomit;
      mmdbomit.export_minimol( molmsk );
      mmdbomit.write_file( prefix + "_OMIT" + ".pdb" );

      // prepare omit models
      double r0 = oc.radius() + rpad;
      for ( int i = 0; i < omits.size(); i++ ) {
	// prepare an omit model
	const clipper::Coord_frac cf0 = omits[i].coord_frac( cell );
	clipper::MiniMol molomit = molwrk;
	for ( int c = 0; c < molomit.size(); c++ )
	  for ( int r = 0; r < molomit[c].size(); r++ )
	    for ( int a = 0; a < molomit[c][r].size(); a++ ) {
	      const clipper::Coord_frac cf1 = molomit[c][r][a].coord_orth()
		.coord_frac( cell )
		.symmetry_copy_near( spgr, cell, cf0 );
	      const double r2 = (cf1-cf0).lengthsq(cell);
	      if ( r2 < r0*r0 ) molomit[c][r][a].set_occupancy( 0.0 );
	    }
	// write model
	clipper::String n = "0000" + clipper::String(i+1,1);
	clipper::MMDBfile mmdbomit;
	mmdbomit.export_minimol( molomit );
	mmdbomit.write_file( prefix + "_" + n.substr(n.length()-4,4) + ".pdb" );
      }

    } else {  // SECOND RUN - COMBINE OMIT MAPS

      // prepare list of omit coordinates
      clipper::MiniMol molmsk;
      clipper::MMDBfile mmdbwrk;
      clipper::MiniMol molwrk;
      mmdbwrk.SetFlag( mmdbflags );
      mmdbwrk.read_file( prefix + "_OMIT" + ".pdb" );
      mmdbwrk.import_minimol( molmsk );

      // get map parameters
      clipper::CCP4MTZfile mtzin;
      mtzin.open_read( prefix + "_0001.mtz" );
      clipper::Spacegroup spgr = mtzin.spacegroup();
      clipper::Cell       cell = mtzin.cell();
      clipper::Resolution reso = mtzin.resolution();
      clipper::Grid_sampling grid( spgr, cell, reso );
      mtzin.close_read();

      // read mtzs and make maps
      clipper::Xmap<float> xrslt( spgr, cell, grid );  // define map
      clipper::Xmap<float> xwght( spgr, cell, grid );  // define map
      xrslt = xwght = 0.0;
      for ( int i = 0; i < molmsk[0].size(); i++ ) {
	const clipper::Coord_orth co = molmsk[0][i][0].coord_orth();
	const double radius = clipper::Util::u2b( molmsk[0][i][0].u_iso() );

	// get mtz
	clipper::String n = "0000" + clipper::String(i+1,1);
	std::cout << prefix + "_" + n.substr(n.length()-4,4) + ".mtz" << "  "
		  << co.format() << "  " << radius << std::endl;
	clipper::CCP4MTZfile mtzin;
	mtzin.open_read( prefix + "_" + n.substr(n.length()-4,4) + ".mtz" );
	clipper::HKL_data<clipper::data32::F_phi> fphi;
	mtzin.import_hkl_data( fphi, "/*/*/[FWT,PHWT]" );
	mtzin.close_read();

	// calculate map
	clipper::Xmap<float> xmap( spgr, cell, grid );
	xmap.fft_from( fphi );

	// now accumulate
	const double r0 = radius - 0.0*rpad;
	const double r1 = radius + rpad;
	MapFilterFn_smooth flt( r0, r1 );
	clipper::Grid_range gm( cell, grid, r1 );
	clipper::Xmap_base::Map_reference_coord i0, iu, iv, iw ;
	clipper::Coord_grid cent = co.coord_frac( cell ).coord_grid( grid );
	clipper::Coord_grid g0 = cent + gm.min();
	clipper::Coord_grid g1 = cent + gm.max();
	i0 = clipper::Xmap_base::Map_reference_coord( xmap, g0 );
	for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() )
	  for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
	    for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() ) {
	      double r2 = ( iw.coord_orth() - co ).lengthsq();
	      float wgt = flt( sqrt( r2 ) );
	      xrslt[iw] += wgt*xmap[iw];
	      xwght[iw] += wgt;
	    }
      }

      // accumulate the omit map
      for ( MRI ix = xwght.first(); !ix.last(); ix.next() )
	xrslt[ix] /= xwght[ix];

      // write map
      if ( opmap != "NONE" ) {
	clipper::CCP4MAPfile mapout;
	mapout.open_write( opmap );
	mapout.export_xmap( xrslt );
	mapout.close_write();
      }

      // write structure factors
      if ( ipmtz != "NONE" && opmtz != "NONE" ) {
	clipper::HKL_data<clipper::data32::F_sigF>  fobs;
	clipper::CCP4MTZfile mtz;
	mtz.set_column_label_mode( clipper::CCP4MTZfile::Legacy );
	mtz.open_read( ipmtz );
	mtz.import_hkl_data( fobs, ipcol_fo );
	if ( opcol[0] != '/' ) opcol = mtz.assigned_paths()[0].notail()+"/"+opcol;
	mtz.close_read();

	// set up sigmaa calc
	clipper::HKL_data<clipper::data32::Flag> modeflag( fobs );
	for ( HRI ih = modeflag.first(); !ih.last(); ih.next() )
	  if ( !fobs[ih].missing() )
	    modeflag[ih].flag() = clipper::SFweight_spline<float>::BOTH;
	  else
	    modeflag[ih].flag() = clipper::SFweight_spline<float>::NONE;

	// calculate
	clipper::HKL_data<clipper::data32::F_phi>   fphi( fobs );
	xrslt.fft_to( fphi );
	clipper::HKL_data<clipper::data32::F_phi> fb( fphi ), fd( fphi );
	clipper::HKL_data<clipper::data32::Phi_fom> phiw( fobs );
	clipper::HKL_data<clipper::data32::ABCD> abcd( fphi );
	clipper::SFweight_spline<float> sfw( n_refln, n_param );
	sfw( fb, fd, phiw, fobs, fphi, modeflag );
	abcd.compute( phiw, clipper::data32::Compute_abcd_from_phifom() );

	// write results
	mtz.open_append( ipmtz, opmtz );
	mtz.export_hkl_data( fphi, opcol );
	mtz.export_hkl_data( abcd, opcol );
	mtz.close_append();
      }

    }

  }

}
