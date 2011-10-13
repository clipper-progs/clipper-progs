/*! \file parrot-ncsaver.cpp parrot library */
/* Copyright 2008 Kevin Cowtan & University of York all rights reserved */

#include "parrot-ncsaver.h"

#include <algorithm>

#ifdef rad2 // defined on Windows
# undef rad2
#endif


/*
void print_stats( const clipper::NXmap<float>& map )
{
  typedef clipper::NXmap<float>::Map_reference_index MRI;
  clipper::Range<double> range;
  for ( MRI ix = map.first(); !ix.last(); ix.next() ) range.include(map[ix]);
  clipper::Histogram hist( range, 50 );
  for ( MRI ix = map.first(); !ix.last(); ix.next() ) hist.accumulate(map[ix]);
  double gridvol = map.operator_grid_orth().rot().det();
  for ( int i = 0; i < hist.size(); i++ ) std::cout << hist.x(i) << " " << gridvol*hist.y(i) << std::endl;
}
*/


void NCSaver::local_correlation( clipper::NXmap<float>& correl, const clipper::NXmap<float>& r0, const clipper::NXmap<float>& r1, const double& local_radius )
{
  typedef clipper::NXmap<float>::Map_reference_index MRI;
  clipper::NXmap<float>
    x0( r0.grid(), r0.operator_orth_grid() ),
    x1( r0.grid(), r0.operator_orth_grid() ),
    x00( r0.grid(), r0.operator_orth_grid() ),
    x11( r0.grid(), r0.operator_orth_grid() ),
    x01( r0.grid(), r0.operator_orth_grid() );
  for ( MRI ix = r0.first(); !ix.last(); ix.next() ) {
    x0[ix] = r0[ix] * r0[ix];
    x1[ix] = r1[ix] * r1[ix];
    x00[ix] = r0[ix] * r1[ix];
  }
  clipper::MapFilterFn_step step( local_radius );
  clipper::MapFilter_fft<float>
    filter( step, 1.0, clipper::MapFilter_fft<float>::Relative );
  filter( x01, x00 );    // cross term
  filter( x00, x0 );     // second moment
  filter( x11, x1 );     // second moment
  filter( x0, r0 );      // first moment
  filter( x1, r1 );      // first moment
  for ( MRI ix = x01.first(); !ix.last(); ix.next() )
    x01[ix] = ( x01[ix] - x0[ix]*x1[ix] ) /
      sqrt( ( x00[ix]-x0[ix]*x0[ix] ) * ( x11[ix]-x1[ix]*x1[ix] ) );
  correl = x01;
}


clipper::Map_stats NCSaver::ncs_stats( const clipper::Xmap<float>& xmap, const double& local_radius )
{
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  clipper::Xmap<float>
    r0( xmap.spacegroup(), xmap.cell(), xmap.grid_sampling() ),
    r1( xmap.spacegroup(), xmap.cell(), xmap.grid_sampling() ),
    x0( xmap.spacegroup(), xmap.cell(), xmap.grid_sampling() ),
    x1( xmap.spacegroup(), xmap.cell(), xmap.grid_sampling() ),
    x00( xmap.spacegroup(), xmap.cell(), xmap.grid_sampling() ),
    x11( xmap.spacegroup(), xmap.cell(), xmap.grid_sampling() ),
    x01( xmap.spacegroup(), xmap.cell(), xmap.grid_sampling() );
  for ( MRI ix = xmap.first(); !ix.last(); ix.next() ) {
    r0[ix] = xmap[ix];
    r1[ix] = xmap.get_data( -ix.coord() );
    x0[ix] = r0[ix] * r0[ix];
    x1[ix] = r1[ix] * r1[ix];
    x00[ix] = r0[ix] * r1[ix];
  }
  clipper::MapFilterFn_step step( local_radius );
  clipper::MapFilter_fft<float>
    filter( step, 1.0, clipper::MapFilter_fft<float>::Relative );
  filter( x01, x00 );    // cross term
  filter( x00, x0 );     // second moment
  filter( x11, x1 );     // second moment
  filter( x0, r0 );      // first moment
  filter( x1, r1 );      // first moment
  for ( MRI ix = xmap.first(); !ix.last(); ix.next() )
    x01[ix] = ( x01[ix] - x0[ix]*x1[ix] ) /
      sqrt( ( x00[ix]-x0[ix]*x0[ix] ) * ( x11[ix]-x1[ix]*x1[ix] ) );
  clipper::Map_stats stats( x01 );
  return stats;
}


void NCSaver::ncs_mask_from_correl( clipper::NXmap<float>& mask, const clipper::NXmap<float>& correl, const double& level )
{
  typedef clipper::NXmap<float>::Map_reference_index MRI;
  clipper::NXmap<int> flag;
  clipper::Grid grid = correl.grid();
  mask.init( grid, correl.operator_orth_grid() );
  flag.init( grid, correl.operator_orth_grid() );

  // set map of flags
  int fval = 0;
  for ( MRI inx = correl.first(); !inx.last(); inx.next() ) {
    if ( correl[inx] > level ) flag[inx] = fval++;
    else                       flag[inx] = -1;
  }

  // now flag regions
  clipper::Coord_grid c0, c1;
  clipper::Coord_grid cn[6];
  for ( int i = 0; i < 6; i++ ) {
    cn[i][0] = cn[i][1] = cn[i][2] = 0;
    cn[i][i%3] = (i<3) ? -1 : 1;
  }
  bool done = false;
  while ( !done ) {
    done = true;
    // forward through map
    clipper::Coord_grid c0;
    for ( c0.w() = 0; c0.w() < grid.nw(); c0.w()++ )
      for ( c0.v() = 0; c0.v() < grid.nv(); c0.v()++ )
	for ( c0.u() = 0; c0.u() < grid.nu(); c0.u()++ ) {
	  int f0 = flag.get_data( c0 );
	  if ( f0 >= 0 ) {
	    for ( int i = 0; i < 6; i++ ) {
	      c1 = c0 + cn[i];
	      if ( grid.in_grid( c1 ) ) {
		int f1 = flag.get_data( c1 );
		if ( f1 >= 0 && f1 < f0 ) f0 = f1;
	      }
	    }
	    if ( f0 != flag.get_data( c0 ) ) {
	      flag.set_data( c0, f0 );
	      done = false;
	    }
	  }
	}
    // backward through map
    for ( c0.w() = grid.nw()-1; c0.w() >= 0; c0.w()-- )
      for ( c0.v() = grid.nv()-1; c0.v() >= 0; c0.v()-- )
	for ( c0.u() = grid.nu()-1; c0.u() >= 0; c0.u()-- ) {
	  int f0 = flag.get_data( c0 );
	  if ( f0 >= 0 ) {
	    for ( int i = 0; i < 6; i++ ) {
	      c1 = c0 + cn[i];
	      if ( grid.in_grid( c1 ) ) {
		int f1 = flag.get_data( c1 );
		if ( f1 >= 0 && f1 < f0 ) f0 = f1;
	      }
	    }
	    if ( f0 != flag.get_data( c0 ) ) {
	      flag.set_data( c0, f0 );
	      done = false;
	    }
	  }
	}
  }

  // extract the peak list
  std::vector<std::pair<int,int> > regions;
  for ( MRI inx = flag.first(); !inx.last(); inx.next() ) {
    int f = flag[inx];
    if ( f >= 0 ) {
      int r;
      for ( r = 0; r < regions.size(); r++ )
	if ( f == regions[r].second ) break;
      if ( r < regions.size() )
	regions[r].first++;
      else
	regions.push_back( std::pair<int,int>( 1, f ) );
    }
  }

  // mask volume stats
  double gridvol = mask.operator_grid_orth().rot().det();

  // make the mask
  if ( regions.size() > 0 ) {
    // sort
    std::sort( regions.begin(), regions.end() );
    std::reverse( regions.begin(), regions.end() );

    // fill mask
    int region0 = regions[0].second;
    for ( MRI inx = flag.first(); !inx.last(); inx.next() ) {
      if ( flag[inx] == region0 ) {
	mask[inx] = tanh((correl[inx]-level)/level);
      } else {
	mask[inx] = 0.0;
      }
    }

    // get stats
    mskvols = std::vector<double>( regions.size(), 0.0 );
    for ( int i = 0; i < regions.size(); i++ )
      mskvols[i] = gridvol * double( regions[i].first );
    mskvol = mskvols[0];
    totvol = 0.0;
    for ( int i = 0; i < regions.size(); i++ ) totvol += mskvols[i];
  } else {
    mskvols = std::vector<double>( 1, 0.0 );
    mskvol = 0.0;
    totvol = 1.0;
  }

  //std::cout << "NCS level: " << level << std::endl;
  //std::cout << "Correlations: " << std::endl;
  //print_stats( correl );
  //std::cout << "Mask: " << std::endl;
  //print_stats( mask );
}


void NCSaver::ncs_mask( clipper::NXmap<float>& mask, const clipper::Xmap<float>& xmap, const Local_rtop& nxop, const double& map_radius, const double& local_radius, const double& level, const int& nscl )
{
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  typedef clipper::NXmap<float>::Map_reference_index NRI;
  typedef float                 TYPE;
  typedef clipper::Xmap<TYPE>   XMAP;
  typedef clipper::Interp_linear INTERP;
  clipper::Coord_grid grid0;
  const clipper::Cell cell = xmap.cell();

  mskmul = 0;
  mskasu = mskovr = correls = 0.0;

  // calc coarse grid
  clipper::Grid_sampling gfine = xmap.grid_sampling();
  clipper::Grid_sampling
    gcoarse( gfine.nu()/nscl, gfine.nv()/nscl, gfine.nw()/nscl );
  // get grid offset for source coordinate
  grid0 = nxop.src().coord_frac( cell ).coord_grid( gcoarse );
  // calc grid containing the desired volume
  clipper::Grid_range gr0( cell, gcoarse, map_radius );
  // and offset by the base coordinate
  clipper::Grid_range gr1( gr0.min() + grid0, gr0.max() + grid0 );

  // construct new operators
  clipper::RTop_orth rtop = nxop.rtop_orth();

  // init 2 nxmaps, one each for unrotated and rotated density
  clipper::NXmap<float> nxmap0( cell, gcoarse, gr1 );
  clipper::NXmap<float> nxmap1( cell, gcoarse, gr1 );

  // populate the unrotated and rotated nxmap
  clipper::Coord_frac cf;
  for ( NRI inx = nxmap0.first(); !inx.last(); inx.next() ) {
    const clipper::Coord_orth co = inx.coord_orth();
    nxmap0[inx] = xmap.interp<INTERP>( xmap.coord_map( co ) );
    nxmap1[inx] = xmap.interp<INTERP>( xmap.coord_map( rtop * co ) );
  }

  // calc correlation
  clipper::NXmap<float> correl, msktmp;
  local_correlation( correl, nxmap0, nxmap1, local_radius );

  // calc mask
  ncs_mask_from_correl( msktmp, correl, level );

  // Calculate final mask extent
  clipper::Range<double> xu, xv, xw;
  for ( NRI inx = msktmp.first(); !inx.last(); inx.next() )
    if ( msktmp[inx] > 0.0 ) {
      const clipper::Coord_map cm = xmap.coord_map( inx.coord_orth() );
      xu.include( cm.u() );
      xv.include( cm.v() );
      xw.include( cm.w() );
    }

  if ( xu.range() <= 0.0 || xv.range() <= 0.0 || xw.range() <= 0.0 ) return;

  // Build final mask on fine grid
  const double d = 1.0;
  const clipper::Coord_map cm0( xu.min()-d, xv.min()-d, xw.min()-d );
  const clipper::Coord_map cm1( xu.max()+d, xv.max()+d, xw.max()+d );
  const clipper::Grid_range gcrop( cm0.coord_grid(), cm1.coord_grid() );
  mask.init( cell, gfine, gcrop );
  for ( NRI inx = mask.first(); !inx.last(); inx.next() ) {
    clipper::Coord_map cm = msktmp.coord_map( inx.coord_orth() );
    if ( INTERP::can_interp( msktmp, cm ) )
      mask[inx] = msktmp.interp<INTERP>( cm );
    else
      mask[inx] = 0.0;
  }

  // Mask normalisation (required for strictly correct weighting)
  // calculate mask self-overlap in the xmap ASU
  clipper::Xmap<int> nmap( xmap.spacegroup(), xmap.cell(), gfine );
  nmap = 0;
  clipper::Xmap<float>::Map_reference_coord ix( nmap );
  clipper::Coord_grid offset = nmap.coord_map( mask.coord_orth( clipper::Coord_map(0.0,0.0,0.0) ) ).coord_grid();
  for ( NRI inx = mask.first(); !inx.last(); inx.next() )
    if ( mask[inx] > 0.0 ) {
      ix.set_coord( inx.coord() + offset );
      nmap[ix] += 1;
    }
  // calculate median self-overlap
  int nmax = 0;
  for ( MRI iy = nmap.first(); !iy.last(); iy.next() )
    if ( nmap[iy] > nmax ) nmax = nmap[iy];
  std::vector<int> hist( nmax+1, 0 );
  for ( MRI iy = nmap.first(); !iy.last(); iy.next() )
    hist[nmap[iy]] += nmap[iy];
  int nmask(0), count(0), median(1);
  for ( int i = 1; i < hist.size(); i++ ) nmask += hist[i];
  for ( median = 1; median < hist.size(); median++ ) {
    count += hist[median];
    if ( count > nmask/2 ) break;
  }
  // weight mask
  double wgt = 1.0 / double(median);
  for ( NRI inx = mask.first(); !inx.last(); inx.next() ) mask[inx] *= wgt;
  // store stats
  mskmul = median;
  double asuvol = cell.volume() / double( nmap.spacegroup().num_symops() );
  mskasu = mskvols[0] / ( mskmul * asuvol );
  mskovr = double( hist[mskmul] ) / double( nmask );

  // Mask statistics
  // now get stats for 10A sphere about highest correl
  clipper::Coord_orth comax(0.0,0.0,0.0);
  double cmax = -1.0;
  for ( NRI inx = nxmap0.first(); !inx.last(); inx.next() )
    if ( correl[inx] > cmax ) { cmax = correl[inx]; comax = inx.coord_orth(); }
  double rad2 = clipper::Util::sqr(10.0);
  double sn(0.0), s0(0.0), s1(0.0), s00(0.0), s11(0.0), s01(0.0);
  for ( NRI inx = nxmap0.first(); !inx.last(); inx.next() ) {
    const clipper::Coord_orth co = inx.coord_orth();
    if ( (co-comax).lengthsq() < rad2 ) {
      const double rho0 = double( nxmap0[inx] );
      const double rho1 = double( nxmap1[inx] );
      sn  += 1.0;
      s0  += rho0;
      s1  += rho1;
      s00 += rho0 * rho0;
      s11 += rho1 * rho1;
      s01 += rho0 * rho1;
    }
  }
  if ( sn > 0.0 ) correls = (sn*s01-s0*s1)/sqrt((sn*s00-s0*s0)*(sn*s11-s1*s1));

  //for ( int i = 1; i < hist.size(); i++ ) std::cout << "overlap: " << i << " " << hist[i] << " " << ((i==median)?'*':' ') << std::endl;
  //std::cout << "Vol/ASU " << double( nmask * xmap.spacegroup().num_symops() ) / double( gfine.size() ) << std::endl;
}


void NCSaver::ncs_refine( Local_rtop& nxop, const clipper::Xmap<float>& xmap, const clipper::NXmap<float>& msk )
{
  if ( msk.is_null() ) return;
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  typedef clipper::NXmap<float>::Map_reference_index NRI;
  clipper::Cell cell = xmap.cell();
  clipper::Coord_grid offset = xmap.coord_map( msk.coord_orth( clipper::Coord_map(0.0,0.0,0.0) ) ).coord_grid();
  clipper::NXmap<float> rho0( msk.grid(), msk.operator_orth_grid() );
  clipper::Coord_orth src( 0.0, 0.0, 0.0 ), tgt( 0.0, 0.0, 0.0 );
  double sw(0.0);
  rho0 = 0.0;
  // get source map
  MRI ix( xmap );
  for ( NRI inx = msk.first(); !inx.last(); inx.next() ) {
    const float w = msk[inx];
    if ( w > 0.0 ) {
      ix.set_coord( inx.coord() + offset );
      rho0[inx] = xmap[ix];
      src += w * inx.coord_orth();
      sw += w;
    }
  }
  // get updated local rtop
  src = ( 1.0/sw ) * src;
  tgt = nxop.rtop_orth() * src;
  Local_rtop nxop0( nxop.rot(), src, tgt );
  // make target fn class
  NCSaver::Target_fn_xmap_mask_rtop tfn( xmap, rho0, msk, 0.05, 0.5, 5 );
  nxop = tfn.refine( nxop0 );
  correl0 = -tfn( nxop0 );
  correl1 = -tfn( nxop  );
}


void NCSaver::ncs_average( clipper::Xmap<float>& xncs, clipper::Xmap<float>& xwgt, const clipper::Xmap<float>& xmap, const clipper::NXmap<float>& msk, const Local_rtop& nxop ) const
{
  if ( msk.is_null() ) return;
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  typedef clipper::NXmap<float>::Map_reference_index NRI;
  typedef clipper::Interp_cubic INTERP;
  clipper::RTop_orth rtop = nxop.rtop_orth();
  clipper::Cell cell = xmap.cell();
  clipper::Coord_grid offset = xncs.coord_map( msk.coord_orth( clipper::Coord_map(0.0,0.0,0.0) ) ).coord_grid();
  MRI ix( xncs );
  for ( NRI inx = msk.first(); !inx.last(); inx.next() ) {
    const float w = msk[inx];
    if ( w > 0.0 ) {
      ix.set_coord( inx.coord() + offset );
      const clipper::Coord_orth co = rtop * inx.coord_orth();
      const clipper::Coord_frac cf = co.coord_frac( cell );
      const float r = xmap.interp<INTERP>( cf );
      xncs[ix] += w * r;
      xwgt[ix] += w;
    }
  }
}


std::vector<Local_rtop> NCSaver::filter_ncs_candidates( std::vector<Local_rtop> rtops, const clipper::Xmap<float>& xmap, const double& local_radius, const double& level, const double& minvol, const int& nscl )
{
  clipper::NXmap<float> mask;
  double asuvol = xmap.cell().volume()/double(xmap.spacegroup().num_symops());
  const double map_radius = pow( asuvol, 0.333 );
  std::vector<Local_rtop> rtopsf;
  for ( int r = 0; r < rtops.size(); r++ ) {
    ncs_mask( mask, xmap, rtops[r],
	      map_radius, local_radius, level, nscl );
    const double mskvol = mask_volumes()[0]/asuvol;
    if ( mskvol > minvol ) rtopsf.push_back( rtops[r] );
  }
  return rtopsf;
}


NCSaver::Target_fn_xmap_mask_rtop::Target_fn_xmap_mask_rtop( const clipper::Xmap<float>& xmap, const clipper::NXmap<float>& src, const clipper::NXmap<float>& msk, const double& rot_step, const double& trn_step, const int& step )
{
  xmap_ = &xmap;
  src_ = &src;
  msk_ = &msk;
  rot_step_ = rot_step;
  trn_step_ = trn_step;
  step_ = step;
}

double NCSaver::Target_fn_xmap_mask_rtop::operator() ( const Local_rtop& rot ) const
{
  typedef clipper::NXmap<float>::Map_reference_index NRI;
  const clipper::Xmap<float>& xmap = *xmap_;
  const clipper::NXmap<float>& src = *src_;
  const clipper::NXmap<float>& msk = *msk_;
  typedef clipper::Interp_cubic INTERP;
  clipper::RTop_orth rtop = rot.rtop_orth();
  clipper::Cell cell = xmap.cell();
  double sn(0.0), s0(0.0), s1(0.0), s00(0.0), s11(0.0), s01(0.0);
  for ( NRI inx = msk.first(); !inx.last(); inx.next() ) {
    clipper::Coord_grid cg = inx.coord();
    const float w = msk[inx];
    if ( cg.u()%step_ == 0 && cg.v()%step_ == 0 && cg.w()%step_ == 0 &&
	 w > 0.0 ) {
      const clipper::Coord_orth co = rtop * inx.coord_orth();
      const clipper::Coord_frac cf = co.coord_frac( cell );
      const float rho0 = src[inx];
      const float rho1 = xmap.interp<INTERP>( cf );
      sn  += double( w );
      s0  += double( w * rho0 );
      s1  += double( w * rho1 );
      s00 += double( w * rho0 * rho0 );
      s11 += double( w * rho1 * rho1 );
      s01 += double( w * rho0 * rho1 );
    }
  }
  return -(sn*s01-s0*s1) / sqrt((sn*s00-s0*s0)*(sn*s11-s1*s1));
}

double NCSaver::Target_fn_xmap_mask_rtop::operator() ( const std::vector<double>& args ) const
{ return (*this)( local_rtop( args ) ); }

Local_rtop NCSaver::Target_fn_xmap_mask_rtop::local_rtop( const std::vector<double>& args ) const
{
  return Local_rtop( clipper::Rotation(clipper::Euler<clipper::Rotation::EulerXYZs>(args[0],args[1],args[2])) * rot_.rot(), rot_.src(), clipper::Coord_orth(args[3],args[4],args[5]) + rot_.tgt() );
}

Local_rtop NCSaver::Target_fn_xmap_mask_rtop::refine( const Local_rtop& rot )
{
  // store initial rotation
  rot_ = rot;

  // calculate initial params
  std::vector<std::vector<double> > args;
  std::vector<double> arg(6,0.0);
  // identity
  clipper::Euler<clipper::Rotation::EulerXYZs> euler( 0.0, 0.0, 0.0 );
  arg[0] = euler.alpha();
  arg[1] = euler.beta();
  arg[2] = euler.gamma();
  args.push_back( arg );
  // rotation steps
  double step = 0.5 * rot_step_;
  euler = clipper::Euler<clipper::Rotation::EulerXYZs>( step, 0.0, 0.0 );
  arg[0] = euler.alpha();
  arg[1] = euler.beta();
  arg[2] = euler.gamma();
  args.push_back( arg );
  euler = clipper::Euler<clipper::Rotation::EulerXYZs>( 0.0, step, 0.0 );
  arg[0] = euler.alpha();
  arg[1] = euler.beta();
  arg[2] = euler.gamma();
  args.push_back( arg );
  euler = clipper::Euler<clipper::Rotation::EulerXYZs>( 0.0, 0.0, step );
  arg[0] = euler.alpha();
  arg[1] = euler.beta();
  arg[2] = euler.gamma();
  args.push_back( arg );
  // translation steps
  step = 0.5 * trn_step_;
  arg = args[0];
  arg[3] = step;
  arg[4] = 0.0;
  arg[5] = 0.0;
  args.push_back( arg );
  arg[3] = 0.0;
  arg[4] = step;
  arg[5] = 0.0;
  args.push_back( arg );
  arg[3] = 0.0;
  arg[4] = 0.0;
  arg[5] = step;
  args.push_back( arg );
  // simple refinement
  Optimiser_simplex os( 0.01, 50 );
  Local_rtop op = local_rtop( os( *this, args ) );
  return Local_rtop( op.rot().norm(), op.src(), op.tgt() );
}
