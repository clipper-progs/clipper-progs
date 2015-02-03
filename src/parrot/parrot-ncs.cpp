/*! \file parrot-ncs.cpp parrot library */
/* Copyright 2008 Kevin Cowtan & University of York all rights reserved */

#include "parrot-ncs.h"

#include <algorithm>

extern "C" {
  #include <stdio.h>
}

/*
 Transform an operator by pre- and post- RToperators.
*/
Local_rtop Local_rtop::transform( const clipper::RTop_orth& r1, const clipper::RTop_orth& r2 ) const
{
  Local_rtop res;
  res.src() = r1 * src();
  res.tgt() = r2 * tgt();
  res.rot() = clipper::Rotation(r2.rot() * rot().matrix() * r1.rot().inverse());
  res.rot().norm();
  return res;
}


/*
 Calculate whether two operators match to within a given tolerance.
*/
std::pair<double,Local_rtop> Local_rtop::symm_match( const Local_rtop& other, const clipper::Spacegroup& spgr, const clipper::Cell& cell, const double& tol_dst, const double& tol_ang ) const
{
  // prepare orthogonal symops
  std::vector<clipper::RTop_orth> symop_orth( spgr.num_symops() );
  for ( int i = 0; i < symop_orth.size(); i++ )
    symop_orth[i] = spgr.symop(i).rtop_orth(cell);

  // find scored symmetry match between two ops
  Local_rtop rcbest( clipper::Rotation::null(), clipper::Coord_orth(), clipper::Coord_orth() );
  double scbest = 1.0e9;
  clipper::Coord_frac cf;
  for ( int sym1 = 0; sym1 < spgr.num_symops(); sym1++ )
    for ( int sym2 = 0; sym2 < spgr.num_symops(); sym2++ ) {
      Local_rtop rc1 = other;
      Local_rtop rc2 = transform( symop_orth[sym1], symop_orth[sym2] );
      rc2.src() = rc2.src().coord_frac(cell).lattice_copy_near(rc1.src().coord_frac(cell)).coord_orth(cell);
      rc2.tgt() = rc2.tgt().coord_frac(cell).lattice_copy_near(rc1.tgt().coord_frac(cell)).coord_orth(cell);
      clipper::Rotation rot = rc1.rot().inverse() * rc2.rot();
      double a = rot.abs_angle();
      double r = sqrt( ( rc2.rtop_orth()*rc1.src() - rc1.tgt() ).lengthsq() );
      double a2 = (a*a)/(tol_ang*tol_ang);
      double r2 = (r*r)/(tol_dst*tol_dst);
      double s2 = r2 + a2;
      if ( s2 < scbest ) {
	scbest = s2;
	rcbest = rc2;
      }
    }
  return std::pair<double,Local_rtop>( scbest, rcbest );
}


/*
 Remove duplicate operators
*/
std::vector<Local_rtop> Local_rtop::exclude_duplicate( const std::vector<Local_rtop>& ncsops, const clipper::Spacegroup& spgr, const clipper::Cell& cell, double tol_dst, double tol_ang )
{
  std::vector<Local_rtop> ncsopsi;
  for ( int i = 0; i < ncsops.size(); i++ ) {
    Local_rtop rtop = ncsops[i];
    int j;
    for ( j = 0; j < ncsopsi.size(); j++ )
      if ( rtop.symm_match( ncsopsi[j], spgr, cell, tol_dst, tol_ang ).first
	   < 2.0 ) break;
    if ( j == ncsopsi.size() ) ncsopsi.push_back( ncsops[i] );
  }
  return ncsopsi;
}


/*
 Remove identity operators
*/
std::vector<Local_rtop> Local_rtop::exclude_identity( const std::vector<Local_rtop>& ncsops, const clipper::Spacegroup& spgr, const clipper::Cell& cell, double tol_dst, double tol_ang )
{
  std::vector<Local_rtop> ncsopsi;
  Local_rtop rtid( clipper::Rotation::zero(),
		   clipper::Coord_orth(0.0,0.0,0.0),
		   clipper::Coord_orth(0.0,0.0,0.0) );
  for ( int i = 0; i < ncsops.size(); i++ ) {
    Local_rtop rtop = ncsops[i];
    if ( !( rtop.symm_match( rtid, spgr, cell, tol_dst, tol_ang ).first
	    < 2.0 ) ) ncsopsi.push_back( ncsops[i] );
  }
  return ncsopsi;
}


/*
 Calculate closest operator to proper form.
*/
Local_rtop Local_rtop::proper( const clipper::Spacegroup& spgr, const clipper::Cell& cell ) const
{
  // prepare orthogonal symops
  std::vector<clipper::RTop_orth> symop_orth( spgr.num_symops() );
  for ( int i = 0; i < symop_orth.size(); i++ )
    symop_orth[i] = spgr.symop(i).rtop_orth(cell);

  // most proper version
  int maxncs = 20;
  double wloop( 1.0 ), wscrw( 1.0 ), wtran( 0.001 );
  // make list of candidate operators
  std::vector<std::pair<double,Local_rtop> > resultsym;
  clipper::RTop_orth rtid = clipper::RTop_orth::identity();
  for ( int j = 0; j < spgr.num_symops(); j++ ) {
    Local_rtop rsym, rcel;
    rsym = transform( rtid, symop_orth[j] );
    clipper::Coord_frac cf = (rsym.tgt()-rsym.src()).coord_frac(cell);
    clipper::Coord_frac df( floor(cf.u()), floor(cf.v()), floor(cf.w()) );
    clipper::Coord_frac d;
    for ( int u = 0; u <= 1; u++ )
      for ( int v = 0; v <= 1; v++ )
	for ( int w = 0; w <= 1; w++ ) {
	  d = df + clipper::Coord_frac(double(u),double(v),double(w));
	  rcel.src() = rsym.src();
	  rcel.tgt() = rsym.tgt() - d.coord_orth(cell);
	  rcel.rot() = rsym.rot();
	  resultsym.push_back(std::pair<double,Local_rtop>(0.0,rcel));
	}
  }
  // score them for properness
  for ( int k = 0; k < resultsym.size(); k++ ) {
    clipper::Mat33<> mat = resultsym[k].second.rot().matrix();
    clipper::Vec3<> v0( mat(0,0) - 1.0, mat(1,0), mat(2,0) );
    clipper::Vec3<> v1( mat(0,1), mat(1,1) - 1.0, mat(2,1) );
    clipper::Vec3<> v2( mat(0,2), mat(1,2), mat(2,2) - 1.0 );
    clipper::Vec3<> v3 = clipper::Vec3<>::cross( v1, v2 );
    clipper::Vec3<> v4 = clipper::Vec3<>::cross( v0, v2 );
    clipper::Vec3<> v5 = clipper::Vec3<>::cross( v0, v1 );
    if ( v3*v3 > v4*v4 && v3*v3 > v5*v5 )
      v0 = v3.unit();
    else if ( v4*v4 > v5*v5 )
      v0 = v4.unit();
    else
      v0 = v5.unit();
    clipper::RTop_orth rtop = resultsym[k].second.rtop_orth();
    double scrwshift = clipper::Util::sqr( v0 * rtop.trn() );
    double transhift = ( resultsym[k].second.tgt() - resultsym[k].second.src() ).lengthsq();
    double loopshift = 1.0e9;
    clipper::Coord_orth init, next;
    init = next = resultsym[k].second.src();
    for ( int r = 0; r < maxncs; r++ ) {
      next = rtop * next;
      loopshift = std::min( loopshift, (next-init).lengthsq() );
    }
    resultsym[k].first = wloop*loopshift+wscrw*scrwshift+wtran*transhift;
  }
  // and pick the best
  std::sort( resultsym.begin(), resultsym.end() );
  return resultsym[0].second;
}


void Local_rtop::print_nxops( const clipper::String msg, const std::vector<Local_rtop>& nxops ) {
  std::cout << "NX operators: " << msg << std::endl;
  for ( int r = 0; r < nxops.size(); r++ )
    printf ("euler,src,tgt: %6.1f %6.1f %6.1f, %6.1f %6.1f %6.1f, %6.1f %6.1f %6.1f\n", clipper::Util::rad2d(nxops[r].rot().euler_ccp4().alpha()), clipper::Util::rad2d(nxops[r].rot().euler_ccp4().beta()), clipper::Util::rad2d(nxops[r].rot().euler_ccp4().gamma()), nxops[r].src().x(), nxops[r].src().y(), nxops[r].src().z(), nxops[r].tgt().x(), nxops[r].tgt().y(), nxops[r].src().z() );
  std::cout << std::endl;
}
