/*! \file parrot-lib.h parrot library */
/* Copyright 2008 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-minimol.h>


#ifndef PARROT_NCS
#define PARROT_NCS


/*! Class for a candidate operator. */
class Local_rtop
{
 public:
  Local_rtop() {}
  Local_rtop( const clipper::Rotation& rot, const clipper::Coord_orth& src, const clipper::Coord_orth& tgt ) : rot_(rot), src_(src), tgt_(tgt) {}
  const clipper::Rotation& rot() const { return rot_; }
  const clipper::Coord_orth& src() const { return src_; }
  const clipper::Coord_orth& tgt() const { return tgt_; }
  clipper::Rotation& rot() { return rot_; }
  clipper::Coord_orth& src() { return src_; }
  clipper::Coord_orth& tgt() { return tgt_; }
  bool is_null() const { return rot().is_null(); }
  // invert
  Local_rtop inverse() const
    { return Local_rtop( rot().inverse(), tgt(), src() ); }
  // convert to RT operator
  clipper::RTop_orth rtop_orth() const
    { return clipper::RTop_orth( rot_.matrix(), tgt_ - rot_.matrix()*src_ ); }

  // apply symmetry operators
  Local_rtop transform( const clipper::RTop_orth& r1, const clipper::RTop_orth& r2 ) const;
  // compare two Local_rtops
  std::pair<double,Local_rtop> symm_match( const Local_rtop& other, const clipper::Spacegroup& spgr, const clipper::Cell& cell, const double& tol_dst, const double& tol_ang ) const;
  // remove duplicates
  static std::vector<Local_rtop> exclude_duplicate( const std::vector<Local_rtop>& ncsops, const clipper::Spacegroup& spgr, const clipper::Cell& cell, double tol_dst = 3.0, double tol_ang = 0.1 );
  // remove identity
  static std::vector<Local_rtop> exclude_identity( const std::vector<Local_rtop>& ncsops, const clipper::Spacegroup& spgr, const clipper::Cell& cell, double tol_dst = 5.0, double tol_ang = 0.2 );
  // convert to closest proper form
  Local_rtop proper( const clipper::Spacegroup& spgr, const clipper::Cell& cell ) const;

  // operator for sorting
  bool operator < ( const Local_rtop& other ) const { return rot().w() < other.rot().w(); }
 private:
  clipper::Rotation rot_; clipper::Coord_orth src_, tgt_;
};


#endif
