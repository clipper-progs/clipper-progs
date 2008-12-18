/*! \file parrot-lib.h parrot library */
/* Copyright 2008 Kevin Cowtan & University of York all rights reserved */

#include "parrot-ncs.h"


#ifndef PARROT_NCSFIND
#define PARROT_NCSFIND


class NCSfind {
 public:
  NCSfind() { debug_ = false; }
  NCSfind( const double& tol_dst, const double& tol_ang ) : tol_dst_(tol_dst), tol_ang_(tol_ang) { debug_ = false; }
  void debug() { debug_ = true; }

  // NCS from partial model
  static clipper::String chain_sequence( const clipper::MPolymer& mp );
  static clipper::RTop_orth superpose( const clipper::MPolymer& mp1, const clipper::MPolymer& mp2, const double& rmsd, const int& nmin );
  static Local_rtop local_rtop( const clipper::MPolymer& mp1, const clipper::MPolymer& mp2, const double& rmsd, const int& nmin );

  // NCS from heavy atoms
  std::vector<Local_rtop> find_ncs_candidates( const clipper::Atom_list& atoms, const std::vector<clipper::String> atominfo, const clipper::Spacegroup& spgr, const clipper::Cell& cell ) const;
  std::vector<Local_rtop> find_ncs_candidates( const clipper::Atom_list& atoms, const clipper::Spacegroup& spgr, const clipper::Cell& cell ) const;

  // INTERNAL
 private:
  class CoordDescr {
  public:
    CoordDescr() {}
    CoordDescr( clipper::Coord_orth co, int index ) : co_(co), i_(index) {}
    const clipper::Coord_orth& coord_orth() const { return co_; }
    const int index() const                       { return i_; }
  private:
    clipper::Coord_orth co_; int i_;
  };

  static std::vector<NCSfind::CoordDescr> find_environ( std::vector<CoordDescr> coords, const clipper::Spacegroup& spgr, const clipper::Cell& cell, const clipper::Coord_orth& near, const int& nnear );
  static clipper::RTop_orth match_rtop( const std::vector<CoordDescr>& near1, const std::vector<CoordDescr>& near2, std::vector<std::pair<int,int> > match );
  static std::vector<double> match_diff( const std::vector<CoordDescr>& near1, const std::vector<CoordDescr>& near2, std::vector<std::pair<int,int> > match, clipper::RTop_orth rtop );
  static Local_rtop local_rtop( const std::vector<std::pair<CoordDescr,CoordDescr> >& coords );

  std::vector<std::vector<std::pair<int,int> > > match_atoms( const std::vector<CoordDescr>& near1, const std::vector<CoordDescr>& near2 ) const;

  double tol_ang_, tol_dst_;
  bool debug_;
};


#endif
