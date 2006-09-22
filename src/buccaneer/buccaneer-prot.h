/*! \file buccaneer-prot.h buccaneer library */
/* (C) 2002-2006 Kevin Cowtan & University of York all rights reserved */


#ifndef BUCCANEER_PROT
#define BUCCANEER_PROT

#include "buccaneer-lib.h"

#include <clipper/clipper-minimol.h>

#include <deque>


//! C-alpha group
/*! The Ca-group class represents a residue by the alpha Carbon and
  its neighbouring N and C main-chain atoms. It has methods to return
  the operator to generate the group from a standard orientation, and
  to generate the next or previous residue given two Ramachandran
  angles (one of this residue and one of the new one). */
class Ca_group {
 public:
  Ca_group() {}  //!< null constructor
  //! constructor: from atom coordinates
  Ca_group( const clipper::Coord_orth& n, const clipper::Coord_orth& ca, const clipper::Coord_orth& c ) : coord_n_(n), coord_ca_(ca), coord_c_(c) {}
  //! get N atom coordinate
  const clipper::Coord_orth& coord_n()  const { return coord_n_; }
  //! get C-alpha atom coordinate
  const clipper::Coord_orth& coord_ca() const { return coord_ca_; }
  //! get C atom coordinate
  const clipper::Coord_orth& coord_c()  const { return coord_c_; }
  //! get C-beta atom coordinate
  clipper::Coord_orth coord_cb() const
    { return rtop_from_std_ori()*std_coord_cb(); }
  //! get operator generating this group from standard orientation
  clipper::RTop_orth rtop_from_std_ori() const;
  //! get operator centering C-beta standard orientation
  clipper::RTop_orth rtop_beta_carbon() const;
  //! generate next Ca_group using Ramachandran angles
  Ca_group next_ca_group( const clipper::ftype& psi, const clipper::ftype& phi ) const;
  //! generate previous Ca_group using Ramachandran angles
  Ca_group prev_ca_group( const clipper::ftype& phi, const clipper::ftype& psi ) const;
  static clipper::Coord_orth std_coord_ca()
    { return clipper::Coord_orth( 0.00, 0.00, 0.00 ); }  //!< std C-a
  static clipper::Coord_orth std_coord_c()
    { return clipper::Coord_orth( 0.87, 0.00, 1.23 ); }  //!< std N
  static clipper::Coord_orth std_coord_n()
    { return clipper::Coord_orth( 0.83, 0.00, -1.18 ); }  //!< std C
  static clipper::Coord_orth std_coord_cb()
    { return clipper::Coord_orth( -1.03, -1.11, 0.00 ); }  //!< std C-b
 private:
  clipper::Coord_orth coord_n_, coord_ca_, coord_c_;
};


//! Chain of Ca-groups
/*! A Ca-chain is a std::deque (double ended queue) of Ca_group. In
  addition it has methods to return the Ramachandran angles of any
  residue. */
class Ca_chain : public std::deque<Ca_group> {
 public:
  //! Return Ramachandran phi for any residue except first in a chain
  clipper::ftype ramachandran_phi( const int& resno ) const;
  //! Return Ramachandran psi for any residue except last in a chain
  clipper::ftype ramachandran_psi( const int& resno ) const;
};


//! Planar residue group
/*! The Planar-residue-group class represents the planar atoms
  surrounding a peptide bond. It is represented by the alpha Carbon, C
  atom, and N atom of the next residue. It may also be described within a
  single residue by referring to the O atom instead of the next N.  It has
  methods to return the operator to generate the group from a standard
  orientation, and to generate the next or previous residue given two
  Ramachandran angles. */
class Pr_group {
 public:
  enum TYPE { CaCN, CaCO };  //!< atom types used in constructor
  Pr_group() {}  //!< null constructor
  //! constructor: from atom coordinates (Ca, C, N[+1] or Ca, C, O)
  Pr_group( const clipper::Coord_orth& ca, const clipper::Coord_orth& c, const clipper::Coord_orth& other, const TYPE& type );
  //! get C-alpha atom coordinate
  const clipper::Coord_orth& coord_ca() const { return coord_ca_; }
  //! get C atom coordinate
  const clipper::Coord_orth& coord_c()  const { return coord_c_; }
  //! get next N atom coordinate
  const clipper::Coord_orth& coord_n_next() const { return coord_n_; }
  //! generate O atom coordinate
  clipper::Coord_orth coord_o() const;
  //! generate next C-alpha atom coordinate
  clipper::Coord_orth coord_ca_next() const;
  //! get operator generating this group from standard orientation
  clipper::RTop_orth rtop_from_std_ori() const;
  //! generate next Pr_group using Ramachandran angles
  Pr_group next_pr_group( const clipper::ftype& phi, const clipper::ftype& psi ) const;
  //! generate previous Pr_group using Ramachandran angles
  Pr_group prev_pr_group( const clipper::ftype& psi, const clipper::ftype& phi ) const;
 private:
  clipper::Coord_orth coord_ca_, coord_c_, coord_n_;
};


//! Usefull tools for manipulating proteins
class ProteinTools {
 public:
  static int residue_index( clipper::String code, bool translate=true );
  static clipper::String residue_code_1( int index );
  static clipper::String residue_code_3( int index );
  static clipper::String residue_code( clipper::String code, bool translate=true );
  static bool chain_tidy( clipper::MiniMol& target, const clipper::MiniMol& source );
  static bool copy_residue_types( clipper::MiniMol& target, const clipper::MiniMol& source );
  static bool globularise( clipper::MiniMol& mol, const clipper::Coord_frac cent );
  static bool globularise( clipper::MiniMol& mol );
 private:
  static const int ntype;
  static const char rtype1[21];
  static const char rtype3[21][4];
};


#endif