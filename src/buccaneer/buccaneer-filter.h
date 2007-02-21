/*! \file buccaneer-filter.h buccaneer library */
/* (C) 2002-2006 Kevin Cowtan & University of York all rights reserved */


#include "buccaneer-prot.h"


//! Class for merging overlapped Ca chains and grouping by symmetry
class Ca_filter {
 public:
  Ca_filter( clipper::ftype sig_cut = 3.0 ) : sigcut(sig_cut) {}
  bool operator() ( clipper::MiniMol& mol2, const clipper::MiniMol& mol1, const
		    clipper::Xmap<float>& xmap );
  int num_filtered() const;

 private:
  double sigcut;
  int num_filter;
};
