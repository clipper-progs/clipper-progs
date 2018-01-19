/*! \file nautilus-util.h nautilus library */
/* (C) 2011 Kevin Cowtan & University of York all rights reserved */


#ifndef NAUTILUS_UTIL_H
#define NAUTILUS_UTIL_H

#include <clipper/clipper-minimol.h>


class NautilusUtil {
 public:
  static void set_reference( clipper::String& pdb );
};


class NautilusLog {
 public:
  NautilusLog( clipper::String& title ) : title_(title), currentcpu(0.0) { log(""); }
  void log( const clipper::String& id );
  void log( const clipper::String& id, const clipper::MiniMol& mol, bool view );
  clipper::String log_info( const clipper::MiniMol& mol, bool summary );
  void xml( const clipper::String& file ) const; //, const clipper::MiniMol& mol ); edited SWH
  void profile();
 private:
  struct cycdat { int nchns, nseq, nres, nmax; }; // added by SWH
  std::vector<cycdat> data; //added by SWH
  std::vector<std::pair<std::string,double> > prof;
  clipper::String title_;
  double currentcpu;
};


#endif
