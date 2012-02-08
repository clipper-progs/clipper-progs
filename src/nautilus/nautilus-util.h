/*! \file nautilus-util.h nautilus library */
/* (C) 2011 Kevin Cowtan & University of York all rights reserved */


#ifndef NAUTILUS_UTIL_H
#define NAUTILUS_UTIL_H

#include <clipper/clipper-minimol.h>


class NautilusLog {
 public:
  NautilusLog() : currentcpu(0.0) { log(""); }
  void log( const clipper::String& id );
  void log( const clipper::String& id, const clipper::MiniMol& mol, bool view );
  void profile();
 private:
  std::vector<std::pair<std::string,double> > prof;
  double currentcpu;
};


#endif
