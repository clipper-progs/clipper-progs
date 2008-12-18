/*! \file parrot-lib.h parrot library */
/* Copyright 2008 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-minimol.h>


#ifndef PARROT_LIB
#define PARROT_LIB


class ParrotUtil {
 public:
  enum TYPE { PROTEIN, NUCLEOTIDE, MIXED };

  class Map_stats {
  public:
    Map_stats() {}
    Map_stats( double min, double max, double mean, double std_dev ) :
      min_(min), max_(max), mean_(mean), std_dev_(std_dev) {}
    double min()     const { return min_;     }
    double max()     const { return max_;     }
    double mean()    const { return mean_;    }
    double std_dev() const { return std_dev_; }
    clipper::Range<> range() const { return clipper::Range<>( min_, max_ ); }
  private:
    double min_, max_, mean_, std_dev_;
  };

  // command input
  static void set_reference( clipper::String& mtz, clipper::String& pdb );

  // reflection data utilities
  static double effective_resolution( const clipper::HKL_data<clipper::data32::Phi_fom>& phiw );

  // coordinate utilities
  static void read_model( clipper::MiniMol& mol, clipper::String file );

  // mask utilities
  static void mask_expand( clipper::Xmap<float>& mskmod,
			   const clipper::Xmap<float>& msk, bool max );
  static void mask_smooth( clipper::Xmap<float>& msk, int ngrw=1, int nshr=1 );
  static Map_stats masked_stats( const clipper::Xmap<float>& map, const clipper::Xmap<float>& msk );
  static clipper::Generic_ordinal masked_ordinal( const clipper::Xmap<float>& map, const clipper::Xmap<float>& msk, const Map_stats& stats );
  static void solvent_mask( clipper::Xmap<float>& msk, const clipper::MiniMol& mol );

  // density modification utilities
  static void density_modify( clipper::Xmap<float>& map_mod, const clipper::Xmap<float>& map_wrk, const clipper::Xmap<float>& map_ncs, const clipper::Xmap<float>& map_nwt, const clipper::Xmap<float>& msk_prt, const clipper::Xmap<float>& msk_sol, const clipper::Xmap<float>& map_ref, const clipper::Xmap<float>& map_sim, const clipper::Xmap<float>& msk_ref, bool dosolv, bool dohist );

  // cell content utilities
  static double matthews_probability( double cmatt, double resol, TYPE type );
  static std::vector<std::pair<double,double> > solvent_probability( clipper::MMoleculeSequence seq, clipper::Spacegroup spgr, clipper::Cell cell, clipper::Resolution reso );

  // other utilites
  template<class T> static T sqr( const T& v ) { return v*v; }
  static double random();
};


#endif
