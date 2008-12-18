/*! \file parrot-lib.h parrot library */
/* Copyright 2008 Kevin Cowtan & University of York all rights reserved */

#include "simplex-lib.h"
#include "parrot-ncs.h"


#ifndef PARROT_NCSAVER
#define PARROT_NCSAVER


class NCSaver {
 public:
  static void local_correlation( clipper::NXmap<float>& correl, const clipper::NXmap<float>& r0, const clipper::NXmap<float>& r1, const double& local_radius );
  static clipper::Map_stats ncs_stats( const clipper::Xmap<float>& xmap, const double& local_radius );
  void ncs_mask_from_correl( clipper::NXmap<float>& mask, const clipper::NXmap<float>& correl, const double& level );
  void ncs_mask( clipper::NXmap<float>& mask, const clipper::Xmap<float>& xmap, const Local_rtop& nxop, const double& map_radius, const double& local_radius, const double& level, const int& nscl );
  void ncs_refine( Local_rtop& nxop, const clipper::Xmap<float>& xmap, const clipper::NXmap<float>& msk, bool refine );
  void ncs_average( clipper::Xmap<float>& xncs, clipper::Xmap<float>& xwgt, const clipper::Xmap<float>& xmap, const clipper::NXmap<float>& src, const Local_rtop& nxop ) const;

  std::vector<Local_rtop> filter_ncs_candidates( std::vector<Local_rtop> rtops, const clipper::Xmap<float>& xmap, const double& local_radius, const double& level, const double& minvol, const int& nscl );

  double correlation_old() const { return correl0; }
  double correlation_new() const { return correl1; }
  double mask_volume_ratio() const { return mskvol/totvol; }
  double mask_volume_total() const { return totvol; }
  const std::vector<double>& mask_volumes() const { return mskvols; }
 private:
  double correl0, correl1, mskvol, totvol;
  std::vector<double> mskvols;

  /*! Target function for nxmap rotation optimisation. */
  class Target_fn_xmap_mask_rtop : public Target_fn_order_zero {
  public:
    enum TYPE { CORREL, RMSD };
    Target_fn_xmap_mask_rtop() {}
    Target_fn_xmap_mask_rtop( const clipper::Xmap<float>& xmap, const clipper::NXmap<float>& src, const clipper::NXmap<float>& msk, const double& rot_step, const double& trn_step, const int& step );
    ~Target_fn_xmap_mask_rtop() {}
    int num_params() const { return 6; }
    //! \internal evaluate target function for EulerXYZr+uvw offset from rot_
    double operator() ( const Local_rtop& rot ) const;
    //! \internal evaluate target function for EulerXYZr+uvw offset from rot_
    double operator() ( const std::vector<double>& args ) const;
    //! \internal convert params to rotation
    Local_rtop local_rtop( const std::vector<double>& args ) const;
    //! refine rotation
    Local_rtop refine( const Local_rtop& rot );
  private:
    const clipper::Xmap<float>* xmap_;
    const clipper::NXmap<float>* src_;
    const clipper::NXmap<float>* msk_;
    double rot_step_, trn_step_;
    int step_;
    Local_rtop rot_;
  };

};


#endif
