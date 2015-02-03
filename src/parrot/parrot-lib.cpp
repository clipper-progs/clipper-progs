/*! \file parrot-lib.cpp parrot library */
/* Copyright 2008 Kevin Cowtan & University of York all rights reserved */

#include "parrot-lib.h"

#include <clipper/clipper-ccp4.h>

#include <fstream>
extern "C" {
#include <stdlib.h>
#include <stdio.h>
}


void ParrotUtil::set_reference( clipper::String& mtz, clipper::String& pdb )
{
  const char* clibdptr = getenv( "CLIBD" );
  if ( clibdptr != NULL ) {
    clipper::String clibd( clibdptr );
    clipper::String path;
    std::ifstream file;
    if ( pdb == "NONE" ) {
      path = clibd+"/reference_structures/reference-1tqw.pdb";
      file.open( path.c_str(), std::ifstream::in ); file.close();
      if ( !file.fail() ) pdb = path;
    }
    if ( pdb == "NONE" ) {
      path = clibd+"\\reference_structures\\reference-1tqw.pdb";
      file.open( path.c_str(), std::ifstream::in ); file.close();
      if ( !file.fail() ) pdb = path;
    }
    if ( mtz == "NONE" ) {
      path = clibd+"/reference_structures/reference-1tqw.mtz";
      file.open( path.c_str(), std::ifstream::in ); file.close();
      if ( !file.fail() ) mtz = path;
    }
    if ( mtz == "NONE" ) {
      path = clibd+"\\reference_structures\\reference-1tqw.mtz";
      file.open( path.c_str(), std::ifstream::in ); file.close();
      if ( !file.fail() ) mtz = path;
    }
    if ( pdb == "NONE" || mtz == "NONE" ) 
      clipper::Message::message( clipper::Message_fatal( "No reference data specified and not in $CLIBD" ) );
  } else {
    clipper::Message::message( clipper::Message_fatal( "No reference data specified and $CLIBD not found" ) );
  }
}


double ParrotUtil::effective_resolution( const clipper::HKL_data<clipper::data32::Phi_fom>& phiw )
{
  typedef clipper::HKL_data_base::HKL_reference_index HRI;
  double smax = 1.0e-6;
  double sw = 0.0;
  double sn = 0.0;
  for ( HRI ih = phiw.first(); !ih.last(); ih.next() ) {
    sn += 1.0;
    if ( !phiw[ih].missing() ) sw += phiw[ih].fom();
    if ( ih.invresolsq() > smax ) smax = ih.invresolsq();
  }
  return ( 1.0 / ( sqrt( smax ) * pow( sw/sn, 0.333 ) ) );
}


void ParrotUtil::read_model( clipper::MiniMol& mol, clipper::String file, bool verbose )
{
  const int mmdbflags = ( ::mmdb::MMDBF_IgnoreBlankLines |
                          ::mmdb::MMDBF_IgnoreDuplSeqNum |
                          ::mmdb::MMDBF_IgnoreNonCoorPDBErrors |
                          ::mmdb::MMDBF_IgnoreRemarks );
  clipper::MMDBfile mmdb;
  mmdb.SetFlag( mmdbflags );
  if ( file != "NONE" ) {
    try {
      mmdb.read_file( file );
      mmdb.import_minimol( mol );
      clipper::Atom_list atoms = mol.atom_list();
      std::cout << "PDB file: " << file << std::endl;
      std::cout << "  Number of atoms read: " << atoms.size() << std::endl;
      if ( verbose ) for ( int i = 0; i < atoms.size(); i += atoms.size()-1 ) printf("%i6  %4s  %8.3f %8.3f %8.3f\n", i, atoms[i].element().c_str(), atoms[i].coord_orth().x(), atoms[i].coord_orth().y(), atoms[i].coord_orth().z() );
    } catch ( clipper::Message_fatal ) {
      std::cout << "FAILED TO READ PDB FILE: " << file << std::endl;
    }
  }
}


void ParrotUtil::output_ncs_mask( clipper::String prefix, const clipper::NXmap<float>& msk, const clipper::Cell& cell, int cyc, int op )
{
  clipper::String nums( "0123456789" );
  int c = cyc+1; int o = op+1;
  clipper::String name = prefix+"_cycle"+nums[c/100%10]+nums[c/10%10]+nums[c%10]+"_operator"+nums[o/100%10]+nums[o/10%10]+nums[o%10]+".map";
  clipper::CCP4MAPfile file;
  file.open_write( name );
  file.set_cell( cell );
  file.export_nxmap( msk );
  file.close_write();
}


void ParrotUtil::mask_expand( clipper::Xmap<float>& mskmod, const clipper::Xmap<float>& msk, bool max )
{
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  float rm, r[6];
  for ( MRI ix = msk.first(); !ix.last(); ix.next() ) {
    r[0] = msk.get_data( ix.index_offset( -1,  0,  0 ) );
    r[1] = msk.get_data( ix.index_offset(  1,  0,  0 ) );
    r[2] = msk.get_data( ix.index_offset(  0, -1,  0 ) );
    r[3] = msk.get_data( ix.index_offset(  0,  1,  0 ) );
    r[4] = msk.get_data( ix.index_offset(  0,  0, -1 ) );
    r[5] = msk.get_data( ix.index_offset(  0,  0,  1 ) );
    rm = r[0];
    if ( max ) {
      for ( int i = 1; i < 6; i++ ) if ( r[i] > rm ) rm = r[i];
    } else {
      for ( int i = 1; i < 6; i++ ) if ( r[i] < rm ) rm = r[i];
    }
    mskmod[ix] = rm;
  }
}


void ParrotUtil::mask_smooth( clipper::Xmap<float>& msk, int ngrw, int nshr )
{
  clipper::Xmap<float> msktmp = msk;
  for ( int i = 0; i < ngrw; i++ ) {
    mask_expand( msktmp, msk, true );
    mask_expand( msk, msktmp, true );
  }
  for ( int i = 0; i < nshr; i++ ) {
    mask_expand( msktmp, msk, false );
    mask_expand( msk, msktmp, false );
  }
}


ParrotUtil::Map_stats ParrotUtil::masked_stats( const clipper::Xmap<float>& map, const clipper::Xmap<float>& msk )
{
  typedef clipper::Xmap_base::Map_reference_index MRI;
  double min, max, s0, s1, s2;
  s0 = s1 = s2 = 0.0;
  min =  1.0e20;
  max = -1.0e20;
  for ( MRI ix = map.first(); !ix.last(); ix.next() ) {
    const double w = msk[ix];
    const double r = map[ix];
    s0 += w;
    s1 += w * r;
    s2 += w * r * r;
    if ( r > max ) max = r;
    if ( r < min ) min = r;
  }
  s1 /= s0;
  s2 /= s0;
  s2 = s2 - s1*s1;
  s2 = sqrt( s2>0.0 ? s2 : 0.0 );
  return Map_stats( min, max, s1, s2 );
}


clipper::Generic_ordinal ParrotUtil::masked_ordinal( const clipper::Xmap<float>& map, const clipper::Xmap<float>& msk, const Map_stats& stats )
{
  typedef clipper::Xmap_base::Map_reference_index MRI;
  clipper::Generic_ordinal result( stats.range(), 200 );
  for ( MRI ix = map.first(); !ix.last(); ix.next() )
    result.accumulate( map[ix], msk[ix] );
  result.prep_ordinal();
  return result;
}


void ParrotUtil::solvent_mask( clipper::Xmap<float>& msk, const clipper::MiniMol& mol )
{
  clipper::EDcalc_mask<float> mskcalc;
  clipper::String sel = "ALA,ARG,ASN,ASP,CYS,GLN,GLU,GLY,HIS,ILE,LEU,LYS,MET,PHE,PRO,SER,THR,TRP,TYR,VAL,MSE";
  clipper::MPolymer mp, mp_wrk;
  for ( int c = 0; c < mol.size(); c++ )
    for ( int r = 0; r < mol[c].size(); r++ )
      if ( sel.find( mol[c][r].type() ) != std::string::npos )
        mp.insert( mol[c][r] );
  mskcalc( msk, mp.atom_list() );
  mask_smooth( msk, 1, 2 );
}


void ParrotUtil::density_modify( clipper::Xmap<float>& map_mod, const clipper::Xmap<float>& map_wrk, const clipper::Xmap<float>& map_ncs, const clipper::Xmap<float>& map_nwt, const clipper::Xmap<float>& msk_prt, const clipper::Xmap<float>& msk_sol, const clipper::Xmap<float>& map_ref, const clipper::Xmap<float>& map_sim, const clipper::Xmap<float>& msk_ref, bool dosolv, bool dohist )
{
  typedef clipper::Xmap_base::Map_reference_index MRI;

  // apply NCS
  clipper::Xmap<float> map_tmp( map_wrk.spacegroup(), map_wrk.cell(),
                                map_wrk.grid_sampling() );
  for ( MRI ix = map_wrk.first(); !ix.last(); ix.next() ) {
    map_tmp[ix] = ( map_wrk[ix] + map_ncs[ix] ) / ( map_nwt[ix] + 1.0 );
  }

  // make data arrays
  ParrotUtil::Map_stats stat_ref, stat_sim, stat_prt, stat_sol;
  stat_ref = ParrotUtil::masked_stats( map_ref, msk_ref );
  stat_sim = ParrotUtil::masked_stats( map_sim, msk_ref );
  stat_prt = ParrotUtil::masked_stats( map_tmp, msk_prt );
  stat_sol = ParrotUtil::masked_stats( map_tmp, msk_sol );

  // make ordinals
  clipper::Generic_ordinal ord_ref, ord_sim, ord_prt, ord_sol;
  ord_ref = ParrotUtil::masked_ordinal( map_ref, msk_ref, stat_ref );
  ord_sim = ParrotUtil::masked_ordinal( map_sim, msk_ref, stat_ref );
  ord_prt = ParrotUtil::masked_ordinal( map_tmp, msk_prt, stat_prt );
  ord_sol = ParrotUtil::masked_ordinal( map_tmp, msk_sol, stat_prt );
  clipper::Generic_ordinal ord_inv = ord_ref;
  ord_inv.invert();

  // now do histogram matching
  for ( MRI ix = map_tmp.first(); !ix.last(); ix.next() ) {
    const double r0 = ord_prt.ordinal( map_tmp[ix] );
    //double x = std::max( std::min( r0, 0.995 ), 0.005 );
    double x = 2.0 * r0 - 1.0;
    x = x / ( sqrt( 1.0 + 0.01 * x * x ) );
    x = ( x + 1.0 ) / 2.0;
    const double r1 = ord_inv.ordinal( x );
    const double r2 = ( r1 - stat_sim.mean() ) / stat_sim.std_dev();
    const double r3 = ( r2 * stat_prt.std_dev() ) + stat_prt.mean();
    map_mod[ix] = r3;
  }

  for ( MRI ix = map_mod.first(); !ix.last(); ix.next() ) {
    double wp = dohist ? msk_prt[ix] : 0.0;
    double ws = dosolv ? msk_sol[ix] : 0.0;
    double wo = ( 1.0 - wp ) * ( 1.0 - ws );
    double sw = wp + ws + wo;
    map_mod[ix] = (wp*map_mod[ix] + wo*map_tmp[ix] + ws*stat_sol.mean()) / sw;
  }
}


// cmatt = volume / mw*z*nsym
double ParrotUtil::matthews_probability( double cmatt, double resol, TYPE type )
{
  const double rbin[] = {1.199, 1.501, 1.650, 1.801, 1.901, 2.001, 2.201, 2.401, 2.601, 2.801, 3.101, 3.501, 5.001, 5.001, 5.001};
  const double p0[] = {0.085, 0.312, 0.400, 0.503, 0.597, 0.729, 1.052, 1.781, 2.852, 3.386, 3.841, 4.281, 4.592, 1.503, 0.257};
  const double vmbar[] = {2.052, 2.102, 2.122, 2.132, 2.140, 2.155, 2.171, 2.182, 2.191, 2.192, 2.205, 2.211, 2.210, 2.256, 2.324};
  const double wcoeff[] = {0.213, 0.214, 0.220, 0.225, 0.226, 0.231, 0.236, 0.241, 0.242, 0.244, 0.244, 0.244, 0.245, 0.446, 0.327};
  const double acoeff[] = {28.38, 102.7, 187.5, 339.3, 434.1, 540.5, 686.2, 767.9, 835.9, 856.9, 854.0, 849.6, 846.7, 136.6, 47.10};
  const double scoeff[] = {0.953, 0.807, 0.775, 0.702, 0.648, 0.640, 0.635, 0.589, 0.584, 0.542, 0.500, 0.485, 0.480, 1.180, 0.466};
  int n;
  if      ( type == PROTEIN ) {
    for ( n = 0; n < 12; n++ )
      if ( resol < rbin[n] ) break;
  } else if ( type == NUCLEOTIDE ) {
    n = 13;
  } else if ( type == MIXED ) {
    n = 14;
  } else {
    n = 14;
  }
  double z=(cmatt-vmbar[n])/wcoeff[n];
  double vm_prob = p0[n]+acoeff[n]*(exp(-exp(-z)-z*scoeff[n]+1.0));
  return vm_prob;
}


std::vector<std::pair<double,double> > ParrotUtil::solvent_probability( clipper::MMoleculeSequence seq, clipper::Spacegroup spgr, clipper::Cell cell, clipper::Resolution reso )
{
  // count sequence elements
  clipper::String codes = "ABCDEFGHIJKLMNOPQRSTUVWYXZ";
  double protmw[] = {  71.0,   0.0, 102.0, 114.0, 128.0,  /* ABCDE */
                      147.0,  57.0, 138.0, 113.0,   0.0,  /* FGHIJ */
                      129.0, 113.0, 131.0, 114.0,   0.0,  /* KLMNO */
                       97.0, 128.0, 157.0,  87.0, 101.0,  /* PQRST */
                        0.0,  99.0, 186.0, 114.0, 163.0,  /* UVWXY */
                        0.0 };                            /* Z */
  double nuclmw[] = { 313.0, -999., 289.0, -999., -999.,  /* ABCDE */
                      -999., 329.0, -999., -999., -999.,  /* FGHIJ */
                      -999., -999., -999., -999., -999.,  /* KLMNO */
                      -999., -999., -999., -999., 304.0,  /* PQRST */
                      290.0, -999., -999., -999., -999.,  /* UVWXY */
                      -999. };                            /* Z */
  double mwp(0.0), mwn(0.0);
  for ( int c = 0; c < seq.size(); c++ ) {
    clipper::String s = seq[c].sequence();
    double n(0.0), p(0.0);
    for ( int m = 0; m < s.size(); m++ ) {
      int t = codes.find( s.substr(m,1) );
      if ( t >= 0 && t < 26 ) {
        p += protmw[t];
        n += nuclmw[t];
      }
    }
    if ( n > p ) mwn += n;
    else         mwp += p;
  }
  double mw  = mwp + mwn;

  // calculate density
  const double prtden = 1.0/0.74;
  const double dnaden = 1.0/0.50;

  double density = ( mwp * prtden + mwn * dnaden ) / mw;

  double z = double( spgr.num_symops() );
  double volume = cell.volume();

  TYPE type = MIXED;
  if ( mwp > 0.90 * mw ) type = PROTEIN;
  if ( mwn > 0.80 * mw ) type = NUCLEOTIDE;
  double matt, prob, solc;
  std::vector<std::pair<double,double> > result;
  for ( int nncs = 0; nncs < 60; nncs++ ) {
    matt = 1.0e6;
    if ( nncs > 0 ) matt = volume / (mw*z*double(nncs));
    prob = matthews_probability( matt, reso.limit(), type );
    solc = 1.0 - 1.0/(0.602*matt*density);
    if ( solc <= 0.0 ) break;
    result.push_back( std::pair<double,double>( prob, solc ) );
  }
  double s = 0.0;
  for ( int n = 0; n < result.size(); n++ ) s += result[n].first;
  for ( int n = 0; n < result.size(); n++ ) result[n].first /= s;
  return result;
}


double ParrotUtil::random() 
{
  return 0.001*(::rand()%2000-1000);
}


// LOGGING

ParrotUtil::ParrotUtil( clipper::String title )
{
  title_ = title;
  cyc = 0;
  ncsdata.resize(1);
  rfldata.resize(1);
}

void ParrotUtil::log_ncs_operator( Local_rtop nxop )
{
      clipper::Polar_ccp4 polar = nxop.rot().polar_ccp4();
      clipper::Euler_ccp4 euler = nxop.rot().euler_ccp4();
      std::cout << " Polar rotation/deg: "
                << clipper::Util::rad2d(polar.omega()) << ","
                << clipper::Util::rad2d(polar.phi() ) << ","
                << clipper::Util::rad2d(polar.kappa()) << std::endl;
      std::cout << " Euler rotation/deg: "
                << clipper::Util::rad2d(euler.alpha()) << ","
                << clipper::Util::rad2d(euler.beta() ) << ","
                << clipper::Util::rad2d(euler.gamma()) << std::endl;
      std::cout << " Source: " << nxop.src().format() << std::endl;
      std::cout << " Target: " << nxop.tgt().format() << std::endl;
}

void ParrotUtil::log_cycle( int c )
{
  cyc = c;
  std::cout << std::endl << "-- Cycle: " << cyc
            << " --------------------------------" << std::endl << std::endl;
  ncsdata.resize(c+1);
  rfldata.resize(c+1);
}

void ParrotUtil::log_histogram_graph( const clipper::Xmap<float>& msk_ref, const clipper::Xmap<float>& msk_prt, const clipper::Xmap<float>& msk_sol, const clipper::Xmap<float>& map_ref, const clipper::Xmap<float>& map_sim, const clipper::Xmap<float>& map_wrk, const clipper::Xmap<float>& map_mod ) const
{
  // make stats
  ParrotUtil::Map_stats stat_ref, stat_sim, stat_prt, stat_sol;
  stat_ref = ParrotUtil::masked_stats( map_ref, msk_ref );
  stat_sim = ParrotUtil::masked_stats( map_sim, msk_ref );
  stat_prt = ParrotUtil::masked_stats( map_wrk, msk_prt );
  stat_sol = ParrotUtil::masked_stats( map_wrk, msk_sol );
  clipper::Generic_ordinal ord_ref, ord_sim, ord_prt, ord_sol, ord_mod_prt, ord_mod_sol;
  ord_ref = ParrotUtil::masked_ordinal( map_ref, msk_ref, stat_ref );
  ord_sim = ParrotUtil::masked_ordinal( map_sim, msk_ref, stat_ref );
  ord_prt = ParrotUtil::masked_ordinal( map_wrk, msk_prt, stat_prt );
  ord_sol = ParrotUtil::masked_ordinal( map_wrk, msk_sol, stat_prt );
  ord_mod_prt = ParrotUtil::masked_ordinal( map_mod, msk_prt, stat_prt );
  ord_mod_sol = ParrotUtil::masked_ordinal( map_mod, msk_sol, stat_prt );

  // print stats
  const int ntab = 20;
  const double rng = stat_prt.max()-stat_prt.min();
  const double scl = stat_sim.std_dev() / stat_prt.std_dev();
  printf("$TABLE :Cycle %i Electron density histograms:\n",cyc);
  printf("$GRAPHS :Protein:N:1,4,5,6::Solvent:N:1,7,8::Simulation:N:1,3,4: $$\n");
  printf("rho_min rho_max   Simulatn P_init P_trgt P_mod  S_init S_mod $$\n");
  printf("$$\n");
  for ( int i = 0; i < ntab; i++ ) {
    double r1 = stat_prt.min() + rng * double(i  ) / double(ntab);
    double r2 = stat_prt.min() + rng * double(i+1) / double(ntab);
    double rr1 = scl * ( r1 - stat_prt.mean() ) + stat_sim.mean();
    double rr2 = scl * ( r2 - stat_prt.mean() ) + stat_sim.mean();
    double rr = ord_ref.ordinal(rr2) - ord_ref.ordinal(rr1);
    double rs = ord_sim.ordinal(rr2) - ord_sim.ordinal(rr1);
    double wp = ord_prt.ordinal(r2) - ord_prt.ordinal(r1);
    double ws = ord_sol.ordinal(r2) - ord_sol.ordinal(r1);
    double mp = ord_mod_prt.ordinal(r2) - ord_mod_prt.ordinal(r1);
    double ms = ord_mod_sol.ordinal(r2) - ord_mod_sol.ordinal(r1);
    printf( "%7.3f %7.3f      %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f\n",
            r1, r2, rs, wp, rr, mp, ws, ms );
  }
  printf("$$\n");
  std::cout << std::endl << std::endl;
}


void ParrotUtil::log_sigmaa_graph( clipper::SFweight_spline<float>& sfw, const clipper::HKL_data<clipper::data32::Flag>& flagwt ) const
{
  printf( "Log likelihood:%14.6e      Log likelihood (free):%14.6e\n\n",
          sfw.log_likelihood_work(), sfw.log_likelihood_free() );
  printf("$TABLE :Cycle %i SigmaA statistics:\n",cyc);
  printf("$GRAPHS :SigmaA statistics:N:1,2,3: $$\n");
  printf(" 1/resol^2  sigmaA(s)  sigmaA(w) $$\n");
  printf("$$\n");
  int npweight = sfw.params_scale().size();
  clipper::Resolution_ordinal ord;
  ord.init( flagwt, flagwt.hkl_info().cell(), 1.0 ); ord.invert();
  for ( int i = 0; i < npweight; i++ ) {
    double err = sfw.params_error()[i];
    double s = ord.ordinal( (double(i)+0.5)/double(npweight) );
    double sigmaa1 = sfw.params_scale()[i];
    double sigmaa2 = ( err < 1.0 ) ? sqrt(1.0-err) : 0.0;
    printf( " %8.3f   %8.3f   %8.3f\n", s, sigmaa1, sigmaa2 );
  }
  printf("$$\n");
  std::cout << std::endl << std::endl;
}


void ParrotUtil::log_rfl_stats( clipper::HKL_data<clipper::data32::F_sigF>& wrk_f, clipper::HKL_data<clipper::data32::F_phi>& wrk_fp, clipper::HKL_data<clipper::data32::Phi_fom>& wrk_pw, const clipper::HKL_data<clipper::data32::Flag>& flagwt )
{
  typedef clipper::HKL_data_base::HKL_reference_index HRI;
  double sn(0.0),sfom(0.0);
  double snw(0.0),s1w(0.0),s2w(0.0),s11w(0.0),s22w(0.0),s12w(0.0);
  double snf(0.0),s1f(0.0),s2f(0.0),s11f(0.0),s22f(0.0),s12f(0.0);
  for ( HRI ih = wrk_f.first(); !ih.last(); ih.next() )
    if ( !wrk_f[ih].missing() && !wrk_fp[ih].missing() && !wrk_pw[ih].missing() ) {
      sn += 1.0;
      sfom += wrk_pw[ih].fom();
      if ( flagwt[ih].flag() == clipper::SFweight_spline<float>::BOTH ) {
        snw  += 1.0;
        s1w  +=  wrk_f[ih].f();
        s2w  += wrk_fp[ih].f();
        s11w +=  wrk_f[ih].f() *  wrk_f[ih].f();
        s22w += wrk_fp[ih].f() * wrk_fp[ih].f();
        s12w +=  wrk_f[ih].f() * wrk_fp[ih].f();
      } else {
        snf  += 1.0;
        s1f  +=  wrk_f[ih].f();
        s2f  += wrk_fp[ih].f();
        s11f +=  wrk_f[ih].f() *  wrk_f[ih].f();
        s22f += wrk_fp[ih].f() * wrk_fp[ih].f();
        s12f +=  wrk_f[ih].f() * wrk_fp[ih].f();
      }
    }
  double m(0.0), cw(0.0), cf(0.0);
  if ( sn > 0.0 )
    m = sfom / sn;
  if ( snw > 0.0 )
    cw = (snw*s12w-s1w*s2w) / sqrt((snw*s11w-s1w*s1w)*(snw*s22w-s2w*s2w));
  if ( snf > 0.0 )
    cf = (snf*s12f-s1f*s2f) / sqrt((snf*s11f-s1f*s1f)*(snf*s22f-s2f*s2f));
  rflcyinf data = { m, cw, cf };
  rfldata[cyc] = data;
}


void ParrotUtil::log_ncs_stats( Local_rtop nxop0, Local_rtop nxop1, double vol, int mult, double correl, double correl0, double correl1, double rcont, double rover )
{
  // store
  ncsopinf data = { vol, correl };
  ncsdata[cyc].push_back( data );

  // output
  printf( "NCS operator: %3i\n", int(ncsdata[cyc].size()) );
  printf( " NCS masking: Mask volume as fraction of ASU: %8.2f   Multiplicity: %i\n", vol, mult );
  printf( "              Contiguity score: %6.3f   Self-overlap score: %6.3f\n", rcont, rover );
  if ( correl0 > -1.0 && correl1 > -1.0 ) {
    clipper::Euler_ccp4 euler1 = nxop0.rot().euler_ccp4();
    clipper::Euler_ccp4 euler2 = nxop1.rot().euler_ccp4();
    printf( " NXop refinement- correlation before: %6.3f, after: %6.3f\n", correl0, correl1 );
    printf( " NXop old: %6.1f %6.1f %6.1f   %6.1f %6.1f %6.1f   %6.1f %6.1f %6.1f\n", clipper::Util::rad2d(euler1.alpha()), clipper::Util::rad2d(euler1.beta() ), clipper::Util::rad2d(euler1.gamma()), nxop0.src().x(), nxop0.src().y(), nxop0.src().z(), nxop0.tgt().x(), nxop0.tgt().y(), nxop0.tgt().z() );
    printf( " NXop new: %6.1f %6.1f %6.1f   %6.1f %6.1f %6.1f   %6.1f %6.1f %6.1f\n", clipper::Util::rad2d(euler2.alpha()), clipper::Util::rad2d(euler2.beta() ), clipper::Util::rad2d(euler2.gamma()), nxop1.src().x(), nxop1.src().y(), nxop1.src().z(), nxop1.tgt().x(), nxop1.tgt().y(), nxop1.tgt().z() );
  }
  std::cout << std::endl;
}


void ParrotUtil::log_ncs_table() const
{
  printf("NCS operator statistics:\n");
  printf(" Operator_number  Mask_volume/ASU  Correlation\n");
  for ( int n = 0; n < ncsdata[cyc].size(); n++ )
    printf( "         %4i         %8.3f     %8.3f\n",
            n+1, ncsdata[cyc][n].ncsvol, ncsdata[cyc][n].ncscor );
}


void ParrotUtil::log_summary_graphs() const
{
  // fom graph
  printf("$TABLE :Reflection statistics:\n");
  printf("$GRAPHS :Reflection statistics:N:1,2,3,4: $$\n");
  printf(" Cycle   FOM     Fcorrel(work)   Fcorrel(free)$$\n");
  printf("$$\n");
  for ( int c = 0; c < rfldata.size(); c++ ) {
    printf( " %4i    %6.3f  %12.3f    %12.3f\n", c,
            rfldata[c].meanfom, rfldata[c].fcorrw, rfldata[c].fcorrf );
  }
  printf("$$\n");
  std::cout << std::endl;

  // ncs graph
  int ncyc = ncsdata.size();
  int nncs = 0;
  std::vector<double> volmin(ncyc,0.0), volmax(ncyc,0.0), volmean(ncyc,0.0);
  std::vector<double> cormean(ncyc,0.0);
  for ( int c = 0; c < ncyc; c++ ) {
    if ( ncsdata[c].size() > 0 ) {
      if ( ncsdata[c].size() > nncs ) nncs = ncsdata[c].size();
      volmin[c] = volmax[c] = ncsdata[c][0].ncsvol;
      for ( int r = 0; r < ncsdata[c].size(); r++ ) {
        cormean[c] += ncsdata[c][r].ncscor;
        volmean[c] += ncsdata[c][r].ncsvol;
        volmin[c] = std::min( ncsdata[c][r].ncsvol, volmin[c] );
        volmax[c] = std::max( ncsdata[c][r].ncsvol, volmax[c] );
      }
      cormean[c] /= double( ncsdata[c].size() );
      volmean[c] /= double( ncsdata[c].size() );
    }
  }
  if ( nncs > 0 ) {
    printf("$TABLE :NCS statistics:\n");
    printf("$GRAPHS :NCS correlation statistics, best 10A sphere:N:1,2: $$\n");
    printf("$GRAPHS :NCS mask volumes:N:1,3,4,5: $$\n");
    printf(" Cycle   Correlation   MaskVol(mean) MaskVol(min)  MaskVol(max)$$\n");
    printf("$$\n");
    for ( int c = 1; c < ncyc; c++ ) {
      printf( " %4i    %12.3f  %12.3f  %12.3f  %12.3f\n", c,
              cormean[c], volmean[c], volmin[c], volmax[c] );
    }
    printf("$$\n");
  }
  std::cout << std::endl;
}

void ParrotUtil::xml( clipper::String xml ) const
{
  // ncs graph
  int ncyc = ncsdata.size();
  int nncs = 0;
  std::vector<double> volmin(ncyc,0.0), volmax(ncyc,0.0), volmean(ncyc,0.0);
  std::vector<double> cormean(ncyc,0.0);
  for ( int c = 0; c < ncyc; c++ ) {
    if ( ncsdata[c].size() > 0 ) {
      if ( ncsdata[c].size() > nncs ) nncs = ncsdata[c].size();
      volmin[c] = volmax[c] = ncsdata[c][0].ncsvol;
      for ( int r = 0; r < ncsdata[c].size(); r++ ) {
        cormean[c] += ncsdata[c][r].ncscor;
        volmean[c] += ncsdata[c][r].ncsvol;
        volmin[c] = std::min( ncsdata[c][r].ncsvol, volmin[c] );
        volmax[c] = std::max( ncsdata[c][r].ncsvol, volmax[c] );
      }
      cormean[c] /= double( ncsdata[c].size() );
      volmean[c] /= double( ncsdata[c].size() );
    }
  }

  // xml output
  clipper::String xmltmp = xml+".tmp";
  FILE *f = fopen( xmltmp.c_str(), "w" );
  if ( f == NULL ) clipper::Message::message( clipper::Message_fatal( "Error: Could not open xml temporary file: "+xmltmp ) );
  fprintf( f, "<ParrotResult>\n" );
  fprintf( f, " <Title>%s</Title>\n", title_.c_str() );
  // initial
  fprintf( f, " <Nncs>%i</Nncs>\n", nncs );
  // by cycle
  fprintf( f, " <Cycles>\n" );
  for ( int c = 0; c < rfldata.size(); c++ ) {
    fprintf( f, "  <Cycle>\n" );
    fprintf( f, "   <Number>%i</Number><MeanFOM>%6.3f</MeanFOM><Fcorrel>%6.3f</Fcorrel><FreeFcorrel>%6.3f</FreeFcorrel><NCScormean>%6.3f</NCScormean><NCSvolmean>%6.3f</NCSvolmean><NCSvolmin>%6.3f</NCSvolmin><NCSvolmax>%6.3f</NCSvolmax>\n", c, rfldata[c].meanfom, rfldata[c].fcorrw, rfldata[c].fcorrf, cormean[c], volmean[c], volmin[c], volmax[c] );
    fprintf( f, "  </Cycle>\n" );
  }
  fprintf( f, " </Cycles>\n" );
  // overall
  int c = rfldata.size()-1;
  fprintf( f, " <Final>\n" );
  fprintf( f, "  <MeanFOM>%6.3f</MeanFOM><Fcorrel>%6.3f</Fcorrel><FreeFcorrel>%6.3f</FreeFcorrel><NCScormean>%6.3f</NCScormean><NCSvolmean>%6.3f</NCSvolmean><NCSvolmin>%6.3f</NCSvolmin><NCSvolmax>%6.3f</NCSvolmax>\n", rfldata[c].meanfom, rfldata[c].fcorrw, rfldata[c].fcorrf, cormean[c], volmean[c], volmin[c], volmax[c] );
  fprintf( f, "  <Operators>\n" );
  for ( int r = 0; r < ncsdata[c].size(); r++ ) fprintf( f, "   <Operator>%i</Operator><NCScorrel>%6.3f</NCScorrel><NCSvolume>%6.3f</NCSvolume>\n", r+1, ncsdata[c][r].ncscor, ncsdata[c][r].ncsvol );
  fprintf( f, "  </Operators>\n" );
  fprintf( f, " </Final>\n" );
  fprintf( f, "</ParrotResult>\n" );
  fclose(f);
  rename( xmltmp.c_str(), xml.c_str() );
}
