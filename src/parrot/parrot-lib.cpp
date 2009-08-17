/*! \file parrot-lib.cpp parrot library */
/* Copyright 2008 Kevin Cowtan & University of York all rights reserved */

#include "parrot-lib.h"

#include <fstream>
extern "C" {
#include <stdlib.h>
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
  const int mmdbflags = ( MMDBF_IgnoreBlankLines |
			  MMDBF_IgnoreDuplSeqNum |
			  MMDBF_IgnoreNonCoorPDBErrors |
			  MMDBF_IgnoreRemarks );
  clipper::MMDBfile mmdb;
  mmdb.SetFlag( mmdbflags );
  if ( file != "NONE" ) {
    try {
      mmdb.read_file( file );
      mmdb.import_minimol( mol );
      std::cout << "Read PDB file: " << file << std::endl;
      if ( verbose ) {
	clipper::Atom_list atoms = mol.atom_list();
	std::cout << "Number of atoms read: " << atoms.size() << std::endl;
	for ( int i = 0; i < atoms.size(); i += atoms.size()-1 ) printf("%i6  %4s  %8.3f %8.3f %8.3f\n", atoms[i].element().c_str(), atoms[i].coord_orth().x(), atoms[i].coord_orth().y(), atoms[i].coord_orth().z() );
      }
    } catch ( clipper::Message_fatal ) {
      std::cout << "FAILED TO READ PDB FILE: " << file << std::endl;
    }
  }
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
    //double x = clipper::Util::max( clipper::Util::min( r0, 0.995 ), 0.005 );
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
#if defined(_MSC_VER) || defined (WIN32)
  return 0.001*(::rand()%2000-1000);
#else
  return 0.001*(::random()%2000-1000);
#endif
}
