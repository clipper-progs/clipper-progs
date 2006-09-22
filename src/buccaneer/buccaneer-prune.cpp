/*! \file buccaneer-prune.cpp buccaneer library */
/* (C) 2006 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-prune.h"

#include <clipper/clipper-contrib.h>



bool Ca_prune::operator() ( clipper::MiniMol& mol2, const clipper::MiniMol& mol1 ) const
{
  typedef clipper::MMonomer Mm;
  clipper::ftype r2cut = rad_*rad_;

  clipper::Cell       cell = mol1.cell();
  clipper::Spacegroup spgr = mol1.spacegroup();

  // split into separate chains
  clipper::MiniMol molx;
  ProteinTools::chain_tidy( molx, mol1 );

  // now loop over chains
  clipper::Coord_frac cf1, cf2;
  for ( int chn1 = 0; chn1 < molx.size()-1; chn1++ ) {
    for ( int chn2 = chn1; chn2 < molx.size(); chn2++ ) {
      // find any clashing residues between these chains
      for ( int res1 = 0; res1 < molx[chn1].size(); res1++ ) {
	for ( int res2 = 0; res2 < molx[chn2].size(); res2++ ) {
	  if ( chn1 != chn2 || res1 != res2 ) {
	    int a1 = molx[chn1][res1].lookup( " CA ", clipper::MM::ANY );
	    int a2 = molx[chn2][res2].lookup( " CA ", clipper::MM::ANY );
	    if ( a1 >= 0 && a2 >= 0 ) {
	      cf1 = molx[chn1][res1][a1].coord_orth().coord_frac(cell);
	      cf2 = molx[chn2][res2][a2].coord_orth().coord_frac(cell);
	      cf2 = cf2.symmetry_copy_near( spgr, cell, cf1 ) - cf1;
	      if ( cf2.lengthsq(cell) < r2cut ) {
		// clash found: if only one is sequenced, keep that,
		//              otherwise keep the longer.
		int scr1 = molx[chn1].size();
		int scr2 = molx[chn2].size();
		if ( molx[chn1][res1].type() == "---" ) scr1 =  -1000;
		if ( molx[chn1][res1].type() == "UNK" ) scr1 += -1000;
		if ( molx[chn2][res2].type() == "---" ) scr2 =  -1000;
		if ( molx[chn2][res2].type() == "UNK" ) scr2 += -1000;
		if ( scr1 > scr2 ) molx[chn2][res2].set_type( "---" );
		else               molx[chn1][res1].set_type( "---" );
	      }
	    }
	  }
	}
      }
    }
  }

  // eliminate any sequences of less than 5 residues
  mol2 = clipper::MiniMol( spgr, cell );
  clipper::MPolymer mp, mpnull;
  for ( int chn = 0; chn < molx.size(); chn++ ) {
    mp = mpnull;
    for ( int res = 0; res < molx[chn].size(); res++ ) {
      if ( molx[chn][res].type() != "---" ) {
	mp.insert( molx[chn][res] );
      } else {
	if ( mp.size() > 5 ) mol2.insert( mp );
	mp = mpnull;
      }
    }
    if ( mp.size() > 5 ) mol2.insert( mp );
  }

  return true;
}