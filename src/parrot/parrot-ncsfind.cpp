/*! \file parrot-ncs.cpp parrot library */
/* Copyright 2008 Kevin Cowtan & University of York all rights reserved */

#include "parrot-ncsfind.h"

extern "C" {
#include <string.h>
}

clipper::String NCSfind::chain_sequence( const clipper::MPolymer& mp )
{
  const int NTYPE = 27;
  const char rtype1[NTYPE] =
    {  'A',  'R',  'N',  'D',  'C',  'Q',  'E',  'G',  'H',  'I',
       'L',  'K',  'M',  'F',  'P',  'S',  'T',  'W',  'Y',  'V',
       'M',  'a',  'c',  'g',  't',  'u',  '?'};
  const char rtype3[NTYPE][4] =
    {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
     "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
     "MSE","  A","  C","  G","  T","  U","UNK"};
  clipper::String seq = "";
  for ( int res = 0; res < mp.size(); res++ ) {
    char c = ' ';
    for ( int t = 0; t < NTYPE; t++ )
      if ( strncmp( mp[res].type().c_str(), rtype3[t], 3 ) == 0 )
	c = rtype1[t];
    if ( c == ' ' )
      c = char( (mp[res].type()[0] + mp[res].type()[1] + mp[res].type()[2])
		% 128 + 128 );  // use dummy sequence symbols for unknown types
    seq += c;
  }
  return seq;
}


clipper::RTop_orth NCSfind::superpose( const clipper::MPolymer& mp1, const clipper::MPolymer& mp2, const double& rmsd, const int& nmin )
{
  clipper::RTop_orth result = clipper::RTop_orth::null();
  clipper::String seq1 = chain_sequence( mp1 );
  clipper::String seq2 = chain_sequence( mp2 );
  // ensure that '?'s don't match
  for ( int i = 0; i < seq1.size(); i++ ) if ( seq1[i] == '?' ) seq1[i] = '1';
  for ( int i = 0; i < seq2.size(); i++ ) if ( seq2[i] == '?' ) seq2[i] = '2';

  // get the sequence alignment
  clipper::MSequenceAlign align( clipper::MSequenceAlign::LOCAL,
                                 1.0, 0.001, -1.0 );
  std::pair<std::vector<int>,std::vector<int> > valign = align( seq1, seq2 );
  const std::vector<int>& v1( valign.first ), v2( valign.second );

  // reject any bad matches
  int nmat, nmis;
  nmat = nmis = 0;
  for ( int i1 = 0; i1 < seq1.size(); i1++ ) {
    int i2 = v1[i1];
    if ( i2 >= 0 && i2 < seq2.size() )
      if ( isalpha(seq1[i1]) && isalpha(seq2[i2]) ) {
	if ( seq1[i1] == seq2[i2] ) nmat++;
	else                        nmis++;
      }
  }
  if ( nmat < nmin ) return result;

  // now get the coordinates
  std::vector<clipper::Coord_orth> c1, c2;
  for ( int i1 = 0; i1 < seq1.size(); i1++ ) {
    int i2 = v1[i1];
    if ( i2 >= 0 && i2 < seq2.size() ) 
      if ( seq1[i1] == seq2[i2] )
	if ( isalpha(seq1[i1]) ) {
	  int a1 = mp1[i1].lookup( " CA ", clipper::MM::ANY );
	  int a2 = mp2[i2].lookup( " CA ", clipper::MM::ANY );
	  if ( a1 < 0 && a2 < 0 ) {
	    a1 = mp1[i1].lookup( " C1*", clipper::MM::ANY );
	    a2 = mp2[i2].lookup( " C1*", clipper::MM::ANY );
	  }
	  if ( a1 >= 0 && a2 >= 0 ) {
	    c1.push_back( mp1[i1][a1].coord_orth() );
	    c2.push_back( mp2[i2][a2].coord_orth() );
	  }
	}
  }

  // refine the alignment
  clipper::RTop_orth rtop_tmp;
  double r2;
  for ( int c = 0; c < 5; c++ ) {
    int nc = c1.size();
    // get transformation
    rtop_tmp = clipper::RTop_orth( c1, c2 );
    // get rmsd
    std::vector<std::pair<double,int> > r2index( nc );
    r2 = 0.0;
    for ( int i = 0; i < nc; i++ ) {
      double d2 = ( rtop_tmp * c1[i] - c2[i] ).lengthsq();
      r2 += d2;
      r2index[i] = std::pair<double,int>( d2, i );
    }
    r2 /= double( nc );
    // prune the list to improve it
    std::sort( r2index.begin(), r2index.end() );
    std::vector<clipper::Coord_orth> t1, t2;
    for ( int i = 0; i < (9*r2index.size())/10; i++ ) {
      t1.push_back( c1[r2index[i].second] );
      t2.push_back( c2[r2index[i].second] );
    }
    c1 = t1;
    c2 = t2;
  }

  // if a close match has been found, return it
  if ( r2 < rmsd*rmsd ) result = rtop_tmp;
  return result;
}


Local_rtop NCSfind::local_rtop( const clipper::MPolymer& mp1, const clipper::MPolymer& mp2, const double& rmsd, const int& nmin )
{
  clipper::RTop_orth rtop = superpose( mp1, mp2, rmsd, nmin );
  if ( rtop.is_null() ) return Local_rtop( clipper::Rotation::null(), clipper::Coord_orth(), clipper::Coord_orth() );

  clipper::Coord_orth co1(0.0,0.0,0.0), co2(0.0,0.0,0.0);
  double n1(0.0), n2(0.0);
  for ( int r = 0; r < mp1.size(); r++ ) 
    for ( int a = 0; a < mp1[r].size(); a++ ) {
      co1 += mp1[r][a].coord_orth();
      n1  += 1.0;
    }
  for ( int r = 0; r < mp2.size(); r++ ) 
    for ( int a = 0; a < mp2[r].size(); a++ ) {
      co2 += mp2[r][a].coord_orth();
      n2  += 1.0;
    }
  co1 = (1.0/n1) * co1;
  co2 = (1.0/n2) * co2;
  co2 = rtop.inverse() * co2;
  co1 = 0.5 * ( co1 + co2 );
  co2 = rtop * co1;
  clipper::Rotation rot( rtop.rot() );
  Local_rtop result = Local_rtop( rot, co1, co2 );
  return result;
}


/*
 Find NCS candidate operators from atoms
*/
std::vector<Local_rtop> NCSfind::find_ncs_candidates( const clipper::Atom_list& atoms, const std::vector<clipper::String> atominfo, const clipper::Spacegroup& spgr, const clipper::Cell& cell ) const
{
  // constants
  const int nnear = 25;

  // collect sites
  std::vector<CoordDescr> coords;
  std::vector<std::vector<CoordDescr> >
    environs;
  for ( int i = 0; i < atoms.size(); i++ )
    coords.push_back( CoordDescr( atoms[i].coord_orth(), i ) );
  for ( int i = 0; i < coords.size(); i++ )
    environs.push_back( find_environ( coords, spgr, cell, coords[i].coord_orth(), nnear ) );

  // find and evaluate density rotation candidates
  std::vector<std::vector<std::pair<int,int> > > matchtmp;
  std::vector<std::vector<std::pair<CoordDescr,CoordDescr> > > matches, matchescut;
  double rmax = 0.0;
  for ( int c1 = 0; c1 < coords.size(); c1++ ) {
    for ( int c2 = c1+1 ; c2 < coords.size(); c2++ ) {
      const std::vector<CoordDescr>& near1 = environs[c1];
      const std::vector<CoordDescr>& near2 = environs[c2];
      matchtmp = match_atoms( near1, near2 );
      for ( int m = 0; m < matchtmp.size(); m++ ) {
	std::vector<std::pair<CoordDescr,CoordDescr> > match;
	for ( int n = 0; n < matchtmp[m].size(); n++ )
	  match.push_back( std::pair<CoordDescr,CoordDescr>(
	       near1[matchtmp[m][n].first], near2[matchtmp[m][n].second] ) );
	matches.push_back( match );
      }
    }
  }

  // sort matches by score
  std::vector<std::pair<double,int> > index(matches.size());
  for ( int m = 0; m < matches.size(); m++ )
    index[m] = std::pair<double,int>( double( matches[m].size() ), m );
  std::sort( index.begin(), index.end() );
  std::reverse( index.begin(), index.end() );

  // make trimmed list
  double scut = 0.2 * index[0].first + 2.4;
  for ( int i = 0; i < index.size(); i++ ) {
    if ( index[i].first < scut-0.01 ) break;
    matchescut.push_back( matches[index[i].second] );
  }

  // assemble results
  std::vector<Local_rtop> results;
  for ( int m = 0; m < matchescut.size(); m++ )
    results.push_back( local_rtop( matchescut[m] ).proper( spgr, cell ) );

  // remove duplicates
  results = Local_rtop::exclude_identity( results, spgr, cell, 10.0, 0.20 );
  results = Local_rtop::exclude_duplicate( results, spgr, cell, 1.5, 0.05 );

  // print out info
  if ( debug_ ) {
    for ( int m = 0; m < matchescut.size(); m++ ) {
      std::cout << m << " " << matchescut[m].size() << std::endl;
      for ( int n = 0; n < matchescut[m].size(); n++ ) {
	std::cout << "   " << atominfo[matchescut[m][n].first.index()] << " " << atominfo[matchescut[m][n].second.index()] << std::endl;
      }
    }
  }

  return results;
}


/*
 Find NCS candidate operators from atoms
*/
std::vector<Local_rtop> NCSfind::find_ncs_candidates( const clipper::Atom_list& atoms, const clipper::Spacegroup& spgr, const clipper::Cell& cell ) const
{
  std::vector<clipper::String> atominfo( atoms.size(), clipper::String("") );
  return find_ncs_candidates( atoms, atominfo, spgr, cell );
}


/*
 Find atoms near some coordinate
*/
std::vector<NCSfind::CoordDescr> NCSfind::find_environ( std::vector<CoordDescr> coords, const clipper::Spacegroup& spgr, const clipper::Cell& cell, const clipper::Coord_orth& near, const int& nnear ) {
  clipper::Coord_frac cf, cfa, cfs, cft;
  int u,v,w;
  std::vector<CoordDescr> all, allcut;
  // generate list of source and target atoms
  cf = near.coord_frac(cell);
  for ( int i = 0; i < coords.size(); i++ ) {
    clipper::Coord_orth coord = coords[i].coord_orth();
    int info = coords[i].index();
    cfa = coord.coord_frac(cell);
    for ( int sym = 0; sym < spgr.num_symops(); sym++ ) {
      cfs = ( spgr.symop(sym) * cfa ).lattice_copy_near(cf);
      for ( u = -1; u <= 1; u++ )
	for ( v = -1; v <= 1; v++ )
	  for ( w = -1; w <= 1; w++ ) {
	    cft = cfs + clipper::Coord_frac( double(u), double(v), double(w) );
	    all.push_back( CoordDescr( cft.coord_orth(cell), info ) );
	  }
    }
  }
  // now sort on distance
  std::vector<std::pair<double,int> > index( all.size() );
  for ( int i = 0; i < all.size(); i++ ) {
    double d2 = (all[i].coord_orth()-near).lengthsq();
    index[i] = std::pair<double,int>( d2, i );
  }
  std::sort( index.begin(), index.end() );
  // make truncated list
  allcut.resize( clipper::Util::min( nnear, int( all.size() ) ) );
  for ( int i = 0; i < allcut.size(); i++ )
    allcut[i] = all[ index[i].second ];
  // return result
  return allcut;
}


/*
 Find RT operator from atoms and matches
*/
clipper::RTop_orth NCSfind::match_rtop( const std::vector<CoordDescr>& near1, const std::vector<CoordDescr>& near2, std::vector<std::pair<int,int> > match )
{
  std::vector<clipper::Coord_orth> src( match.size() ), tgt( match.size() );
  for ( int i = 0; i < match.size(); i++ ) {
    src[i] = near1[match[i].first ].coord_orth();
    tgt[i] = near2[match[i].second].coord_orth();
  }
  return clipper::RTop_orth( src, tgt );
}


/*
 Find coordinate disagreements from atoms and matches
*/
std::vector<double> NCSfind::match_diff( const std::vector<CoordDescr>& near1, const std::vector<CoordDescr>& near2, std::vector<std::pair<int,int> > match, clipper::RTop_orth rtop )
{
  std::vector<double> diff( match.size() );
  for ( int i = 0; i < match.size(); i++ )
    diff[i] = ( rtop * near1[match[i].first ].coord_orth() - near2[match[i].second].coord_orth() ).lengthsq();
  return diff;
}


/*
 Find NCS operator from pre-matched atoms
*/
Local_rtop NCSfind::local_rtop( const std::vector<std::pair<CoordDescr,CoordDescr> >& coords )
{
  std::vector<clipper::Coord_orth> src( coords.size() ), tgt( coords.size() );
  clipper::Coord_orth com1(0.0,0.0,0.0), com2(0.0,0.0,0.0);
  for ( int i = 0; i < coords.size(); i++ ) {
    src[i] = coords[i].first .coord_orth();
    tgt[i] = coords[i].second.coord_orth();
    com1 += src[i];
    com2 += tgt[i];
  }
  com1 = ( 1.0 / double(coords.size()) ) * com1;
  com2 = ( 1.0 / double(coords.size()) ) * com2;
  clipper::RTop_orth rtop( src, tgt );
  clipper::Rotation rot( rtop.rot() );
  com2 = rtop.inverse() * com2;
  com1 = 0.5 * ( com1 + com2 );
  com2 = rtop * com1;
  return Local_rtop( rot, com1, com2 );
}


/*
Superpose atoms starting from an initial atom
*/
std::vector<std::vector<std::pair<int,int> > > NCSfind::match_atoms( const std::vector<CoordDescr>& near1, const std::vector<CoordDescr>& near2 ) const
{
  const double dcut = 2.0;
  const int max_match = 8;

  // limit maximum number of atoms for initial match
  int nmatch  = clipper::Util::min( near1.size(), near2.size() );
  int nmax = clipper::Util::min( nmatch, max_match );

  // first make a list of rotation candidates
  std::vector<std::pair<int,int> > match(3);
  std::vector<std::vector<std::pair<int,int> > > matches;
  match[0] = std::pair<int,int>( 0, 0 );
  // loop over first source and target atom
  for ( int s1 = 1; s1 < nmax-1; s1++ ) {
    bool uniqs1 = near1[s1].index() != near1[0].index();
    double ls1 = sqrt((near1[s1].coord_orth()-
		       near1[0].coord_orth()).lengthsq());
    for ( int t1 = 1; t1 < nmax; t1++ ) {
      bool uniqt1 = near2[t1].index() != near2[0].index();
      double lt1 = sqrt((near2[t1].coord_orth()-
			 near2[0].coord_orth()).lengthsq());
      if ( uniqs1 && uniqt1 && fabs(lt1-ls1) < dcut ) {
	// loop over second source and target atom
	for ( int s2 = 1; s2 < nmax; s2++ ) if ( s2 >  s1 ) {
	  bool uniqs2 = near1[s2].index() != near1[0].index();
	  double ls2 = sqrt((near1[s2].coord_orth()-
			     near1[0].coord_orth()).lengthsq());
	  for ( int t2 = 1; t2 < nmax; t2++ ) if ( t2 != t1 ) {
	    bool uniqt2 = near2[t2].index() != near2[0].index();
	    double lt2 = sqrt((near2[t2].coord_orth()-
			       near2[0].coord_orth()).lengthsq());
	    if ( uniqs2 && uniqt2 && fabs(lt2-ls2) < dcut ) {
	      // check third edge of triangle
	      bool uniqs12 = near1[s2].index() != near1[s1].index();
	      bool uniqt12 = near2[t2].index() != near2[t1].index();
	      double ls12 = sqrt((near1[s2].coord_orth()-
				  near1[s1].coord_orth()).lengthsq());
	      double lt12 = sqrt((near2[t2].coord_orth()-
				  near2[t1].coord_orth()).lengthsq());
	      if ( uniqs12 && uniqt12 && fabs(ls12-lt12) < dcut ) {
		match[1] = std::pair<int,int>( s1, t1 );
		match[2] = std::pair<int,int>( s2, t2 );
		matches.push_back( match );
	      }
	    }
	  }
	}
      }
    }
  }

  // find biggest atom index...
  int maxindex = 0;
  for ( int i = 0; i < near1.size(); i++ )
    if ( near1[i].index() > maxindex ) maxindex = near1[i].index();
  for ( int i = 0; i < near2.size(); i++ )
    if ( near2[i].index() > maxindex ) maxindex = near2[i].index();

  // now try to extend each match in turn
  for ( int m = 0; m < matches.size(); m++ ) {
    std::vector<std::pair<int,int> > match = matches[m];
    for ( int n = 3; n < nmatch; n++ ) {
      // mark used atoms
      std::vector<bool>	used1( maxindex, false ), used2( maxindex, false );
      for ( int n = 0; n < match.size(); n++ ) {
	used1[near1[match[n].first ].index()] = true;
	used2[near2[match[n].second].index()] = true;
      }
      // get initial superposition
      clipper::RTop_orth rtop = match_rtop( near1, near2, match );
      // try other atoms
      std::pair<int,int> newpair( -1, -1 );
      double d2min = dcut*dcut;
      for ( int n1 = 0; n1 < near1.size(); n1++ )
	if ( !used1[near1[n1].index()] )
	  for ( int n2 = 0; n2 < near2.size(); n2++ )
	    if ( !used2[near2[n2].index()] ) {
	      double d2 = ( rtop * near1[n1].coord_orth() -
			    near2[n2].coord_orth() ).lengthsq();
	      if ( d2 < d2min ) {
		newpair = std::pair<int,int>(n1,n2);
		d2min = d2;
	      }
	    }
      if ( newpair.first < 0 ) break;
      // try the new match
      std::vector<std::pair<int,int> > newmatch = match;
      newmatch.push_back( newpair );
      clipper::RTop_orth newop = match_rtop( near1, near2, match );
      std::vector<double> diff = match_diff( near1, near2, match, newop );
      double mdiff = 0.0;
      for ( int i = 0; i < diff.size()-1; i++ )	mdiff += diff[i];
      mdiff = mdiff / double( diff.size()-1 );
      if ( diff.back() > 2.0 * mdiff ) break;
      // good match found, save it
      match.push_back( newpair );
    }
    // store result
    matches[m] = match;
  }

  // sort matches by score
  std::vector<std::pair<double,int> > index(matches.size());
  for ( int m = 0; m < matches.size(); m++ )
    index[m] = std::pair<double,int>( double( matches[m].size() ), m );
  std::sort( index.begin(), index.end() );
  std::reverse( index.begin(), index.end() );

  // make final list
  std::vector<std::vector<std::pair<int,int> > > result;
  double scut = 0.2 * index[0].first + 2.4;
  for ( int i = 0; i < index.size(); i++ ) {
    if ( index[i].first < scut-0.01 ) break;
    result.push_back( matches[index[i].second] );
  }

  return result;
}
