/*! \file buccaneer-prot.cpp buccaneer library */
/* (C) 2002-2006 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-prot.h"
#include "buccaneer-sequence.h"

#include <clipper/clipper-contrib.h>
extern "C" {
#include <string.h>
}


Ca_group::Ca_group( const clipper::MMonomer& mm )
{
  int index_n  = mm.lookup( " N  ", clipper::MM::ANY );
  int index_ca = mm.lookup( " CA ", clipper::MM::ANY );
  int index_c  = mm.lookup( " C  ", clipper::MM::ANY );
  if ( index_n >= 0 && index_ca >= 0 && index_c >= 0 ) {
    coord_n_  = mm[index_n ].coord_orth();
    coord_ca_ = mm[index_ca].coord_orth();
    coord_c_  = mm[index_c ].coord_orth();
  } else {
    coord_n_ = coord_ca_ = coord_c_ =
      clipper::Coord_orth( clipper::Coord_orth::null() );
  }
}

clipper::RTop_orth Ca_group::rtop_from_std_ori() const
{
  const clipper::Coord_orth c1 = coord_c() - coord_ca();
  const clipper::Coord_orth c2 = coord_n() - coord_ca();
  const clipper::Vec3<> v1( (c1.unit()+c2.unit()).unit() );
  const clipper::Vec3<> v2( clipper::Vec3<>::cross(c1,c2).unit() );
  const clipper::Vec3<> v3( clipper::Vec3<>::cross(v1,v2).unit() );
  return
    clipper::RTop_orth( clipper::Mat33<>( v1[0], v2[0], v3[0],
					  v1[1], v2[1], v3[1],
					  v1[2], v2[2], v3[2] ), coord_ca() );
}

clipper::RTop_orth Ca_group::rtop_beta_carbon() const
{
  clipper::RTop_orth rtop = rtop_from_std_ori();
  return clipper::RTop_orth( rtop.rot(), rtop.trn() + coord_cb() - coord_ca() );
}

Ca_group Ca_group::next_ca_group(const clipper::ftype& psi, const clipper::ftype& phi ) const
{
  const clipper::ftype pi = clipper::Util::pi();
  clipper::Coord_orth n ( coord_n(),  coord_ca(), coord_c(), 1.32, 1.99, psi );
  clipper::Coord_orth ca( coord_ca(), coord_c(),  n,         1.47, 2.15, pi  );
  clipper::Coord_orth c ( coord_c(),  n,          ca,        1.53, 1.92, phi );
  return Ca_group( n, ca, c );
}

Ca_group Ca_group::prev_ca_group(const clipper::ftype& phi, const clipper::ftype& psi ) const
{
  const clipper::ftype pi = clipper::Util::pi();
  clipper::Coord_orth c ( coord_c(),  coord_ca(), coord_n(), 1.32, 2.15, phi );
  clipper::Coord_orth ca( coord_ca(), coord_n(),  c,         1.53, 1.99, pi  );
  clipper::Coord_orth n ( coord_n(),  c,          ca,        1.47, 1.92, psi );
  return Ca_group( n, ca, c );
}

clipper::ftype Ca_chain::ramachandran_phi( const int& resno ) const
{
  const Ca_chain& chain = (*this);
  if ( resno-1 >= 0 && resno < chain.size() )
    return clipper::Coord_orth::torsion( chain[resno-1].coord_c(), chain[resno].coord_n(), chain[resno].coord_ca(), chain[resno].coord_c() );
  else
    return clipper::Util::nan();
}

clipper::ftype Ca_chain::ramachandran_psi( const int& resno ) const
{
  const Ca_chain& chain = (*this);
  if ( resno >= 0 && resno+1 < chain.size() )
    return clipper::Coord_orth::torsion( chain[resno].coord_n(), chain[resno].coord_ca(), chain[resno].coord_c(), chain[resno+1].coord_n() );
  else
    return clipper::Util::nan();
}


Pr_group::Pr_group( const clipper::Coord_orth& ca, const clipper::Coord_orth& c, const clipper::Coord_orth& other, const TYPE& type )
{
  coord_ca_ = ca;
  coord_c_  = c;
  if ( type == CaCN ) {
    coord_n_ = other;
  } else {
    const clipper::Vec3<> v0 = ( ca - c ).unit();
    const clipper::Vec3<> v1 = ( clipper::Vec3<>::cross( v0, clipper::Vec3<>::cross( v0, other - c ) ) ).unit();
    // length 1.32 angle 1.99
    coord_n_ = c + clipper::Coord_orth( -0.537*v0 + 1.206*v1 );
  }
}

clipper::Coord_orth Pr_group::coord_o() const
{
  const clipper::Vec3<> v0 = ( coord_ca() - coord_c() ).unit();
  const clipper::Vec3<> v1 = ( clipper::Vec3<>::cross( v0, clipper::Vec3<>::cross( v0, coord_n_next() - coord_c() ) ) ).unit();
  // length 1.24 angle 2.11
  return coord_c() + clipper::Coord_orth( -0.637*v0 + 1.064*v1 );
}

clipper::Coord_orth Pr_group::coord_ca_next() const
{
  const clipper::Vec3<> v0 = ( coord_c() - coord_n_next() ).unit();
  const clipper::Vec3<> v1 = ( clipper::Vec3<>::cross( v0, clipper::Vec3<>::cross( v0, coord_ca() - coord_c() ) ) ).unit();
  // length 1.47 angle 2.15
  return coord_n_next() + clipper::Coord_orth( -0.805*v0 + 1.230*v1 );
}

clipper::RTop_orth Pr_group::rtop_from_std_ori() const
{
  clipper::Coord_orth ca = coord_ca_next();
  const clipper::Coord_orth cen = 0.5*( coord_ca() + ca );
  const clipper::Coord_orth c1 = ca - cen;
  const clipper::Coord_orth c2 = coord_c() - cen;
  const clipper::Vec3<> v1( c1.unit() );
  const clipper::Vec3<> v2( clipper::Vec3<>::cross(c1,c2).unit() );
  const clipper::Vec3<> v3( clipper::Vec3<>::cross(v1,v2).unit() );
  return
    clipper::RTop_orth( clipper::Mat33<>( v1[0], v2[0], v3[0],
					  v1[1], v2[1], v3[1],
					  v1[2], v2[2], v3[2] ), cen );
}

Pr_group Pr_group::next_pr_group( const clipper::ftype& phi, const clipper::ftype& psi ) const
{
  const clipper::ftype pi = clipper::Util::pi();
  clipper::Coord_orth ca(coord_ca(), coord_c(), coord_n_next(), 1.47,2.15,pi);
  clipper::Coord_orth c (coord_c(), coord_n_next(), ca,         1.53,1.92,phi);
  clipper::Coord_orth n (coord_n_next(), ca, c,                 1.32,1.99,psi);
  return Pr_group( ca, c, n, CaCN );
}

Pr_group Pr_group::prev_pr_group( const clipper::ftype& psi, const clipper::ftype& phi ) const
{
  const clipper::ftype pi = clipper::Util::pi();
  clipper::Coord_orth n (coord_n_next(), coord_c(), coord_ca(), 1.47,1.92,psi);
  clipper::Coord_orth c (coord_c(), coord_ca(), n,              1.32,2.15,phi);
  clipper::Coord_orth ca(coord_ca(), n, c,                      1.53,1.99,pi);
  return Pr_group( ca, c, n, CaCN );
}


// Protein loops

ProteinLoop::ProteinLoop( int torsion_sampling )
{
  rama = clipper::Ramachandran( clipper::Ramachandran::All );
  ntor = torsion_sampling;
}


clipper::Coord_orth ProteinLoop::Coord_O( const clipper::Coord_orth ca0, const clipper::Coord_orth c0, const clipper::Coord_orth n1 ) const
{
  clipper::Vec3<> v = ( (c0-n1).unit() + (c0-ca0).unit() ).unit();
  return clipper::Coord_orth( c0 + 1.23*v );
}


clipper::Coord_orth ProteinLoop::Coord_Cb( const clipper::Coord_orth n0, const clipper::Coord_orth ca0, const clipper::Coord_orth c0 ) const
{
  clipper::Vec3<> v1 = ( (ca0-n0).unit() + (ca0-c0).unit() ).unit();
  clipper::Vec3<> v2 = clipper::Coord_orth::cross( (ca0-n0), (ca0-c0) );
  return clipper::Coord_orth( ca0 + 1.04*v1 + 0.53*v2 );
}


/* torsion-rebuild 5 atoms */
std::vector<ProteinLoop::CoordList<5> > ProteinLoop::rebuild5atoms( const clipper::Coord_orth c0, const clipper::Coord_orth n1, const clipper::Coord_orth ca1, const clipper::Coord_orth ca3, const clipper::Coord_orth c3, const clipper::Coord_orth n4 ) const
{
  using clipper::Coord_orth;
  const double dtor = clipper::Util::twopi() / double( ntor );
  std::vector<CoordList<5> > result;
  CoordList<5> r;
  for ( int iphi1 = 0; iphi1 < ntor; iphi1++ ) {
    // build c1
    double phi1 = dtor * double( iphi1 );
    r[0] = Coord_orth( c0, n1, ca1, 1.53, 1.92, phi1 );
    // find 2 possible positions for ca2
    std::vector<Coord_orth> ca2s = constrained_coords( ca1, r[0]-ca1, 3.8, 0.360, ca3, 3.8 );
    for ( int ipsi1 = 0; ipsi1 < ca2s.size(); ipsi1++ ) {
      r[2] = ca2s[ipsi1];
      double psi1 = Coord_orth::torsion( n1, ca1, r[0], r[2] );
      if ( rama.allowed( phi1, psi1 ) ) {
	// build n2
	r[1] = Coord_orth( n1, ca1, r[0], 1.33, 1.99, psi1 );
	// find 2 possible positions for c2
	std::vector<Coord_orth> cc2s = constrained_coords( r[2], r[2]-r[1], 1.53, 1.22, ca3, 2.43 );
	for ( int iphi2 = 0; iphi2 < cc2s.size(); iphi2++ ) {
	  // build c2, n3
	  r[3] = cc2s[iphi2];
	  double phi2 = Coord_orth::torsion( r[0], r[1], r[2], r[3] );
	  double psi2 = Coord_orth::torsion( r[1], r[2], r[3], ca3 );
	  if ( rama.allowed( phi2, psi2 ) ) {
	    r[4] = Coord_orth( r[1], r[2], r[3], 1.33, 1.99, psi2 );
	    double ang = Coord_orth::angle( r[4], ca3, c3 );
	    if ( ang > 1.75 && ang < 2.10 ) result.push_back( r );
	  }
	}
      }
    }
  }
  return result;
}


/* torsion-rebuild 8 atoms */
std::vector<ProteinLoop::CoordList<8> > ProteinLoop::rebuild8atoms( const clipper::Coord_orth c0, const clipper::Coord_orth n1, const clipper::Coord_orth ca1, const clipper::Coord_orth ca4, const clipper::Coord_orth c4, const clipper::Coord_orth n5 ) const
{
  using clipper::Coord_orth;
  const double dtor = clipper::Util::twopi() / double( ntor );
  const double pi = clipper::Util::pi();
  std::vector<CoordList<8> > result;
  CoordList<8> r;
  std::vector<CoordList<5> > r5;
  for ( int iphi1 = 0; iphi1 < ntor; iphi1++ ) {
    // build c1
    double phi1 = dtor * double( iphi1 );
    for ( int ipsi1 = 0; ipsi1 < ntor; ipsi1++ ) {
      double psi1 = dtor * double( ipsi1 );
      if ( rama.allowed( phi1, psi1 ) ) {
	r[0] = Coord_orth(  c0,   n1,  ca1, 1.53, 1.92, phi1 );
	r[1] = Coord_orth(  n1,  ca1, r[0], 1.33, 1.99, psi1 );
	r[2] = Coord_orth( ca1, r[0], r[1], 1.47, 2.15, pi   );
	r5 = rebuild5atoms( r[0], r[1], r[2], ca4, c4, n5 );
	for ( int i = 0; i < r5.size(); i++ ) {
	  for ( int j = 0; j < 5; j++ ) r[j+3] = r5[i][j];
	  result.push_back(r);
	}
      }
    }
  }
  return result;
}


/* find up to 2 possible positions obtained by rotating an arm of
   given base, length and angle to the axis about a rotation axis, to
   give a point a given distance from a target point */
std::vector<clipper::Coord_orth> ProteinLoop::constrained_coords( const clipper::Coord_orth& srcpos, const clipper::Coord_orth& rtnvec, const double& length, const double& angle, const clipper::Coord_orth& tgtpos, const double& tgtdst ) const
{
  using clipper::Coord_orth;
  std::vector<Coord_orth> result;
  Coord_orth v0( rtnvec.unit() );
  Coord_orth v1( Coord_orth::cross( v0, tgtpos-srcpos ).unit() );
  Coord_orth v2( Coord_orth::cross( v1, v0 ).unit() );
  double lcos = length * cos(angle);
  double lsin = length * sin(angle);
  Coord_orth p0 = lcos*v0 + srcpos;
  double z = fabs( (tgtpos-p0) * v0 );
  double x = fabs( (tgtpos-p0) * v2 );
  if ( tgtdst <= z ) return result;
  double d2 = tgtdst*tgtdst - z*z;
  double cost = (lsin*lsin + x*x - d2)/(2.0*lsin*x);
  if ( cost*cost <= 0.995 ) {
    double sint = sqrt(1.0-cost*cost);
    result.push_back( p0 + (lsin*cost)*v2 - (lsin*sint)*v1 );
    result.push_back( p0 + (lsin*cost)*v2 + (lsin*sint)*v1 );
  } else if ( cost*cost <= 1.1 ) {
    result.push_back( p0 + (lsin*cost)*v2 );
  }
  return result;
}


// Protein tools

clipper::String ProteinTools::chain_sequence( const clipper::MPolymer& mp )
{
  clipper::String seq = "";
  for ( int res = 0; res < mp.size(); res++ ) {
    if        ( mp[res].type() == "+++" ) {
      seq += "+";
    } else if ( mp[res].type() == "---" ) {
      seq += "-";
    } else {
      int t = ProteinTools::residue_index_3( mp[res].type() );
      if ( t >= 0 ) seq += ProteinTools::residue_code_1( t );
      else          seq += "?";
    }
  }
  return seq;
}


std::pair<int,int> ProteinTools::chain_sequence_match( const clipper::String& chnseq, const clipper::MMoleculeSequence& seq )
{
  // convert sequences to unique strings
  std::vector<clipper::String> seqs( seq.size() );
  for ( int chn = 0; chn < seq.size(); chn++ ) {
    clipper::String s = "";
    for ( int res = 0; res < seq[chn].sequence().length(); res++ )
      s += residue_code_1(residue_index_1(seq[chn].sequence().substr(res,1)));
    seqs[chn] = s;
  }
  // set minimum match threshold
  int minscr = 0;
  for ( int i = 0; i < chnseq.size(); i++ )
    if ( isupper(chnseq[i]) ) minscr++;
  minscr = minscr/3+4;
  // now find best match
  int bestchn = -1;
  int bestoff = -1;
  int bestscr = minscr;
  int lenc = chnseq.length();
  for ( int chn = 0; chn < seqs.size(); chn++ ) {
    int lens = seqs[chn].length();
    for ( int off = -lenc+bestscr; off < lens-bestscr; off++ ) {
      int scr = 0;
      for ( int i = 0; i < seqs[chn].length(); i++ )
	if ( i-off >= 0 && i-off < chnseq.length() )
	  if ( seqs[chn][i] == chnseq[i-off] )
	    if ( isupper(chnseq[i-off]) )
	      scr++;
      if ( scr > bestscr ) {
	bestchn = chn;
	bestoff = off;
	bestscr = scr;
      }
    }
  }
  // if ( bestchn >= 0 ) { clipper::String tmp = "???????????????????????????????????????????????????????????????????????????????????????????????????"+seqs[bestchn]+"???????????????????????????????????????????????????????????????????????????????????????????????????"; std::cout << tmp.substr( bestoff+99, chnseq.length() ) << "\n"; }
  // return best match
  return std::pair<int,int>( bestchn, bestoff );
}


bool ProteinTools::chain_renumber( clipper::MPolymer& pol, const clipper::MMoleculeSequence& seq )
{
  // convert sequences to unique strings
  clipper::String chnseq = chain_sequence( pol );
  std::vector<clipper::String> seqs( seq.size() );
  for ( int chn = 0; chn < seq.size(); chn++ ) {
    clipper::String s = "";
    for ( int res = 0; res < seq[chn].sequence().length(); res++ )
      s += residue_code_1(residue_index_1(seq[chn].sequence().substr(res,1)));
    seqs[chn] = s;
  }

  // now find best match
  int bestchn = -1;
  int bestscr = 0;
  clipper::MSequenceAlign align( clipper::MSequenceAlign::LOCAL,
				 1.0, -1000.0, -4.0 );
  std::pair<std::vector<int>,std::vector<int> > result;
  for ( int chn = 0; chn < seqs.size(); chn++ ) {
    const clipper::String& seqseq = seqs[chn];
    result = align( chnseq, seqseq );
    int scr = 0;
    for ( int i = 0; i < result.first.size(); i++ ) {
      if ( result.first[i] >= 0 )
	if ( chnseq[i] == seqseq[result.first[i]] )
	  scr++;
	else
	  scr--;
    }
    if ( scr > bestscr ) {
      bestchn = chn;
      bestscr = scr;
    }
  }
  if ( bestchn < 0 ) return false;

  // now number residues
  clipper::String truseq = seqs[bestchn];
  result = align( chnseq, truseq );
  std::vector<int> nums = result.first;

  /*
  std::cout << bestchn << " " << bestscr << std::endl;
  for ( int i = 0; i < chnseq.size(); i++ )
    std::cout << chnseq[i];
  std::cout << std::endl;
  for ( int i = 0; i < chnseq.size(); i++ )
    std::cout << (nums[i]>=0&&nums[i]<truseq.length() ? truseq[nums[i]] : '-');
  std::cout << std::endl;
  */

  // find bounds of sequenced region
  int i0, i1, i;
  for ( i0 = 0;    i0 < nums.size(); i0++ ) if ( nums[i0] >= 0 ) break;
  for ( i1 = nums.size()-1; i1 >= 0; i1-- ) if ( nums[i1] >= 0 ) break;
  if ( i0 < nums.size() )
    for ( i = i0 - 1; i >= 0;          i-- ) nums[i] = nums[i+1] - 1;
  if ( i1 >= 0 )
    for ( i = i1 + 1; i < nums.size(); i++ ) nums[i] = nums[i-1] + 1;

  // renumber the model, with inscodes if necessary
  const clipper::String inscodes = " ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  const int l = inscodes.length();
  for ( i = 0; i < nums.size(); i++ ) nums[i] = nums[i] * l;
  for ( i = 1; i < nums.size(); i++ )
    if ( nums[i] <= nums[i-1] ) nums[i] = nums[i-1]+1;
  for ( i = 0; i < pol.size(); i++ ) {
    const int nres = (nums[i]/l) + 1;
    const int nins = clipper::Util::mod( nums[i], l );
    if ( nins == 0 ) pol[i].set_seqnum( nres );
    else             pol[i].set_seqnum( nres, inscodes.substr( nins, 1 ) );
  }

  /*
  for ( int i = 0; i < chnseq.size()-1; i++ )
    if ( nums[i+1] != nums[i]+1 ) std::cout << "! " << pol.id() << " " << nums[i] << " " << nums[i+1] << std::endl;
  */

  return true;
}


clipper::RTop_orth ProteinTools::superpose( const clipper::MPolymer& mp1, const clipper::MPolymer& mp2, const double& rmsd, const int& nmatch, const int& nmismatch )
{
  clipper::RTop_orth result = clipper::RTop_orth::null();
  clipper::String seq1 = ProteinTools::chain_sequence( mp1 );
  clipper::String seq2 = ProteinTools::chain_sequence( mp2 );
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
      if ( isupper(seq1[i1]) && isupper(seq2[i2]) ) {
	if ( seq1[i1] == seq2[i2] ) nmat++;
	else                        nmis++;
      }
  }
  if ( nmat < nmatch || nmis > nmismatch ) return result;

  /*
  std::cout << seq1 << std::endl;
  for ( int i = 0; i < seq1.size(); i++ ) if ( v1[i] >= 0 ) std::cout << "|"; else std::cout << " "; 
  std::cout << std::endl;
  for ( int i = 0; i < seq2.size(); i++ ) if ( v2[i] >= 0 ) std::cout << "|"; else std::cout << " "; 
  std::cout << std::endl;
  std::cout << seq2 << std::endl;
  */

  // now get the coordinates
  std::vector<clipper::Coord_orth> c1, c2;
  for ( int i1 = 0; i1 < seq1.size(); i1++ ) {
    int i2 = v1[i1];
    if ( i2 >= 0 && i2 < seq2.size() )
      if ( isupper(seq1[i1]) && isupper(seq2[i2]) && seq1[i1] == seq2[i2] ) {
	int a1 = mp1[i1].lookup( " CA ", clipper::MM::ANY );
	int a2 = mp2[i2].lookup( " CA ", clipper::MM::ANY );
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


bool ProteinTools::chain_tidy( clipper::MiniMol& mol )
{
  typedef clipper::MMonomer Mm;
  // create new minimol
  clipper::MiniMol target( mol.spacegroup(), mol.cell() );
  // now separate unlinked fragments into separate chains
  clipper::MPolymer mp, mpnull;
  for ( int chn = 0; chn < mol.size(); chn++ ) {
    mp = mpnull;
    for ( int res = 0; res < mol[chn].size(); res++ ) {
      mp.insert( mol[chn][res] );
      if ( res < mol[chn].size()-1 )
	if ( !Mm::protein_peptide_bond(mol[chn][res],mol[chn][res+1]) ) {
	  target.insert( mp );
	  mp = mpnull;
	}
    }
    if ( mp.size() > 0 ) target.insert( mp );
  }

  // sort the chains by size
  std::vector<std::pair<int,int> > chnsiz( target.size() ); 
  for ( int chn = 0; chn < target.size(); chn++ )
    chnsiz[chn] = std::pair<int,int>( -target[chn].size(), chn );
  std::sort( chnsiz.begin(), chnsiz.end() );
  clipper::MiniMol temp = target;
  for ( int chn = 0; chn < target.size(); chn++ )
    target[chn] = temp[chnsiz[chn].second];

  // now relabel the chains and residues
  clipper::String labels1 = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  clipper::String labels2 = "abcdefghijklmnopqrstuvwxyz";
  for ( int chn = 0; chn < target.size(); chn++ ) {
    if ( chn < labels1.length() ) {
      target[chn].set_id( labels1.substr( chn, 1 ) );
      for ( int res = 0; res < target[chn].size(); res++ )
	target[chn][res].set_seqnum( res + 1 );
    } else {
      int c = chn - labels1.length();
      int c1 = c % labels2.length();
      int c2 = c / labels2.length();
      target[chn].set_id( labels2.substr( c1, 1 ) );
      for ( int res = 0; res < target[chn].size(); res++ )
	target[chn][res].set_seqnum( res + 1 + 1000*c2 );
    }
  }

  mol = target;
  return true;
}


bool ProteinTools::copy_residue_types( clipper::MiniMol& target, const clipper::MiniMol& source )
{
  const clipper::Spacegroup& spgr = target.spacegroup();
  const clipper::Cell&       cell = target.cell();

  // make Ca-only copy of source molecule
  clipper::MiniMol sourca( source.spacegroup(), source.cell() );
  sourca.model() = source.select( "*/*/ CA ", clipper::MM::ANY );

  // non-bond search
  const double dmin  = 1.0;
  const double d2min = 1.0;
  clipper::MAtomNonBond nb( sourca, 3.0 );
  std::vector<clipper::MAtomIndexSymmetry> atoms;
  clipper::Coord_frac f1, f2;
  // search over target atoms
  for ( int c2 = 0; c2 < target.size(); c2++ )
    for ( int r2 = 0; r2 < target[c2].size(); r2++ ) {
      int a2 = target[c2][r2].lookup( " CA ", clipper::MM::ANY );
      if ( a2 >= 0 ) {
	// find source atom near target atom
	f2 = target[c2][r2][a2].coord_orth().coord_frac(cell);
	atoms = nb.atoms_near( f2.coord_orth(cell), dmin );
	for ( int i = 0; i < atoms.size(); i++ ) {
	  const int c1 = atoms[i].polymer();
	  const int r1 = atoms[i].monomer();
	  const int a1 = atoms[i].atom();
	  if ( sourca[c1][r1][a1].id() == " CA " ) {
	    f1 = sourca[c1][r1][a1].coord_orth().coord_frac(cell);
	    f1 = spgr.symop(atoms[i].symmetry()) * f1;
	    f1 = f1.lattice_copy_near( f2 );
	    double d2 = ( f2 - f1 ).lengthsq( cell );
	    // if distance < 1A
	    if ( d2 < d2min ) {
	      // copy residue type
	      if ( sourca[c1][r1].type() != "UNK" &&
		   target[c2][r2].type() == "UNK" )
		target[c2][r2].set_type( sourca[c1][r1].type() );
	      // copy sequence data
	      if ( sourca[c1][r1].exists_property( "SEQDAT" ) &&
		   !target[c2][r2].exists_property( "SEQDAT" ) )
		target[c2][r2].set_property( "SEQDAT", static_cast<const clipper::Property<Ca_sequence::Sequence_data>&>( sourca[c1][r1].get_property( "SEQDAT" ) ) );
	    }
	  }
	}
      }
    }

  // remove any short sequence fragments
  for ( int chn = 0; chn < target.size(); chn++ ) {
    int res = 0;
    while( res < target[chn].size() ) {
      while( res < target[chn].size() ) {
	if ( target[chn][res].type() != "UNK" ) break;
	res++;
      }
      int res0 = res;
      while( res < target[chn].size() ) {
	if ( target[chn][res].type() == "UNK" ) break;
	res++;
      }
      int res1 = res;
      if ( res1 - res0 < 6 )
	for ( int r = res0; r < res1; r++ )
	  target[chn][r].set_type( "UNK" );
    }
  }

  return true;
}


bool ProteinTools::globularise( clipper::MiniMol& mol, const clipper::Coord_frac cent )
{
  typedef clipper::MMonomer Mm;
  const clipper::Spacegroup& spgr = mol.spacegroup();
  const clipper::Cell&       cell = mol.cell();
  clipper::Coord_frac cf1, cf2, cf3, cfmin;

  for ( int chn = 0; chn < mol.size(); chn++ ) {
    int res0 = 0;
    for ( int res = 0; res < mol[chn].size(); res++ ) {
      bool chnbrk;
      if ( res == mol[chn].size()-1 )
	chnbrk = true;
      else
	chnbrk = !Mm::protein_peptide_bond( mol[chn][res], mol[chn][res+1] );
      if ( chnbrk ) {
	// process res0...res
	// get COM of current chain
	clipper::Coord_orth co( 0.0, 0.0, 0.0 );
	double no = 0.0;
	for ( int r1 = res0; r1 <= res; r1++ )
	  for ( int a1 = 0; a1 < mol[chn][r1].size(); a1++ ) {
	    co += mol[chn][r1][a1].coord_orth();
	    no += 1.0;
	  }
	co = (1.0/no) * co;
	cf1 = co.coord_frac( cell );
	// find symop and lattice shift which brings it close to cent
	double r2min = 1.0e9;
	int symin = 0;
	cfmin = clipper::Coord_frac( 0.0, 0.0, 0.0 );
	for ( int s = 0; s < spgr.num_symops(); s++ ) {
	  cf2 = spgr.symop(s) * cf1;
	  cf3 = cf2.lattice_copy_near( cent );
	  double r2 = (cf3-cent).lengthsq(cell);
	  if ( r2 < r2min ) {
	    r2min = r2;
	    symin = s;
	    cfmin = cf3 - cf2;
	  }
	}
	// apply the shifts
	for ( int r1 = res0; r1 <= res; r1++ )
	  for ( int a1 = 0; a1 < mol[chn][r1].size(); a1++ ) {
	    cf1 = mol[chn][r1][a1].coord_orth().coord_frac( cell );
	    cf1 = spgr.symop(symin) * cf1 + cfmin;
	    mol[chn][r1][a1].set_coord_orth( cf1.coord_orth( cell ) );
	  }
	// done processing
	res0 = res+1;
      }
    }
  }
  return true;
}

bool ProteinTools::globularise( clipper::MiniMol& mol )
{
  // iteratively globularise the model
  for ( int i = 0; i < 3; i++ ) {
    // find centre of mass
    clipper::Coord_orth co( 0.0, 0.0, 0.0 );
    double no = 0.0;
    for ( int chn = 0; chn < mol.size(); chn++ )
      for ( int res = 0; res < mol[chn].size(); res++ )
	for ( int atm = 0; atm < mol[chn][res].size(); atm++ ) {
	  co += mol[chn][res][atm].coord_orth();
	  no += 1.0;
	}
    co = (1.0/no) * co;
    clipper::Coord_frac cf = co.coord_frac( mol.cell() );
    // globularise about it
    globularise( mol, cf );
  }
  return true;
}

bool ProteinTools::symm_match( clipper::MiniMol& molwrk, const clipper::MiniMol& molref )
{
  clipper::Spacegroup spg1 = clipper::Spacegroup(clipper::Spacegroup::P1);
  clipper::Spacegroup spgr = molwrk.spacegroup();
  clipper::Cell       cell = molwrk.cell();

  // calculate extent of model
  clipper::Atom_list atomr = molref.atom_list();
  clipper::Range<clipper::ftype> urange, vrange, wrange;
  clipper::Coord_frac cfr( 0.0, 0.0, 0.0 );
  for ( int i = 0; i < atomr.size(); i++ ) {
    clipper::Coord_frac cf = atomr[i].coord_orth().coord_frac( cell );
    cfr += cf;
    urange.include( cf.u() );
    vrange.include( cf.v() );
    wrange.include( cf.w() );
  }
  clipper::Coord_frac cf0( urange.min(), vrange.min(), wrange.min() );
  clipper::Coord_frac cf1( urange.max(), vrange.max(), wrange.max() );
  cfr = (1.0/double(atomr.size())) * cfr;

  // calculate mask using wrk cell and ref atoms
  clipper::Resolution reso( 5.0 );
  clipper::Grid_sampling grid( spg1, cell, reso );
  clipper::Grid_range    grng( grid,  cf0,  cf1 );
  grng.add_border(4);
  clipper::NXmap<float> nxmap( cell, grid, grng ), nxflt( cell, grid, grng );
  clipper::EDcalc_mask<float> maskcalc( 2.0 );
  nxmap = 0.0;
  maskcalc( nxmap, atomr );
  MapFilterFn_g5 fn;
  clipper::MapFilter_fft<float>
    fltr( fn, 1.0, clipper::MapFilter_fft<float>::Relative );
  fltr( nxflt, nxmap );

  // now score each chain, symmetry and offset in turn
  for ( int c = 0; c < molwrk.size(); c++ ) {
    double              bestscr = 0.0;
    int                 bestsym = 0;
    clipper::Coord_frac bestoff( 0.0, 0.0, 0.0 );
    const clipper::Coord_frac cfh( 0.5, 0.5, 0.5 );
    for ( int sym = 0; sym < spgr.num_symops(); sym++ ) {
      clipper::Atom_list atomw = molwrk[c].atom_list();
      clipper::RTop_orth rtop = spgr.symop(sym).rtop_orth( cell );
      clipper::Coord_orth cow( 0.0, 0.0, 0.0 );
      for ( int a = 0; a < atomw.size(); a++ ) {
	atomw[a].transform( rtop );
	cow += atomw[a].coord_orth();
      }
      if ( atomw.size() > 0 ) cow = (1.0/double(atomw.size())) * cow;
      clipper::Coord_frac cfw = cow.coord_frac( cell );
      clipper::Coord_frac cfwt = cfw.lattice_copy_near( cfr - cfh );
      clipper::Coord_frac off0 = cfwt - cfw;

      // try offsets
      for ( double du = 0.0; du <= 1.01; du += 1.0 )
	for ( double dv = 0.0; dv < 1.01; dv += 1.0 )
	  for ( double dw = 0.0; dw < 1.01; dw += 1.0 ) {
	    clipper::Coord_frac off( rint( off0.u() ) + du,
				     rint( off0.v() ) + dv,
				     rint( off0.w() ) + dw );
	    clipper::Coord_orth ofo = off.coord_orth( cell );
	    double scr = 0.0;
	    for ( int a = 0; a < atomw.size(); a++ ) {
	      clipper::Coord_orth coa = atomw[a].coord_orth() + ofo;
	      clipper::Coord_grid cga = nxflt.coord_map( coa ).coord_grid();
	      if ( nxflt.in_map( cga ) ) scr += nxflt.get_data( cga );
	    }
	    if ( scr > bestscr ) {
	      bestscr = scr;
	      bestsym = sym;
	      bestoff = off;
	    }
	  }
    }
    // now transform using the best operator
    clipper::Coord_orth cot = bestoff.coord_orth( cell );
    clipper::RTop_orth rtop = spgr.symop(bestsym).rtop_orth( cell );
    rtop = clipper::RTop_orth( rtop.rot(), rtop.trn()+cot );
    molwrk[c].transform( rtop );
  }

  return true;
}


std::vector<float> ProteinTools::main_chain_densities( const clipper::MPolymer& mp, const clipper::Xmap<float>& xmap, int nsmooth )
{
  const char atms[3][5] = { " N  ", " CA ", " C  " };
  // calculate mean density for each residue
  std::vector<float> scores( mp.size(), 0.0 );
  double t0, t1;
  for ( int res = 0; res < mp.size(); res++ ) {
    t0 = t1 = 0.0;
    for ( int i = 0; i < 3; i++ ) {
      int atm = mp[res].lookup( atms[i], clipper::MM::ANY );
      if ( atm >= 0 ) {
	clipper::Coord_orth co = mp[res][atm].coord_orth();
	clipper::Coord_frac cf = co.coord_frac( xmap.cell() );
	double r = xmap.interp<clipper::Interp_cubic>( cf );
	t0 += 1.0;
	t1 += r;
      }
    }
    if ( t0 > 0.5 ) scores[res] = t1 / t0;
  }

  // smooth if required
  for ( int i = 0; i < nsmooth; i++ ) {
    int n = scores.size();
    std::vector<float> smooth( n );
    smooth[ 0 ] = 0.25*( 3.0*scores[0] + scores[1] );
    for ( int i = 1; i < n-1; i++ )
      smooth[i] = 0.25*( scores[i-1] + 2.0*scores[i] + scores[i+1] );
    smooth[n-1] = 0.25*( scores[n-2] + 3.0*scores[n-1] );
    scores = smooth;
  }

  return scores;
}


bool ProteinTools::break_chains( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap )
{
  // separate broken chains
  for ( int chn = 0; chn < mol.size(); chn++ ) {
    // score the residues
    std::vector<float> scores = main_chain_densities( mol[chn], xmap, 5 );
    // mark residues in sequence breaks
    int res = 0;
    while( res < mol[chn].size() ) {
      if ( mol[chn][res].type() != "UNK" ) break;
      res++;
    }
    while( res < mol[chn].size() ) {
      while( res < mol[chn].size() ) {
	if ( mol[chn][res].type() == "UNK" ) break;
	res++;
      }
      int res0 = res;
      while( res < mol[chn].size() ) {
	if ( mol[chn][res].type() != "UNK" ) break;
	res++;
      }
      int res1 = res;
      if ( res1 < mol[chn].size() ) {
	int resm = res0;  // find the weakest residues to delete
	for ( int r = res0; r < res1; r++ )
	  if ( scores[r] < scores[resm] ) resm = r;
	mol[chn][resm].set_type( "~~~" );
	if ( resm-1 >= res0 ) mol[chn][resm-1].set_type( "~~~" );
	if ( resm+1 <  res1 ) mol[chn][resm+1].set_type( "~~~" );
      }
    }
  }
  // eliminate any sequences of less than 6 residues
  clipper::MiniMol mol2( mol.spacegroup(), mol.cell() );
  clipper::MPolymer mp, mpnull;
  for ( int chn = 0; chn < mol.size(); chn++ ) {
    mp = mpnull;
    for ( int res = 0; res < mol[chn].size(); res++ ) {
      if ( mol[chn][res].type() != "~~~" ) {
	mp.insert( mol[chn][res] );
      } else {
	if ( mp.size() > 5 ) mol2.insert( mp );
	mp = mpnull;
      }
    }
    if ( mp.size() > 5 ) mol2.insert( mp );
  }
  mol = mol2;

  return true;
}


const int ProteinTools::ntype = 21;
const char ProteinTools::rtype1[21] =
 {  'A',  'R',  'N',  'D',  'C',  'Q',  'E',  'G',  'H',  'I',
    'L',  'K',  'M',  'F',  'P',  'S',  'T',  'W',  'Y',  'V',
    'M'};
const char ProteinTools::rtype3[21][4] =
 {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
  "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
  "MSE"};
const int ProteinTools::tindex[21] =
 {    0,    1,    2,    3,    4,    5,    6,    7,    8,    9,
     10,   11,   12,   13,   14,   15,   16,   17,   18,   19,
     12};
int ProteinTools::rindex[256], ProteinTools::rindext[256];

ProteinTools::ProteinTools() {
  for ( int i = 0; i < 256; i++ ) {
    clipper::String s = "  ";
    s[0] = char(i);
    rindex [i] = residue_index_1( s, false );
    rindext[i] = residue_index_1( s, true );
  }
}

int ProteinTools::residue_index_1( clipper::String code, bool translate )
{
  int r = -1;
  for ( r = 0; r < ntype; r++ )
    if ( code[0] == rtype1[r] ) break;
  if ( r == ntype ) return -1;
  if ( translate ) r = tindex[r];
  return r;
}

int ProteinTools::residue_index_3( clipper::String code, bool translate )
{
  int r = -1;
  for ( r = 0; r < ntype; r++ )
    if ( strncmp( code.c_str(), rtype3[r], 3 ) == 0 ) break;
  if ( r == ntype ) return -1;
  if ( translate ) r = tindex[r];
  return r;
}

int ProteinTools::residue_index( clipper::String code, bool translate )
{
  if ( code.length() == 3 )
    return residue_index_3( code, translate );
  if ( code.length() == 1 )
    return residue_index_1( code, translate );
  return -1;
}

clipper::String ProteinTools::residue_code_1( int index )
{
  clipper::String result;
  if ( index >= 0 && index < ntype ) result += rtype1[index];
  return result;
}

clipper::String ProteinTools::residue_code_3( int index )
{
  clipper::String result;
  if ( index >= 0 && index < ntype ) result = rtype3[index];
  return result;
}

clipper::String ProteinTools::residue_code( clipper::String code, bool translate )
{
  if ( code.length() == 1 )
    return residue_code_3( residue_index_1( code, translate ) );
  if ( code.length() == 3 )
    return residue_code_1( residue_index_3( code, translate ) );
  return "";
}
