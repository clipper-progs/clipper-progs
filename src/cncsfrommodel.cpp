// Clipper app to get NCS from model
/* Copyright 2008 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-minimol.h>
#include <clipper/clipper-ccp4.h>


clipper::String chain_sequence( const clipper::MPolymer& mp )
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


clipper::RTop_orth superpose( const clipper::MPolymer& mp1, const clipper::MPolymer& mp2, const double& rmsd, const int& nmin )
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


int main( int argc, char** argv )
{
  CCP4Program prog( "cncsfrommodel", "0.1", "$Date: 2008/03/31" );

  // defaults
  clipper::String title;
  clipper::String ippdb    = "NONE";
  double rmsd = 1.0;
  int over = 12;

  // command input
  CCP4CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if ( args[arg] == "-title" ) {
      if ( ++arg < args.size() ) title = args[arg];
    } else if ( args[arg] == "-pdbin" ) {
      if ( ++arg < args.size() ) ippdb = args[arg];
    } else if ( args[arg] == "-max-rmsd" ) {
      if ( ++arg < args.size() ) rmsd = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-min-overlap" ) {
      if ( ++arg < args.size() ) over = clipper::String(args[arg]).i();
    } else {
      std::cout << "Unrecognized:\t" << args[arg] << "\n";
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "Usage: cncsfrompdb\n\t-pdbin <filename>\n\t-max-rmsd <rmsd/A>\n\t-min-overlap <residues>\nDetermine NCS from a protein model.\n";
    exit(1);
  }

  // atomic models
  const int mmdbflags = MMDBF_IgnoreBlankLines | MMDBF_IgnoreDuplSeqNum | MMDBF_IgnoreNonCoorPDBErrors | MMDBF_IgnoreRemarks;
  clipper::MMDBfile mmdb;
  clipper::MiniMol mol;
  mmdb.SetFlag( mmdbflags );
  mmdb.read_file( ippdb );
  mmdb.import_minimol( mol );

  for ( int c1 = 0; c1 < mol.size(); c1++ )
    for ( int c2 = 0; c2 < mol.size(); c2++ )
      if ( c1 != c2 ) {
	clipper::RTop_orth rtop = superpose( mol[c1], mol[c2], rmsd, over );
	if ( !rtop.is_null() ) {
	  clipper::Rotation rot( rtop.rot() );
	  clipper::Euler_ccp4 euler = rot.euler_ccp4();
	  clipper::Coord_orth coord( rtop.trn() );
	  std::cout << std::endl << "NCS operator found relating chains " << mol[c1].id() << " and " <<  mol[c2].id() << std::endl;
	  std::cout << "Euler rotation/deg: " << clipper::Util::rad2d(euler.alpha()) << "," << clipper::Util::rad2d(euler.beta()) << "," << clipper::Util::rad2d(euler.gamma()) << std::endl;
	  std::cout << "Orth translation/A: " << coord.x() << "," << coord.y() << "," << coord.z() << std::endl;
	}
      }
}
