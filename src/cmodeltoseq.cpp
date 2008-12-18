// Clipper app to get sequence from model
/* Copyright 2008 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-minimol.h>
#include <fstream>

extern "C" {
#include <stdlib.h>
}
 

int main( int argc, char** argv )
{
  CCP4Program prog( "cmodeltoseq", "0.1", "$Date: 2008/06/09" );

  // defaults
  clipper::String title;
  clipper::String ippdb    = "NONE";
  clipper::String opseq    = "modeltoseq.seq";

  // command input
  CCP4CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if ( args[arg] == "-title" ) {
      if ( ++arg < args.size() ) title = args[arg];
    } else if ( args[arg] == "-pdbin" ) {
      if ( ++arg < args.size() ) ippdb = args[arg];
    } else if ( args[arg] == "-seqout" ) {
      if ( ++arg < args.size() ) opseq = args[arg];
    } else {
      std::cout << "Unrecognized:\t" << args[arg] << "\n";
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "Usage: cmodeltoseq\n\t-pdbin <filename>\n\t-seqout <filename>\nGet sequence from model\n";
    exit(1);
  }

  // atomic models
  const int mmdbflags = MMDBF_IgnoreBlankLines | MMDBF_IgnoreDuplSeqNum | MMDBF_IgnoreNonCoorPDBErrors | MMDBF_IgnoreRemarks;
  clipper::MMDBfile mmdbwrk;
  clipper::MiniMol molwrk;
  mmdbwrk.SetFlag( mmdbflags );
  mmdbwrk.read_file( ippdb    );
  mmdbwrk.import_minimol( molwrk );

  // assemble sequences
  const char rtype1[21] =
    {  'A',  'R',  'N',  'D',  'C',  'Q',  'E',  'G',  'H',  'I',
       'L',  'K',  'M',  'F',  'P',  'S',  'T',  'W',  'Y',  'V',
       'M'};
  const char rtype3[21][4] =
    {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
     "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
     "MSE"};
  const int ntype = sizeof( rtype1 ) / sizeof( rtype1[0] );
  clipper::String seqfile = "";
  for ( int c = 0; c < molwrk.size(); c++ ) {
    clipper::String id = molwrk[c].id();
    clipper::String seq = "";
    for ( int r = 0; r < molwrk[c].size(); r++ ) {
      char symbol = ' ';
      for ( int t = 0; t < ntype; t++ )
	if ( molwrk[c][r].type() == rtype3[t] )
	  symbol = rtype1[t];
      if ( symbol != ' ' ) seq = seq + symbol;
    }
    if ( seq.length() > 0 )
      seqfile = seqfile + "> " + id + "\n\n" + seq + "\n\n";
  }

  // write file
  std::ofstream file( opseq.c_str() );
  file << seqfile;
  file.close();
}
