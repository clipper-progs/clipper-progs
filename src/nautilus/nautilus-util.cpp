/*! \file nautilus-util.cpp nautilus library */
/* (C) 2011 Kevin Cowtan & University of York all rights reserved */

#include "nautilus-util.h"

extern "C" {
#include <stdlib.h>
}


#ifdef NAUTILUS_PROFILE
#include <sys/times.h>
void NautilusLog::log( const clipper::String& id )
{
  int i;
  tms tmst;
  times( &tmst );
  long ticks = sysconf(_SC_CLK_TCK);
  double cpu = double( tmst.tms_utime ) / double( ticks );
  double elapsed = cpu - currentcpu;
  if ( id != "" ) {
    for ( i = 0; i < prof.size(); i++ )
      if ( id == prof[i].first ) break;
    if ( i < prof.size() )
      prof[i].second += elapsed;
    else
      prof.push_back( std::pair<std::string,double>( id, elapsed ) );
  }
  currentcpu = cpu;
}
#else
void NautilusLog::log( const clipper::String& id ) {}
#endif


void NautilusLog::log( const clipper::String& id, const clipper::MiniMol& mol, bool view )
{
  int nnc(0), nna(0);
  for ( int c = 0; c < mol.size(); c++ )
    if ( !mol[c].exists_property( "NON-NA" ) ) {
      nnc += 1;
      nna += mol[c].size();
    }
  std::cout << id << ": " << nna << " nucleic acids built in " << nnc << " chains." << std::endl;
  if ( view ) {
    for ( int c = 0; c < mol.size(); c++ )
      std::cout << mol[c].size() << " ";
    std::cout << std::endl;
    for ( int c = 0; c < mol.size(); c++ ) {
      if ( !mol[c].exists_property( "NON-NA" ) ) {
	for ( int r1 = 0; r1 < mol[c].size()-1; r1++ ) {
	  int r2 = r1 + 1;
	  int a1 = mol[c][r1].lookup( " O3'", clipper::MM::ANY );
	  int a2 = mol[c][r2].lookup( " O5'", clipper::MM::ANY );
	  if ( a1 >= 0 && a2 >= 0 ) {
	    double r = sqrt( ( mol[c][r1][a1].coord_orth() -
			       mol[c][r2][a2].coord_orth() ).lengthsq() );
	    if ( r > 5.0 ) std::cout << "BREAK " << c << " " << r1 << " " << r2 << " " << r << std::endl;
	  }
	}
      }
    }
  }
  log( id );
}


/*
void NautilusLog::xml( const clipper::String& file, const clipper::MiniMol& mol )
{
  int nres, nseq, nchn, nmax;
  nchn = mol.size();
  nres = nseq = nmax = 0;
  for ( int c = 0; c < mol.size(); c++ ) {
    if ( mol[c].size() > nmax ) nmax = mol[c].size();
    for ( int r = 0; r < mol[c].size(); r++ ) {
      if ( mol[c][r].lookup( " CA ", clipper::MM::ANY ) >= 0 ) nres++;
      if ( ProteinTools::residue_index_3( mol[c][r].type() ) >= 0 ) nseq++;
    }
  }

  std::ofstream f;
  f.open( file.c_str(), std::ios::out );
  f << "<NautilusResult>" << std::endl;
  f << "<FragmentsBuilt>" << nchn << "</FragmentsBuilt>" << std::endl;
  f << "<ResiduesBuilt>" << nres << "</ResiduesBuilt>" << std::endl;
  f << "<ResiduesSequenced>" << nseq << "</ResiduesSequenced>" << std::endl;
  f << "<ResiduesLongestFragment>" << nmax << "</ResiduesLongestFragment>" << std::endl;
  f << "</NautilusResult>" << std::endl;
  f.close();
}
*/


void NautilusLog::profile()
{
  if ( prof.size() > 0 ) {
    std::cout << std::endl << "Profile:" << std::endl;
    for ( int i = 0; i < prof.size(); i++ )
      std::cout << prof[i].first << ": " << clipper::String( prof[i].second, 8 ) << " s" << std::endl;
  }
}
