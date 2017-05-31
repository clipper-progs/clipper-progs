/*! \file nautilus-util.cpp nautilus library */
/* (C) 2011 Kevin Cowtan & University of York all rights reserved */

#include "nautilus-util.h"

#include <fstream>
extern "C" {
#include <stdlib.h>
}


void NautilusUtil::set_reference( clipper::String& pdb )
{
  const char* clibdptr = getenv( "CLIBD" );
  const char* ccp4ptr = getenv( "CCP4" );
  if ( clibdptr != NULL ) {
    clipper::String clibd( clibdptr );
    clipper::String ccp4( ccp4ptr );
    clipper::String path;
    std::ifstream file;
    if ( pdb == "NONE" ) {
      path = clibd+"/nautilus_lib.pdb";
      file.open( path.c_str(), std::ifstream::in ); file.close();
      if ( !file.fail() ) pdb = path;
    }
    if ( pdb == "NONE" ) {
      path = clibd+"\\nautilus_lib.pdb";
      file.open( path.c_str(), std::ifstream::in ); file.close();
      if ( !file.fail() ) pdb = path;
    }
    if ( pdb == "NONE" ) {
      path = ccp4+"/share/nautilus_lib.pdb";
      file.open( path.c_str(), std::ifstream::in ); file.close();
           if ( !file.fail() ) pdb = path;
    }
    if ( pdb == "NONE" ) {
      path = ccp4+"\\share\\nautilus_lib.pdb";
      file.open( path.c_str(), std::ifstream::in ); file.close();
      if ( !file.fail() ) pdb = path;
    }
    if ( pdb == "NONE" ) 
      clipper::Message::message( clipper::Message_fatal( "No reference data specified and not in $CLIBD" ) );
  } else {
    clipper::Message::message( clipper::Message_fatal( "No reference data specified and $CLIBD not found" ) );
  }
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
  clipper::String msg = log_info( mol );
  std::cout << id << ": " << msg << std::endl;
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


clipper::String NautilusLog::log_info( const clipper::MiniMol& mol )
{
  int nnc(0), nna(0);
  for ( int c = 0; c < mol.size(); c++ )
    if ( !mol[c].exists_property( "NON-NA" ) ) {
      nnc += 1;
      nna += mol[c].size();
    }
  return clipper::String(nna,4) + " nucleic acids built in " +
         clipper::String(nnc,3) + " chains.";
}


void NautilusLog::xml( const clipper::String& file, const clipper::MiniMol& mol )
{
  int nres, nseq, nchn, nmax;
  nchn = mol.size();
  nres = nseq = nmax = 0;
  for ( int c = 0; c < mol.size(); c++ ) {
    if ( !mol[c].exists_property( "NON-NA" ) ) {
      if ( mol[c].size() > nmax ) nmax = mol[c].size();
      for ( int r = 0; r < mol[c].size(); r++ ) {
        if ( mol[c][r].lookup( " C4 ", clipper::MM::ANY ) >= 0 ) nres++;
        if ( mol[c][r].lookup( " C4 ", clipper::MM::ANY ) >= 0 && mol[c][r].lookup( " O4 ", clipper::MM::ANY ) >= 0 ) nseq++;
      }
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


void NautilusLog::profile()
{
  if ( prof.size() > 0 ) {
    std::cout << std::endl << "Profile:" << std::endl;
    for ( int i = 0; i < prof.size(); i++ )
      std::cout << prof[i].first << ": " << clipper::String( prof[i].second, 8 ) << " s" << std::endl;
  }
}
