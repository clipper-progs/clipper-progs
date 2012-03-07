// Clipper privateer
/* Copyright 2003-2009 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-minimol.h>

#include <algorithm>

#include "protein_db.h"
#include "protein_db_utils.h"

int main( int argc, char** argv )
{
  CCP4Program prog( "csloop", "0.5", "$Date: 2010/10/05" );

  std::cout << "\nCopyright 2009 Kevin Cowtan and University of York\n";
  std::cout << "All rights reserved. Please reference:\n";
  std::cout << " Cowtan K. (2009) 'Sloop loop-building software'.\n";

  // defaults
  clipper::String title;
  clipper::String ippdb_wrk = "NONE";
  clipper::String ipmap_wrk = "NONE";
  clipper::String ipmtz_wrk = "NONE";
  clipper::String ipcol_wrk_fo = "NONE";
  clipper::String ipcol_wrk_hl = "NONE";
  clipper::String ipcol_wrk_fc = "NONE";
  clipper::String ipcol_wrk_fr = "NONE";
  clipper::String oppdb = "sloop.pdb";
  clipper::String opfrg = "NONE";
  clipper::String opmap = "NONE";
  clipper::String prtdb = "NONE";
  double res_in = 1.0;         // Resolution limit
  double clashr = 2.0;         // clash radius
  double clashp = 1.0;         // clash penalty
  int verbose = 0;

  // command input
  CCP4CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if        ( args[arg] == "-title" ) {
      if ( ++arg < args.size() ) title = args[arg];
    } else if ( args[arg] == "-protein-db" ) {
      if ( ++arg < args.size() ) prtdb = args[arg];
    } else if ( args[arg] == "-pdbin" ) {
      if ( ++arg < args.size() ) ippdb_wrk = args[arg];
    } else if ( args[arg] == "-mapin" ) {
      if ( ++arg < args.size() ) ipmap_wrk = args[arg];
    } else if ( args[arg] == "-mtzin" ) {
      if ( ++arg < args.size() ) ipmtz_wrk = args[arg];
    } else if ( args[arg] == "-pdbout" ) {
      if ( ++arg < args.size() ) oppdb = args[arg];
    } else if ( args[arg] == "-pdbout-loops" ) {
      if ( ++arg < args.size() ) opfrg = args[arg];
    } else if ( args[arg] == "-mapout" ) {
      if ( ++arg < args.size() ) opmap  = args[arg];
    } else if ( args[arg] == "-colin-fo" ) {
      if ( ++arg < args.size() ) ipcol_wrk_fo = args[arg];
    } else if ( args[arg] == "-colin-hl" ) {
      if ( ++arg < args.size() ) ipcol_wrk_hl = args[arg];
    } else if ( args[arg] == "-colin-fc" ) {
      if ( ++arg < args.size() ) ipcol_wrk_fc = args[arg];
    } else if ( args[arg] == "-colin-free" ) {
      if ( ++arg < args.size() ) ipcol_wrk_fr = args[arg];
    } else if ( args[arg] == "-resolution" ) {
      if ( ++arg < args.size() ) res_in = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-clash-radius" ) {
      if ( ++arg < args.size() ) clashr = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-clash-penalty" ) {
      if ( ++arg < args.size() ) clashp = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-verbose" ) {
      if ( ++arg < args.size() ) verbose = clipper::String(args[arg]).i();
    } else {
      std::cout << "\nUnrecognized:\t" << args[arg] << std::endl;
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "\nUsage: csloop\n\t-pdbin <filename>\t\tCOMPULSORY\n\t-mapin <filename>\n\t-mtzin <filename>\n\t-pdbout <filename>\n\t-pdbout-loops <filename>\n\t-colin-fo <colpath>\n\t-colin-hl <colpath>\n\t-colin-fc <colpath>\n\t-colin-free <colpath>\n\t-resolution <resolution/A>\n\t-clash-radius <radius/A>\n\t-clash-penalty <penalty>\n\t-protein-db <filename>\nAn input pdb and mtz file or map are required. Missing residue ranges from correctly numbered chains will be inserted if possible. The model will be written to the output pdb file. Optionally, the loops can also be written to output file, one per chain.\n";
    clipper::Message::message(clipper::Message_fatal("Invalid option."));
  }

  using clipper::data32::Compute_phifom_from_abcd;
  using clipper::data32::Compute_fphi_from_fsigf_phifom;
  using clipper::data32::F_sigF;
  using clipper::data32::F_phi;
  using clipper::data32::Phi_fom;
  using clipper::data32::ABCD;
  using clipper::data32::Flag;

  // numbers to output
  clipper::Resolution resol;
  clipper::CCP4MTZfile mtzfile;
  mtzfile.set_column_label_mode( clipper::CCP4MTZfile::Legacy );
 
  // Get resolution for calculation
  mtzfile.open_read( ipmtz_wrk );
  double res_wrk = clipper::Util::max( mtzfile.resolution().limit(), res_in );
  mtzfile.close_read();
  resol = clipper::Resolution( res_wrk );

  // get work map
  clipper::Xmap<float> xwrk;

  // fill it
  if ( ipmtz_wrk != "NONE" ) {
    // Get work reflection data
    clipper::HKL_info hkls_wrk;
    mtzfile.open_read( ipmtz_wrk );
    hkls_wrk.init( mtzfile.spacegroup(), mtzfile.cell(), resol, true );
    clipper::HKL_data<F_sigF>  wrk_f( hkls_wrk );
    clipper::HKL_data<ABCD>    wrk_hl( hkls_wrk );
    clipper::HKL_data<Phi_fom> wrk_pw( hkls_wrk );
    clipper::HKL_data<F_phi>   wrk_fp( hkls_wrk );
    clipper::HKL_data<Flag>    flag( hkls_wrk );
    if ( ipcol_wrk_fo != "NONE" ) mtzfile.import_hkl_data(wrk_f,ipcol_wrk_fo);
    if ( ipcol_wrk_hl != "NONE" ) mtzfile.import_hkl_data(wrk_hl,ipcol_wrk_hl);
    if ( ipcol_wrk_fc != "NONE" ) mtzfile.import_hkl_data(wrk_fp,ipcol_wrk_fc);
    if ( ipcol_wrk_fr != "NONE" ) mtzfile.import_hkl_data(flag,ipcol_wrk_fr);
    mtzfile.close_read();

    // apply free flag
    clipper::HKL_data<F_sigF> wrk_fwrk = wrk_f;
    wrk_fwrk.mask( flag != 0 );

    // work map
    wrk_pw.compute( wrk_hl, Compute_phifom_from_abcd() );
    if ( ipcol_wrk_fc == "NONE" )
      wrk_fp.compute( wrk_fwrk, wrk_pw, Compute_fphi_from_fsigf_phifom() );
    clipper::Grid_sampling grid( hkls_wrk.spacegroup(),
				 hkls_wrk.cell(), hkls_wrk.resolution() );
    xwrk.init( hkls_wrk.spacegroup(), hkls_wrk.cell(), grid );
    xwrk.fft_from( wrk_fp );
  } else {
    clipper::CCP4MAPfile mapfile;
    mapfile.open_read( ipmap_wrk );
    mapfile.import_xmap( xwrk );
    mapfile.close_read();
  }

  // optionlly write work map
  if ( opmap != "NONE" ) {
    clipper::CCP4MAPfile mapfile;
    mapfile.open_write( opmap );
    mapfile.export_xmap( xwrk );
    mapfile.close_write();
  }

  // Get work model
  clipper::MiniMol mol_wrk( xwrk.spacegroup(), xwrk.cell() );
  {
    clipper::MiniMol mol_tmp;
    clipper::MMDBfile mmdb;
    mmdb.read_file( ippdb_wrk );
    mmdb.import_minimol( mol_tmp );
    mol_wrk.copy( mol_tmp, clipper::MM::COPY_MPC );
  }

  // fetch db
  if ( prtdb == "NONE" ) {
    const char* clibdptr = getenv( "CLIBD" );
    if ( clibdptr != NULL ) prtdb = clipper::String( clibdptr ) + "/prot500.db";
  }
  std::cout << std::endl << "Reading Top500 protein DB" << std::endl
	    << " See http://kinemage.biochem.duke.edu/databases/top500.php"
	    << std::endl << std::endl;
  const ProteinDB::ChainDB chaindb( prtdb );
  if ( chaindb.size() == 0 ) clipper::Message::message( clipper::Message_fatal( "Could not find chain DB: check or specify --protein-db\n" ) );
  if ( verbose > 0 ) chaindb.debug();

  // Do loop build
  {
    clipper::Spacegroup    spgr = xwrk.spacegroup();
    clipper::Cell          cell = xwrk.cell();
    clipper::Grid_sampling grid = xwrk.grid_sampling();
    clipper::MiniMol mol_new( spgr, cell );

    // density scoring class
    ProteinDB::ScoreDensity score_rho( xwrk, 0.0, 2.0 );

    // clash scoring class
    std::vector<clipper::Coord_orth> nnbcoords;
    for ( int c1 = 0; c1 < mol_wrk.size(); c1++ )
      for ( int r1 = 0; r1 < mol_wrk[c1].size(); r1++ )
	for ( int a1 = 0; a1 < mol_wrk[c1][r1].size(); a1++ ) {
	  clipper::String t = mol_wrk[c1][r1][a1].id();
	  if ( t == " N  " || t == " CA " || t == " C  " )
	    nnbcoords.push_back( mol_wrk[c1][r1][a1].coord_orth() );
	}
    ProteinDB::ScoreClashes score_cls( nnbcoords, spgr, cell, clashr );

    // search for chain breaks
    for ( int c = 0; c < mol_wrk.size(); c++ ) {
      clipper::MiniMol mol_frg( spgr, cell );
      clipper::MPolymer mp;
      mp.set_id( mol_wrk[c].id() );
      for ( int r = 0; r < mol_wrk[c].size(); r++ ) {
	clipper::MMonomer mm = mol_wrk[c][r];
	if ( r >= 2 && r < mol_wrk[c].size()-2 ) {
	  int seq1 = mol_wrk[c][r-1].seqnum();
	  int seq2 = mol_wrk[c][r+0].seqnum();
	  ProteinDB::Residue r1( mol_wrk[c][r-2] );
	  ProteinDB::Residue r2( mol_wrk[c][r-1] );
	  ProteinDB::Residue r3( mol_wrk[c][r+0] );
	  ProteinDB::Residue r4( mol_wrk[c][r+1] );
	  ProteinDB::Residue::FLAG normal = ProteinDB::Residue::NORMAL;
	  if ( seq2 > seq1+1 &&
	       r1.flag() == normal && r2.flag() == normal &&
	       r3.flag() == normal && r4.flag() == normal ) {
	    // break found
	    std::cout << "Chain " << mol_wrk[c].id() << ", "
		      << "Residues " << seq1 << "-" << seq2 << std::endl;
	    // construct search fragment
	    ProteinDB::Chain frag;
	    frag.add_residue( r1 );
	    frag.add_residue( r2 );
	    for ( int i = seq1+1; i < seq2; i++ )
	      frag.add_residue( ProteinDB::Residue() );
	    frag.add_residue( r3 );
	    frag.add_residue( r4 );

	    // /*
	    // search for fragment
	    ProteinDB::ChainDB fragdb( frag );
	    std::vector<ProteinDB::Chain> frags =
	      chaindb.match_fragment( fragdb, 50 );

	    if ( frags.size() > 0 ) {

	      // exclude the initial fragment environment from clash scoring
	      score_cls.set_exclude( fragdb );

	      // score vs density and clashes (omitting first and last residues)
	      std::vector<std::pair<double,int> > fragscore( frags.size() );
	      for ( int f = 0; f < frags.size(); f++ ) {
		// density score
		double scr_rho = score_rho.score( frags[f] );
		// clash score
		double scr_cls = score_cls.score( frags[f] ); 
		// store
		double s = scr_rho + clashp * scr_cls;
		//std::cout << scr_rho << " " << scr_cls << std::endl;
		fragscore[f] = std::pair<double,int>( -s, f );
	      }

	      // sort
	      std::sort( fragscore.begin(), fragscore.end() );

	      // update the C and N atom coords either end of the fragment
	      const int f = fragscore[0].second;
	      const int r1 = 1;
	      const int r2 = frags[f].size()-2;
	      clipper::MMonomer& mm1 = mp[mp.size()-1];
	      const int cc1 = mm1.lookup( " C  ", clipper::MM::ANY );
	      if ( cc1 >= 0 ) mm1[cc1].set_coord_orth( frags[f][r1].coord_c() );
	      clipper::MMonomer& mm2 = mm;
	      const int cn2 = mm2.lookup( " N  ", clipper::MM::ANY );
	      if ( cn2 >= 0 ) mm2[cn2].set_coord_orth( frags[f][r2].coord_n() );
	      // and add best fragment to the model
	      for ( int r = 2; r < frags[f].size()-2; r++ ) {
		clipper::MMonomer mm = frags[f][r].mmonomer();
		mm.set_seqnum( seq1+r-1 );
		mm.set_type( "UNK" );
		mp.insert( mm );
	      }

	      // store fragments for output
	      clipper::String ids = "123456789";
	      for ( int i = 0; i < fragscore.size(); i++ ) {
		if ( i < ids.size() ) {
		  clipper::MPolymer fp;
		  fp.set_id( ids.substr(i,1) );
		  int f = fragscore[i].second;
		  for ( int r = 0; r < frags[f].size(); r++ ) {
		    clipper::MMonomer fm = frags[f][r].mmonomer();
		    fm.set_seqnum( seq1+r );
		    fm.set_type( "UNK" );
		    fp.insert( fm );
		  }
		  mol_frg.insert( fp );
		  std::cout << " Fragment: " << fp.id() << "  " 
			    << "LLK score: " << -fragscore[i].first
			    << std::endl;
		}
	      }
	    }
	    // */

	    /*
	    ProteinDB::ProteinDBSearch pdbsrch( prtdb );
	    std::vector<ProteinDB::Chain> frags =
	      pdbsrch.search( frag, 50, score_rho, score_cls );
	    std::vector<double> scores = pdbsrch.scores();
	    if ( frags.size() > 0 ) {
	      // update the C and N atom coords either end of the fragment
	      const int r1 = 1;
	      const int r2 = frags[0].size()-2;
	      clipper::MMonomer& mm1 = mp[mp.size()-1];
	      const int cc1 = mm1.lookup( " C  ", clipper::MM::ANY );
	      if ( cc1 >= 0 ) mm1[cc1].set_coord_orth( frags[0][r1].coord_c() );
	      clipper::MMonomer& mm2 = mm;
	      const int cn2 = mm2.lookup( " N  ", clipper::MM::ANY );
	      if ( cn2 >= 0 ) mm2[cn2].set_coord_orth( frags[0][r2].coord_n() );
	      // and add best fragment to the model
	      for ( int r = 2; r < frags[0].size()-2; r++ ) {
		clipper::MMonomer mm = frags[0][r].mmonomer();
		mm.set_seqnum( seq1+r-1 );
		mm.set_type( "UNK" );
		mp.insert( mm );
	      }

	      // store fragments for output
	      clipper::String ids = "123456789";
	      for ( int f = 0; f < frags.size(); f++ ) {
		if ( f < ids.size() ) {
		  clipper::MPolymer fp;
		  fp.set_id( ids.substr(f,1) );
		  for ( int r = 0; r < frags[f].size(); r++ ) {
		    clipper::MMonomer fm = frags[f][r].mmonomer();
		    fm.set_seqnum( seq1+r );
		    fm.set_type( "UNK" );
		    fp.insert( fm );
		  }
		  mol_frg.insert( fp );
		  std::cout << " Fragment: " << fp.id() << "  " 
			    << "LLK score: " << scores[f]
			    << std::endl;
		}
	      }	      
	    }
	    */

	  }
	}
	mp.insert( mm );
      }
      mol_new.insert( mp );

      // write fragment molecule
      if ( opfrg != "NONE" && mol_frg.size() > 0 ) {
	int i  = opfrg.rfind( "." );
	if ( i != std::string::npos ) opfrg = opfrg.substr( 0, i );
	clipper::MMDBfile mmdbfrag;
	mmdbfrag.export_minimol( mol_frg );
	mmdbfrag.write_file( opfrg+"."+mol_wrk[c].id()+".pdb" );
      }
    }

    // write answers
    clipper::MMDBfile mmdb;
    mmdb.export_minimol( mol_new );
    mmdb.write_file( oppdb );
  }
}
