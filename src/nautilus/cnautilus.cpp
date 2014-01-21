// Clipper ssfind
/* Copyright 2008 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-minimol.h>
#include <algorithm>

#include "nucleicacid_db.h"
#include "nautilus-tools.h"
#include "nautilus-ss-find.h"
#include "nautilus-target.h"
#include "nautilus-join.h"
#include "nautilus-sequence.h"
#include "nautilus-rebuild-bases.h"
#include "nautilus-tidy.h"
#include "nautilus-util.h"


int main( int argc, char** argv )
{
  CCP4Program prog( "cnautilus", "0.4", "$Date: 2014/01/10" );
  prog.set_termination_message( "Failed" );

  // defaults
  clipper::String title;
  clipper::String ipmtz = "NONE";
  clipper::String ipcol_fo = "NONE";
  clipper::String ipcol_hl = "NONE";
  clipper::String ipcol_pw = "NONE";
  clipper::String ipcol_fc = "NONE";
  clipper::String ipcol_fr = "NONE";
  clipper::String ipseq = "NONE";
  clipper::String ippdb = "NONE";
  clipper::String ippdb_ref = "NONE";
  clipper::String oppdb = "nautilus.pdb";
  clipper::String opmap = "NONE";
  int ncyc = 3;
  bool doanis = false;
  int nhit = 100;
  double res_in = 2.0;         // Resolution limit
  double srchst = 18.0;        // Search angle step
  int verbose = 5;

  // command input
  CCP4CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if        ( args[arg] == "-title" ) {
      if ( ++arg < args.size() ) title = args[arg];
    } else if ( args[arg] == "-pdbin-ref" ) {
      if ( ++arg < args.size() ) ippdb_ref = args[arg];
    } else if ( args[arg] == "-mtzin" ) {
      if ( ++arg < args.size() ) ipmtz = args[arg];
    } else if ( args[arg] == "-seqin" ) {
      if ( ++arg < args.size() ) ipseq = args[arg];
    } else if ( args[arg] == "-pdbin" ) {
      if ( ++arg < args.size() ) ippdb = args[arg];
    } else if ( args[arg] == "-pdbout" ) {
      if ( ++arg < args.size() ) oppdb = args[arg];
    } else if ( args[arg] == "-mapout" ) {
      if ( ++arg < args.size() ) opmap  = args[arg];
    } else if ( args[arg] == "-colin-fo" ) {
      if ( ++arg < args.size() ) ipcol_fo = args[arg];
    } else if ( args[arg] == "-colin-hl" ) {
      if ( ++arg < args.size() ) ipcol_hl = args[arg];
    } else if ( args[arg] == "-colin-phifom" ) {
      if ( ++arg < args.size() ) ipcol_pw = args[arg];
    } else if ( args[arg] == "-colin-fc" ) {
      if ( ++arg < args.size() ) ipcol_fc = args[arg];
    } else if ( args[arg] == "-colin-free" ) {
      if ( ++arg < args.size() ) ipcol_fr = args[arg];
    } else if ( args[arg] == "-cycles" ) {
      if ( ++arg < args.size() ) ncyc  = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-anisotropy-correction" ) {
      doanis = true;
    } else if ( args[arg] == "-fragments" ) {
      if ( ++arg < args.size() ) nhit = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-resolution" ) {
      if ( ++arg < args.size() ) res_in = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-search-step" ) {
      if ( ++arg < args.size() ) srchst = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-verbose" ) {
      if ( ++arg < args.size() ) verbose = clipper::String(args[arg]).i();
    } else {
      std::cout << "\nUnrecognized:\t" << args[arg] << std::endl;
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "\nUsage: cnautilus\n\t-mtzin <filename>\t\tCOMPULSORY\n\t-seqin <filename>\n\t-pdbin <filename>\n\t-pdbout <filename>\n\t-colin-fo <colpath>\n\t-colin-hl <colpath> or -colin-phifom <colpath>\n\t-colin-fc <colpath>\n\t-colin-free <colpath>\n\t-cycles <number>\n\t-anisotropy-correction <number>\n\t-fragments <number>\n\t-resolution <resolution/A>\n\t-pdbin-ref <filename>\n.\n";
    return 1;
  }

  // check data present
  if ( ipcol_fc == "NONE" && ipcol_fo == "NONE" )
    { std::cerr << "No F's provided." << std::endl; return 1; }
  if ( ipcol_fc == "NONE" && ipcol_hl == "NONE" && ipcol_pw == "NONE" )
    { std::cerr << "No phases provided." << std::endl; return 1; }

  // other initialisations
  typedef clipper::HKL_data_base::HKL_reference_index HRI;
  using clipper::data32::Compute_fphi_from_fsigf_phifom;
  using clipper::data32::Compute_scale_u_aniso_fphi;
  clipper::Resolution resol;
  clipper::CCP4MTZfile mtzfile;
  mtzfile.set_column_label_mode( clipper::CCP4MTZfile::Legacy );
  const int mmdbflags = MMDBF_IgnoreBlankLines | MMDBF_IgnoreDuplSeqNum | MMDBF_IgnoreNonCoorPDBErrors | MMDBF_IgnoreRemarks;
  typedef clipper::HKL_data_base::HKL_reference_index HRI;

  // Get work reflection data
  clipper::HKL_info hkls;
  mtzfile.open_read( ipmtz );
  double res = clipper::Util::max( mtzfile.resolution().limit(), res_in );
  resol = clipper::Resolution( res );
  hkls.init( mtzfile.spacegroup(), mtzfile.cell(), resol, true );
  clipper::HKL_data<clipper::data32::F_sigF>  wrk_f ( hkls );
  clipper::HKL_data<clipper::data32::ABCD>    wrk_hl( hkls );
  clipper::HKL_data<clipper::data32::Phi_fom> wrk_pw( hkls );
  clipper::HKL_data<clipper::data32::F_phi>   fphi( hkls );
  clipper::HKL_data<clipper::data32::Flag>    flag( hkls );
  if ( ipcol_fo != "NONE" ) mtzfile.import_hkl_data( wrk_f ,ipcol_fo );
  if ( ipcol_hl != "NONE" ) mtzfile.import_hkl_data( wrk_hl,ipcol_hl );
  if ( ipcol_pw != "NONE" ) mtzfile.import_hkl_data( wrk_pw,ipcol_pw );
  if ( ipcol_fc != "NONE" ) mtzfile.import_hkl_data( fphi,  ipcol_fc );
  if ( ipcol_fr != "NONE" ) mtzfile.import_hkl_data( flag,  ipcol_fr );
  mtzfile.close_read();

  // do anisotropy correction
  clipper::U_aniso_orth uaniso( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
  if ( doanis ) {
    if ( ipcol_fo == "NONE" ) for ( HRI ih = hkls.first(); !ih.last(); ih.next() ) wrk_f[ih].f() = fphi[ih].f();
    // scale obs data
    typedef clipper::SFscale_aniso<float> SFscale;
    SFscale sfscl( 3.0, SFscale::SHARPEN );
    sfscl( wrk_f );
    uaniso = sfscl.u_aniso_orth( SFscale::F );
    // scale map coeffs
    Compute_scale_u_aniso_fphi compute_aniso( 1.0, -uaniso );
    if ( ipcol_fc != "NONE" ) fphi.compute( fphi, compute_aniso );    
    // output
    std::cout << std::endl << "Applying anisotropy correction:"
	      << std::endl << uaniso.format() << std::endl << std::endl;
  }

  // apply free flag
  clipper::HKL_data<clipper::data32::F_sigF> wrk_f1 = wrk_f;
  //wrk_f1.mask( flag != 0 );
  for ( HRI ih = hkls.first(); !ih.last(); ih.next() ) if ( flag[ih].flag() == 0 ) wrk_f1[ih] = clipper::data32::F_sigF();  //ugly hack for broken SGI compilers

  // Get reference model
  if ( ippdb_ref == "NONE" ) NautilusUtil::set_reference( ippdb_ref );
  NucleicAcidTargets natools;
  natools.add_pdb( ippdb_ref );
  NucleicAcidTools tools;

  // Get sequence
  clipper::MMoleculeSequence seq_wrk;
  if ( ipseq != "NONE" ) {
    clipper::SEQfile seqf_wrk;
    seqf_wrk.read_file( ipseq );
    seqf_wrk.import_molecule_sequence( seq_wrk );
  }

  // Get model
  clipper::MiniMol mol_wrk( hkls.spacegroup(), hkls.cell() );
  if ( ippdb != "NONE" ) {
    clipper::MiniMol mol_tmp;
    clipper::MMDBfile mmdb;
    mmdb.SetFlag( mmdbflags );
    mmdb.read_file( ippdb );
    mmdb.import_minimol( mol_tmp );
    std::cout << mol_tmp.spacegroup().symbol_hm() << " " << mol_tmp.cell().format() << " " << mol_tmp.atom_list().size() << std::endl;
    for ( int c = 0; c < mol_tmp.size(); c++ ) mol_wrk.insert( mol_tmp[c] );
  }
  
  // work map
  if ( ipcol_hl == "NONE" )
    wrk_hl.compute( wrk_pw, clipper::data32::Compute_abcd_from_phifom() );
  if ( ipcol_pw == "NONE" )
    wrk_pw.compute( wrk_hl, clipper::data32::Compute_phifom_from_abcd() );
  if ( ipcol_fc == "NONE" )
    fphi.compute( wrk_f1, wrk_pw, Compute_fphi_from_fsigf_phifom() );
  clipper::Spacegroup cspg = hkls.spacegroup();
  clipper::Cell       cxtl = hkls.cell();
  clipper::Grid_sampling grid( cspg, cxtl, hkls.resolution() );
  clipper::Xmap<float>   xwrk( cspg, cxtl, grid );
  xwrk.fft_from( fphi );

  // output some statistics
  std::cout << std::endl;
  std::cout << " Spgr " << hkls.spacegroup().symbol_xhm() << std::endl;
  std::cout << hkls.cell().format() << std::endl;
  std::cout << " Nref " << hkls.num_reflections() << " " << fphi.num_obs() << std::endl;
  double smax = 0.0;
  for ( HRI ih = fphi.first(); !ih.last(); ih.next() )
    if ( !fphi[ih].missing() )
      if ( fphi[ih].f() > 0.0 )
	smax = std::max( smax, double(hkls.invresolsq(ih.index())) );
  std::cout << " Reso " << hkls.resolution().limit() << " " << 1.0/std::max(sqrt(smax),1.0e-3) << std::endl;
  if ( ipcol_fo != "NONE" ) {
    double sf(0.0), sw(0.0);
    for ( HRI ih = fphi.first(); !ih.last(); ih.next() )
      if ( !wrk_f1[ih].missing() && !fphi[ih].missing() ) {
	sf += wrk_f1[ih].f();
	sw += fphi[ih].f();
      }
    std::cout << " Fw/Fo " << sw/sf << std::endl;
  }
  if ( smax == 0.0 )
    { std::cerr << "No density provided." << std::endl; return 1; }

  // store copy of input model
  clipper::MiniMol mol_wrk_in = mol_wrk;

  // map stats
  natools.init_stats( xwrk );
  NautilusLog log;
  std::cout << std::endl;

  for ( int cyc = 0; cyc < ncyc; cyc++ ) {
    std::cout << "Internal cycle " << clipper::String( cyc+1, 3 ) << std::endl;

    // adjust labels and label non-NA chains to keep
    mol_wrk = NucleicAcidTools::flag_chains( mol_wrk );

    // find chains
    mol_wrk = natools.find( xwrk, mol_wrk, nhit/2, nhit/2, srchst );
    log.log( "FIND", mol_wrk, verbose >= 5 );

    // grow chains
    mol_wrk = natools.grow( xwrk, mol_wrk, 25, 0.001 );
    log.log( "GROW", mol_wrk, verbose >= 5 );

    // join
    NucleicAcidJoin na_join;
    mol_wrk = na_join.join( mol_wrk );
    log.log( "JOIN", mol_wrk, verbose >= 5 );
    //for ( int c = 0; c < mol_wrk.size(); c++ ) { for ( int r = 0; r < mol_wrk[c].size(); r++ ) std::cout << mol_wrk[c][r].type().trim(); std::cout << std::endl; }

    // link
    mol_wrk = natools.link( xwrk, mol_wrk );
    log.log( "LINK", mol_wrk, verbose >= 5 );
    //for ( int c = 0; c < mol_wrk.size(); c++ ) { for ( int r = 0; r < mol_wrk[c].size(); r++ ) std::cout << mol_wrk[c][r].type().trim(); std::cout << std::endl; }

    // prune
    mol_wrk = natools.prune( mol_wrk );
    log.log( "PRUNE", mol_wrk, verbose >= 5 );
    //for ( int c = 0; c < mol_wrk.size(); c++ ) { for ( int r = 0; r < mol_wrk[c].size(); r++ ) std::cout << mol_wrk[c][r].type().trim(); std::cout << std::endl; }

    // rebuild
    mol_wrk = natools.rebuild_chain( xwrk, mol_wrk );
    log.log( "CHAIN", mol_wrk, verbose >= 5 );
    //for ( int c = 0; c < mol_wrk.size(); c++ ) { for ( int r = 0; r < mol_wrk[c].size(); r++ ) std::cout << mol_wrk[c][r].type().trim(); std::cout << std::endl; }

    // sequence
    NucleicAcidSequence na_seqnc;
    mol_wrk = na_seqnc.sequence( xwrk, mol_wrk, seq_wrk );
    log.log( "SEQNC", mol_wrk, verbose >= 5 );
    //for ( int c = 0; c < mol_wrk.size(); c++ ) { for ( int r = 0; r < mol_wrk[c].size(); r++ ) std::cout << mol_wrk[c][r].type().trim(); std::cout << std::endl; }

    // rebuild
    NucleicAcidRebuildBases na_bases;
    mol_wrk = na_bases.rebuild_bases( xwrk, mol_wrk );
    log.log( "BASES", mol_wrk, verbose >= 5 );
    for ( int c = 0; c < mol_wrk.size(); c++ ) { for ( int r = 0; r < mol_wrk[c].size(); r++ ) std::cout << mol_wrk[c][r].type().trim(); std::cout << std::endl; }

    prog.summary_beg();
    clipper::String msg = log.log_info( mol_wrk );
    std::cout << "Internal cycle " << clipper::String( cyc+1, 3 ) << std::endl << msg << std::endl;
    prog.summary_end();
  }

  // move to match input model
  if ( mol_wrk_in.size() > 0 )
    NucleicAcidTools::symm_match( mol_wrk, mol_wrk_in );

  // output
  const clipper::String basetypes = "ACGTU";
  clipper::MiniMol mol_new( xwrk.spacegroup(), xwrk.cell() );
  const clipper::String chainid1 = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  const clipper::String chainid2 = "abcdefghijklmnopqrstuvwxyz";
  for ( int i = 0; i < mol_wrk.size(); i++ ) {
    clipper::String id; int roff;
    if ( i < 26 ) {
      id = chainid1.substr( i, 1 ); roff = 0;
    } else {
      id = chainid2.substr( (i-26)/100, 1 ); roff = 10*(i%100);
    }
    clipper::MPolymer mpx = mol_wrk[i];
    mpx.set_id( id );
    if ( !mpx.exists_property( "NON-NA" ) ) {
      for ( int r = 0; r < mpx.size(); r++ ) {
	const clipper::String type = mpx[r].type().trim()+" ";
	const char ctype = type[0];
	int t = NucleicAcidTools::base_index( ctype );
	if ( t >= 0 ) mpx[r].set_type( "  "+basetypes.substr(t,1) );
	else          mpx[r].set_type( "  U" );
	mpx[r].set_seqnum( roff+r+1 );
      }
    }
    mol_new.insert( mpx );
  }
  ModelTidy::chain_renumber( mol_new, seq_wrk );
  clipper::MMDBfile pdbfile;
  pdbfile.export_minimol( mol_new );
  pdbfile.write_file( oppdb );

  log.profile();
  prog.set_termination_message( "Normal termination" );
}
