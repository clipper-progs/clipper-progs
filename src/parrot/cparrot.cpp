// Clipper parrot
/* Copyright 2008 Kevin Cowtan & University of York all rights reserved */


#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>

#include "simulate-lib.h"
#include "parrot-lib.h"
#include "parrot-ncsfind.h"
#include "parrot-ncsaver.h"
extern "C" {
#include <stdio.h>
#include <stdlib.h>
}


int main( int argc, char** argv )
{
  CCP4Program prog( "cparrot", "1.0.4", "$Date: 2015/05/20" );
  prog.set_termination_message( "Failed" );

  std::cout << std::endl << "Copyright 2008-2010 Kevin Cowtan and University of York." << std::endl << std::endl;
  prog.summary_beg();
  std::cout << "$TEXT:Reference: $$ Please reference $$" << std::endl << std::endl << " 'Recent developments in classical density modification.'" << std::endl << " Cowtan K. (2010) Acta Cryst. D66, 470-478." << std::endl << std::endl << "$$";
  prog.summary_end();

  // defaults
  clipper::String title;
  clipper::String ipmtz_ref = "NONE";
  clipper::String ipmtz_wrk = "NONE";
  clipper::String ippdb_ref = "NONE";
  clipper::String ippdb_wrk = "NONE";
  clipper::String ippdb_wrk_ha = "NONE";
  clipper::String ippdb_wrk_mr = "NONE";
  clipper::String ipseq_wrk = "NONE";
  clipper::String ipcol_ref_fo = "/*/*/FP";
  clipper::String ipcol_ref_hl = "/*/*/FC";
  clipper::String ipcol_wrk_fo = "NONE";
  clipper::String ipcol_wrk_hl = "NONE";
  clipper::String ipcol_wrk_pw = "NONE";
  clipper::String ipcol_wrk_fc = "NONE";
  clipper::String ipcol_wrk_fr = "NONE";
  clipper::String opfile = "parrot.mtz";
  clipper::String opcol = "parrot";
  clipper::String opcol_hl = "NONE";
  clipper::String opcol_fc = "NONE";
  clipper::String opncsm = "NONE";
  clipper::String opxml = "NONE";
  double res_in = 1.0;         // Resolution limit
  double solrad = 0.0;
  double solc   = 0.0;         // solvent content
  double usharp = 0.0;
  double ncs_radius = 6.00;    // Radius of NCS local correlation
  double ncs_asufrc = 0.10;    // min fraction of ASU for NCS mask
  std::vector<Local_rtop> nxops;  // ncs operators
  bool dosolv = false;
  bool dohist = false;
  bool doncsa = false;
  bool domask = false;
  bool dorice = false;
  bool doanis = false;
  bool dodump = false;
  int nncs = 1;
  int ncycles = 3;
  int n_refln = 2000;
  int n_param = 10;
  int verbose = 0;

  // command input
  CCP4CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    std::string key = args[arg];
    if        ( key == "-title" ) {
      if ( ++arg < args.size() ) title = args[arg];
    } else if ( key == "-mtzin-ref" ) {
      if ( ++arg < args.size() ) ipmtz_ref = args[arg];
    } else if ( key == "-pdbin-ref" ) {
      if ( ++arg < args.size() ) ippdb_ref = args[arg];
    } else if ( key == "-mtzin"        || key == "-mtzin-wrk" ) {
      if ( ++arg < args.size() ) ipmtz_wrk = args[arg];
    } else if ( key == "-seqin"        || key == "-seqin-wrk" ) {
      if ( ++arg < args.size() ) ipseq_wrk = args[arg];
    } else if ( key == "-pdbin"        || key == "-pdbin-wrk" ) {
      if ( ++arg < args.size() ) ippdb_wrk = args[arg];
    } else if ( key == "-pdbin-ha"     || key == "-pdbin-wrk-ha" ) {
      if ( ++arg < args.size() ) ippdb_wrk_ha  = args[arg];
    } else if ( key == "-pdbin-mr"     || key == "-pdbin-wrk-mr" ) {
      if ( ++arg < args.size() ) ippdb_wrk_mr  = args[arg];
    } else if ( key == "-colin-ref-fo" ) {
      if ( ++arg < args.size() ) ipcol_ref_fo = args[arg];
    } else if ( key == "-colin-ref-hl" ) {
      if ( ++arg < args.size() ) ipcol_ref_hl = args[arg];
    } else if ( key == "-colin-fo"     || key == "-colin-wrk-fo" ) {
      if ( ++arg < args.size() ) ipcol_wrk_fo = args[arg];
    } else if ( key == "-colin-hl"     || key == "-colin-wrk-hl" ) {
      if ( ++arg < args.size() ) ipcol_wrk_hl = args[arg];
    } else if ( key == "-colin-phifom" || key == "-colin-wrk-phifom" ) {
      if ( ++arg < args.size() ) ipcol_wrk_pw = args[arg];
    } else if ( key == "-colin-fc"     || key == "-colin-wrk-fc" ) {
      if ( ++arg < args.size() ) ipcol_wrk_fc = args[arg];
    } else if ( key == "-colin-free"   || key == "-colin-wrk-free" ) {
      if ( ++arg < args.size() ) ipcol_wrk_fr = args[arg];
    } else if ( key == "-mtzout" ) {
      if ( ++arg < args.size() ) opfile = args[arg];
    } else if ( key == "-colout" ) {
      if ( ++arg < args.size() ) opcol = args[arg];
    } else if ( key == "-colout-hl" ) {
      if ( ++arg < args.size() ) opcol_hl = args[arg];
    } else if ( key == "-colout-fc" ) {
      if ( ++arg < args.size() ) opcol_fc = args[arg];
    } else if ( key == "-mapout-ncs" ) {
      if ( ++arg < args.size() ) opncsm = args[arg];
    } else if ( key == "-xmlout" ) {
      if ( ++arg < args.size() ) opxml  = args[arg];
    } else if ( key == "-solvent-flatten" ) {
      dosolv = true;
    } else if ( key == "-histogram-match" ) {
      dohist = true;
    } else if ( key == "-ncs-average" ) {
      doncsa = true;
    } else if ( key == "-force-solvent-mask-calculation" ) {
      domask = true;
    } else if ( key == "-rice-probability" ) {
      dorice = true;
    } else if ( key == "-modify-map-and-terminate" ) {
      dodump = true;
    } else if ( key == "-resolution" ) {
      if ( ++arg < args.size() ) res_in = clipper::String(args[arg]).f();
    } else if ( key == "-cycles" ) {
      if ( ++arg < args.size() ) ncycles = clipper::String(args[arg]).i();
    } else if ( key == "-anisotropy-correction" ) {
      doanis = true;
    } else if ( key == "-solvent-content" ) {
      if ( ++arg < args.size() ) solc   = clipper::String(args[arg]).f();
    } else if ( key == "-solvent-mask-filter-radius" ) {
      if ( ++arg < args.size() ) solrad = clipper::String(args[arg]).f();
    } else if ( key == "-ncs-mask-filter-radius" ) {
      if ( ++arg < args.size() ) ncs_radius = clipper::String(args[arg]).f();
    } else if ( key == "-ncs-asu-fraction" ) {
      if ( ++arg < args.size() ) ncs_asufrc = clipper::String(args[arg]).f();
    } else if ( key == "-ncs-operator" ) {
      if ( ++arg < args.size() ) {
        std::vector<clipper::String> op=clipper::String(args[arg]).split(", ");
        if ( op.size() == 9 ) {
          clipper::Euler_ccp4 eul( clipper::Util::d2rad(op[0].f64()),
                                   clipper::Util::d2rad(op[1].f64()),
                                   clipper::Util::d2rad(op[2].f64()) );
          clipper::Coord_orth src( op[3].f64(), op[4].f64(), op[5].f64() );
          clipper::Coord_orth tgt( op[6].f64(), op[7].f64(), op[8].f64() );
          clipper::Rotation rot( eul );
          nxops.push_back( Local_rtop( rot, src, tgt ) );
        } else {
          std::cout << "\nInvalid ncs operator:\t" << args[arg] << "\n";
          args.clear();
        }
      }
    } else if ( key == "-sharpen" ) {
      if ( ++arg < args.size() ) usharp = clipper::String(args[arg]).f();
    } else if ( key == "-verbose" ) {
      if ( ++arg < args.size() ) verbose = clipper::String(args[arg]).i();
    } else {
      std::cout << "\nUnrecognized:\t" << args[arg] << std::endl;
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "\nUsage: cparrot\n\t-mtzin-ref <filename>\n\t-pdbin-ref <filename>\n\t-mtzin <filename>\t\tCOMPULSORY\n\t-seqin <filename>\n\t-pdbin <filename>\n\t-pdbin-ha <filename>\n\t-pdbin-mr <filename>\n\t-colin-ref-fo <colpath>\n\t-colin-ref-hl <colpath>\n\t-colin-fo <colpath>\t\tCOMPULSORY\n\t-colin-hl <colpath> or -colin-phifom <colpath>\tCOMPULSORY\n\t-colin-fc <colpath>\n\t-colin-free <colpath>\n\t-mtzout <filename>\n\t-colout <colpath>\n\t-colout-hl <colpath>\n\t-colout-fc <colpath>\n\t-mapout-ncs <filename prefix>\n\t-solvent-flatten\n\t-histogram-match\n\t-ncs-average\n\t-rice-probability\n\t-anisotropy-correction\n\t-cycles <cycles>\n\t-resolution <resolution/A>\n\t-solvent-content <fraction>\n\t-solvent-mask-filter-radius <radius>\n\t-ncs-mask-filter-radius <radius>\n\t-ncs-asu-fraction <fraction>\n\t-ncs-operator <alpha>,<beta>,<gamma>,<x>,<y>,<z>,<x>,<y>,<z>\n\t-xmlout <filename>\nAn input mtz is specified, F's and HL coefficients are required.\n";
    exit(1);
  }

  // other initialisations
  typedef clipper::HKL_data_base::HKL_reference_index HRI;
  typedef clipper::Xmap_base::Map_reference_index MRI;
  typedef clipper::NXmap_base::Map_reference_index NRI;
  using clipper::data32::Compute_abcd_from_phifom;
  using clipper::data32::Compute_phifom_from_abcd;
  using clipper::data32::Compute_fphi_from_fsigf_phifom;
  using clipper::data32::Compute_scale_u_iso_fsigf;
  using clipper::data32::Compute_scale_u_aniso_fphi;
  using clipper::data32::F_sigF;
  using clipper::data32::F_phi;
  using clipper::data32::Phi_fom;
  using clipper::data32::ABCD;
  using clipper::data32::Flag;
  clipper::Resolution resol;
  clipper::CCP4MTZfile mtzfile;
  mtzfile.set_column_label_mode( clipper::CCP4MTZfile::Legacy );
  std::string msg;
  ParrotUtil util( title );
  if ( opcol_hl == "NONE" ) opcol_hl = opcol;
  if ( opcol_fc == "NONE" ) opcol_fc = opcol;
  if ( !( dosolv || dohist || doncsa ) ) dosolv = dohist = doncsa = true;
  if ( ipmtz_ref == "NONE" || ippdb_ref == "NONE" )
    ParrotUtil::set_reference( ipmtz_ref, ippdb_ref );
  std::cout << std::endl << std::endl;

  // Get coordinate models
  clipper::MiniMol mol_ref, mol_wrk, mol_wrk_ha, mol_wrk_mr;
  ParrotUtil::read_model( mol_ref, ippdb_ref, verbose>1 );
  ParrotUtil::read_model( mol_wrk, ippdb_wrk, verbose>1 );
  ParrotUtil::read_model( mol_wrk_ha, ippdb_wrk_ha, verbose>1 );
  ParrotUtil::read_model( mol_wrk_mr, ippdb_wrk_mr, verbose>1 );
  if ( mol_wrk.size() == 0 ) domask = true;
  std::cout << std::endl;

  // Get resolution for calculation
  mtzfile.open_read( ipmtz_ref );
  double res_ref = std::max( mtzfile.resolution().limit(), res_in );
  mtzfile.close_read();
  mtzfile.open_read( ipmtz_wrk );
  double res_wrk = std::max( mtzfile.resolution().limit(), res_in );
  mtzfile.close_read();
  resol = clipper::Resolution( std::max( res_ref, res_wrk ) );
  if ( res_ref > res_wrk ) std::cout << "\nWARNING: resolution of work structure truncated to reference:\n Ref: " << res_ref << " Wrk: " << res_wrk << std::endl;

  // Get reference reflection data
  clipper::HKL_info hkls_ref;
  mtzfile.open_read( ipmtz_ref );
  hkls_ref.init( mtzfile.spacegroup(), mtzfile.cell(), resol, true );
  clipper::HKL_data<F_sigF> ref_f( hkls_ref );
  clipper::HKL_data<ABCD>   ref_hl( hkls_ref );
  mtzfile.import_hkl_data( ref_f,  ipcol_ref_fo );
  mtzfile.import_hkl_data( ref_hl, ipcol_ref_hl );
  mtzfile.close_read();

  // Get work reflection data
  clipper::MTZcrystal cxtl;
  clipper::HKL_info hkls_wrk;
  mtzfile.set_verbose( (verbose>0) ? 3 : 2 );
  mtzfile.open_read( ipmtz_wrk );
  hkls_wrk.init( mtzfile.spacegroup(), mtzfile.cell(), resol, true );
  mtzfile.import_crystal( cxtl, ipcol_wrk_fo+".F_sigF.F" );
  clipper::HKL_data<F_sigF>  wrk_f ( hkls_wrk, cxtl );
  clipper::HKL_data<ABCD>    wrk_hl( hkls_wrk, cxtl );
  clipper::HKL_data<Phi_fom> wrk_pw( hkls_wrk, cxtl );
  clipper::HKL_data<F_phi>   wrk_fp( hkls_wrk, cxtl );
  clipper::HKL_data<Flag>    flag( hkls_wrk, cxtl );
  mtzfile.import_hkl_data( wrk_f , ipcol_wrk_fo );
  if ( ipcol_wrk_hl != "NONE" ) mtzfile.import_hkl_data( wrk_hl,ipcol_wrk_hl );
  if ( ipcol_wrk_pw != "NONE" ) mtzfile.import_hkl_data( wrk_pw,ipcol_wrk_pw );
  if ( ipcol_wrk_fc != "NONE" ) mtzfile.import_hkl_data( wrk_fp,ipcol_wrk_fc );
  if ( ipcol_wrk_fr != "NONE" ) mtzfile.import_hkl_data( flag,  ipcol_wrk_fr );
  clipper::String oppath = mtzfile.assigned_paths()[0].notail() + "/";
  mtzfile.close_read();

  // do anisotropy correction
  clipper::U_aniso_orth uaniso( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
  if ( doanis ) {
    // scale obs data
    typedef clipper::SFscale_aniso<float> SFscale;
    SFscale sfscl( 3.0, SFscale::SHARPEN );
    sfscl( wrk_f );
    uaniso = sfscl.u_aniso_orth( SFscale::F );
    // scale map coeffs
    Compute_scale_u_aniso_fphi compute_aniso( 1.0, -uaniso );
    if ( ipcol_wrk_fc != "NONE" ) wrk_fp.compute( wrk_fp, compute_aniso );    
    // output
    std::cout << std::endl << "Applying anisotropy correction:"
              << std::endl << uaniso.format() << std::endl;
  }

  if ( usharp != 0.0 )
    wrk_f.compute( wrk_f, Compute_scale_u_iso_fsigf( 1.0, usharp ) );
  if ( ipcol_wrk_hl == "NONE" )
    wrk_hl.compute( wrk_pw, Compute_abcd_from_phifom() );

  // Get work sequence (optional)
  clipper::MMoleculeSequence seq_wrk;
  if ( ipseq_wrk != "NONE" ) {
    clipper::SEQfile seqf_wrk;
    seqf_wrk.read_file( ipseq_wrk );
    seqf_wrk.import_molecule_sequence( seq_wrk );
  }

  // crystal data
  clipper::Spacegroup spgr_ref = hkls_ref.spacegroup();
  clipper::Cell cell_ref = hkls_ref.cell();
  clipper::Grid_sampling grid_ref( spgr_ref, cell_ref, resol );
  clipper::Spacegroup spgr_wrk = hkls_wrk.spacegroup();
  clipper::Cell cell_wrk = cxtl;
  clipper::Grid_sampling grid_wrk( spgr_wrk, cell_wrk, resol );

  // calculate reference mask
  clipper::Xmap<float> msk_ref( spgr_ref, cell_ref, grid_ref );
  ParrotUtil::solvent_mask( msk_ref, mol_ref );

  // calculate work mask (if required)
  clipper::Xmap<float> msk_prt( spgr_wrk, cell_wrk, grid_wrk );
  clipper::Xmap<float> msk_sol( spgr_wrk, cell_wrk, grid_wrk );
  if ( mol_wrk.size() > 0 ) {
    ParrotUtil::solvent_mask( msk_prt, mol_wrk );
    for ( MRI ix = msk_prt.first(); !ix.last(); ix.next() )
      msk_sol[ix] = 1.0 - msk_prt[ix];
    double solest = 1.0 - clipper::Map_stats( msk_prt ).mean();
    if ( solc <= 0.0 ) solc = solest;
    // OUTPUT
    prog.summary_beg();
    std::cout << std::endl << "Solvent content from model: " << solest
              << std::endl << std::endl;
    prog.summary_end(); std::cout << std::endl;
  }

  // solvent content
  if ( seq_wrk.size() > 0 ) {
    std::vector<std::pair<double,double> > solcs = ParrotUtil::solvent_probability( seq_wrk, hkls_wrk.spacegroup(), hkls_wrk.cell(), resol );
    nncs = 0;
    if ( solc <= 0.0 ) {
      for ( int n = 0; n < solcs.size(); n++ )
        if ( solcs[n].first > solcs[nncs].first ) nncs = n;
      solc = solcs[nncs].second;
    }
    // OUTPUT
    solcs[0].first = double(nncs);
    solcs[0].second = solc;
    prog.summary_beg();
    util.log_solvent_content( solcs );
    prog.summary_end(); std::cout << std::endl;
  }

  // other data objects
  clipper::Xmap<float> map_ref( spgr_ref, cell_ref, grid_ref );
  clipper::Xmap<float> map_sim( spgr_ref, cell_ref, grid_ref );
  clipper::Xmap<float> map_wrk( spgr_wrk, cell_wrk, grid_wrk );
  clipper::Xmap<float> map_mod( spgr_wrk, cell_wrk, grid_wrk );
  clipper::Xmap<float> map_wrk_gamma( spgr_wrk, cell_wrk, grid_wrk );
  clipper::Xmap<float> map_mod_gamma( spgr_wrk, cell_wrk, grid_wrk );
  clipper::Xmap<float> map_ncs( spgr_wrk, cell_wrk, grid_wrk );
  clipper::Xmap<float> map_nwt( spgr_wrk, cell_wrk, grid_wrk );
  clipper::HKL_data<Phi_fom> sim_pw( hkls_ref ), ref_pw( hkls_ref );
  clipper::HKL_data<F_phi>   sim_fp( hkls_ref ), ref_fp( hkls_ref );
  ref_pw.compute( ref_hl, Compute_phifom_from_abcd() );
  ref_fp.compute( ref_f, ref_pw, Compute_fphi_from_fsigf_phifom() );
  clipper::HKL_data<F_phi> fbest( hkls_wrk ), fdiff( hkls_wrk );
  clipper::HKL_data<Flag> flagwt( hkls_wrk );
  clipper::HKL_data<ABCD> mod_hl( hkls_wrk );
  for ( HRI ih = flagwt.first(); !ih.last(); ih.next() )
    if ( !wrk_f[ih].missing() )
      flagwt[ih].flag() = clipper::SFweight_spline<float>::BOTH;

  // work map coefficients
  wrk_pw.compute( wrk_hl, Compute_phifom_from_abcd() );
  if ( ipcol_wrk_fc == "NONE" )
    wrk_fp.compute( wrk_f, wrk_pw, Compute_fphi_from_fsigf_phifom() );

  // get NCS stats
  map_wrk.fft_from( wrk_fp );
  const clipper::Map_stats ncsstats = NCSaver::ncs_stats( map_wrk, ncs_radius );
  const double ncs_level = 4.0*ncsstats.std_dev();

  // get NCS operators
  if ( doncsa && mol_wrk_ha.size() > 0 ) {
    std::cout << std::endl << "NCS from heavy atoms: " << std::endl;
    std::vector<Local_rtop> atmops, newops;
    const double tol_dst = 3.0;
    const double tol_ang = 0.1;
    NCSfind ncsatom( tol_dst, tol_ang );
    NCSaver ncsaver;
    atmops = ncsatom.find_ncs_candidates( mol_wrk_ha.atom_list(),
                                          spgr_wrk, cell_wrk );
    newops = ncsaver.filter_ncs_candidates( atmops, map_wrk, ncs_radius,
                                            ncs_level, ncs_asufrc, 3 );
    if ( verbose > 0 ) Local_rtop::print_nxops( "from heavy atoms", atmops );
    if ( verbose > 0 ) Local_rtop::print_nxops( "density filtered", newops );
    for ( int r = 0; r < newops.size(); r++ ) {
      std::cout << std::endl << "NCS operator found relating heavy atoms: "
                << std::endl;
      util.log_ncs_operator( newops[r] );
    }
    if ( newops.size() == 0 ) std::cout << "$TEXT:Warning: $$ $$\nWARNING: No NCS found from heavy atoms.\n$$" << std::endl;
    else nxops.insert( nxops.end(), newops.begin(), newops.end() );
    std::cout << std::endl;
    util.log_ncs_ha( newops.size() );
  }
  if ( doncsa && mol_wrk_mr.size() > 0 ) {
    std::cout << std::endl << "NCS from atomic model: " << std::endl;
    std::vector<Local_rtop> newops;
    for ( int c1 = 0; c1 < mol_wrk_mr.size(); c1++ )
      for ( int c2 = 0; c2 < mol_wrk_mr.size(); c2++ )
        if ( c1 != c2 ) {
          const double rmsd = 1.0;
          Local_rtop nxop =
            NCSfind::local_rtop( mol_wrk_mr[c1], mol_wrk_mr[c2], rmsd, 16 );
          if ( !nxop.is_null() ) {
            newops.push_back( nxop );
            std::cout << std::endl << "NCS operator found relating chains "
                      << mol_wrk_mr[c1].id() << " and "
                      << mol_wrk_mr[c2].id() << std::endl;
            util.log_ncs_operator( nxop );
          }
        }
    if ( newops.size() == 0 ) std::cout << "$TEXT:Warning: $$ $$\nWARNING: No NCS found from atomic model.\n$$" << std::endl;
    else nxops.insert( nxops.end(), newops.begin(), newops.end() );
    std::cout << std::endl;
    util.log_ncs_mr( newops.size() );
  }

  // store initial stats
  util.log_rfl_stats( wrk_f, wrk_fp, wrk_pw, flagwt );

  // ----------------------------------------------------------------------
  // Main program loop
  for ( int cycle = 0; cycle < ncycles; cycle++ ) {
    util.log_cycle( cycle+1 );

    // reference map coefficients
    clipper::HKL_data<F_sigF> sim_f( hkls_ref );
    clipper::HKL_data<ABCD> sim_hl( hkls_ref );
    MapSimulate mapsim( 100, 20 );
    mapsim( sim_f, sim_hl, ref_f, ref_hl, wrk_f, wrk_hl );
    sim_pw.compute( sim_hl, Compute_phifom_from_abcd() );
    sim_fp.compute( sim_f, sim_pw, Compute_fphi_from_fsigf_phifom() );
    for ( HRI ih = ref_fp.first(); !ih.last(); ih.next() )
      ref_fp[ih] = F_phi( sim_f[ih].f(), ref_fp[ih].phi() );

    // map calculation
    map_sim.fft_from( sim_fp );
    map_ref.fft_from( ref_fp );
    map_wrk.fft_from( wrk_fp );

    // mask radius
    double autorad = 1.0 * ParrotUtil::effective_resolution( wrk_pw );
    std::cout << "Suggested radius for solvent mask determination: "
              << autorad << std::endl << std::endl;

    // mask calculation from density
    if ( domask ) {
      clipper::Xmap<float> lmom1( map_wrk ), lmom2( map_wrk );
      for ( MRI ix = lmom2.first(); !ix.last(); ix.next() )
        lmom2[ix] = clipper::Util::sqr( lmom2[ix] );
      // now calculate local mom1, local mom1 squared
      clipper::MapFilterFn_step fn( solrad > 0.0 ? solrad : autorad );
      clipper::MapFilter_fft<float>
        fltr( fn, 1.0, clipper::MapFilter_fft<float>::Relative );
      fltr( lmom1, lmom1 );
      fltr( lmom2, lmom2 );
      // calculate std deviation
      for ( MRI ix = lmom1.first(); !ix.last(); ix.next() )
        lmom2[ix] = sqrt( lmom2[ix] - clipper::Util::sqr( lmom1[ix] ) );
      // now make ordinal and cutoff
      clipper::Map_stats stat_sd( lmom2 );
      clipper::Generic_ordinal ord_msk( stat_sd.range(), 200 );
      for ( MRI ix = lmom2.first(); !ix.last(); ix.next() )
        ord_msk.accumulate( lmom2[ix] );
      ord_msk.prep_ordinal();
      ord_msk.invert();
      double solcut = ord_msk.ordinal( solc );
      for ( MRI ix = msk_prt.first(); !ix.last(); ix.next() ) {
        msk_prt[ix] = ( lmom2[ix] > solcut ) ? 1.0 : 0.0;
        msk_sol[ix] = 1.0 - msk_prt[ix];
      }
    }

    // do NCS calc
    map_ncs = 0.0;
    map_nwt = 0.0;
    bool ncsref = ( cycle <= 1 ) || ( cycle%3 == 0 );
    NCSaver ncsaver;
    if ( doncsa && nxops.size() > 0 ) {
      std::vector<double> ncsvols( nxops.size() ), ncscors( nxops.size() );
      for ( int r = 0; r < nxops.size(); r++ ) {
        Local_rtop nxop = nxops[r];
        Local_rtop nxop_old = nxop;
        clipper::NXmap<float> ncsmsk;
        const double asuvol = map_wrk.cell().volume() /
          double( map_wrk.spacegroup().num_symops() );
        const double map_radius = pow( asuvol, 0.333 );
        ncsaver.ncs_mask( ncsmsk, map_wrk, nxop,
                          map_radius, ncs_radius, ncs_level, 3 );
        if ( ncsref ) ncsaver.ncs_refine( nxop, map_wrk, ncsmsk );
        ncsaver.ncs_average( map_ncs, map_nwt, map_wrk, ncsmsk, nxop );
        nxops[r] = nxop;

        // store and log output
        util.log_ncs_stats
          ( nxop_old, nxop, ncsaver.mask_volume_asu(),
            ncsaver.mask_multiplicity(), ncsaver.correlation_sphere(),
            ncsaver.correlation_old(), ncsaver.correlation_new(),
            ncsaver.mask_volume_ratio(), ncsaver.mask_overlap_ratio() );
        if ( opncsm != "NONE" )
          util.output_ncs_mask( opncsm, ncsmsk, cell_wrk, cycle, r );
      }

      // output
      prog.summary_beg();
      util.log_ncs_table();
      prog.summary_end(); std::cout << std::endl;
    }

    // density modify
    ParrotUtil::density_modify( map_mod, map_wrk, map_ncs, map_nwt, msk_prt, msk_sol, map_ref, map_sim, msk_ref, dosolv, dohist );

    // gamma correct
    ParrotUtil::Map_stats statp = ParrotUtil::masked_stats( map_wrk, msk_prt );
    double sclrnd = 0.01 * statp.std_dev();
    for ( MRI ix = msk_prt.first(); !ix.last(); ix.next() )
      map_wrk_gamma[ix] = map_wrk[ix] + sclrnd * ParrotUtil::random();
    ParrotUtil::density_modify( map_mod_gamma, map_wrk_gamma, map_ncs, map_nwt, msk_prt, msk_sol, map_ref, map_sim, msk_ref, dosolv, dohist );
    double s11(0.0), s12(0.0);
    for ( MRI ix = msk_prt.first(); !ix.last(); ix.next() ) {
      map_wrk_gamma[ix] -= map_wrk[ix];
      map_mod_gamma[ix] -= map_mod[ix];
      s11 += map_wrk_gamma[ix] * map_wrk_gamma[ix];
      s12 += map_mod_gamma[ix] * map_wrk_gamma[ix];
    }
    double gamma = s12 / s11;

    // display stats
    util.log_histogram_graph( msk_ref, msk_prt, msk_sol, map_ref, map_sim, map_wrk, map_mod );

    // apply gamma
    std::cout << "Gamma " << gamma << std::endl << std::endl;
    for ( MRI ix = msk_prt.first(); !ix.last(); ix.next() )
      map_mod[ix] -= gamma * map_wrk[ix];

    // back transform
    map_mod.fft_to( wrk_fp );
    if ( dodump ) break;

    // weighting
    clipper::SFweight_spline<float> sfw( n_refln, n_param );
    if ( !dorice ) {
      // mlhl phase weighting calc
      sfw( fbest, fdiff, wrk_pw, mod_hl, wrk_f, wrk_hl, wrk_fp, flagwt );
      wrk_fp = fbest;
    } else {
      // rice phase weighting calc
      sfw( fbest, fdiff, wrk_pw, wrk_f, wrk_fp, flagwt );
      mod_hl.compute( wrk_pw, Compute_abcd_from_phifom() );
      mod_hl = wrk_hl + mod_hl;
      // work map coefficients
      wrk_pw.compute( mod_hl, Compute_phifom_from_abcd() );
      wrk_fp.compute( wrk_f, wrk_pw, Compute_fphi_from_fsigf_phifom() );
      // restore missing data
      for ( HRI ih = wrk_fp.first(); !ih.last(); ih.next() )
        if ( wrk_fp[ih].missing() ) wrk_fp[ih] = fbest[ih];
    }

    // display stats
    util.log_sigmaa_graph( sfw, flagwt );
    util.log_rfl_stats( wrk_f, wrk_fp, wrk_pw, flagwt );
    if ( opxml != "NONE" ) util.xml( opxml );
  } // cycle loop
  // ----------------------------------------------------------------------

  // output new results
  Compute_scale_u_aniso_fphi compute_aniso( 1.0, uaniso );
  wrk_fp.compute( wrk_fp, compute_aniso );
  if ( opcol_hl[0] != '/' ) opcol_hl = oppath + opcol_hl;
  if ( opcol_fc[0] != '/' ) opcol_fc = oppath + opcol_fc;
  mtzfile.open_append( ipmtz_wrk, opfile );
  if ( !dodump ) {
    mtzfile.export_hkl_data( mod_hl, opcol_hl );
    mtzfile.export_hkl_data( wrk_fp, opcol_fc );
  } else {
    mtzfile.export_hkl_data( wrk_fp, opcol_fc );
  }
  mtzfile.close_append();

  // output NCS operators
  if ( nxops.size() > 0 ) {
    std::cout << "Non-crystallographic operators:" << std::endl;
    for ( int r = 0; r < nxops.size(); r++ ) {
      Local_rtop nxop = nxops[r];
      clipper::Rotation rot( nxop.rtop_orth().rot() );
      clipper::Euler_ccp4 euler = rot.euler_ccp4();
      std::cout << " -ncs-operator "
                << clipper::Util::rad2d(euler.alpha()) << ","
                << clipper::Util::rad2d(euler.beta() ) << ","
                << clipper::Util::rad2d(euler.gamma()) << ","
                << nxop.src().x() << "," << nxop.src().y() << ","
                << nxop.src().z() << "," << nxop.tgt().x() << ","
                << nxop.tgt().y() << "," << nxop.tgt().z()
                << std::endl;
    }
  }
  std::cout << std::endl;

  util.log_summary_graphs();

  prog.set_termination_message( "Normal termination" );
}
