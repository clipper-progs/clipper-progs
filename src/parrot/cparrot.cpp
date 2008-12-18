// Clipper parrot
/* Copyright 2008 Kevin Cowtan & University of York all rights reserved */


#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-minimol.h>
#include <clipper/clipper-contrib.h>

#include "simulate-lib.h"
#include "parrot-lib.h"
#include "parrot-ncsfind.h"
#include "parrot-ncsaver.h"


int main( int argc, char** argv )
{
  CCP4Program prog( "cparrot", "1.0.0", "$Date: 2008/10/30" );
  prog.set_termination_message( "Failed" );

  std::cout << std::endl << "Copyright 2008 Kevin Cowtan and University of York." << std::endl << std::endl;
  prog.summary_beg();
  std::cout << "$TEXT:Reference: $$ Please reference $$" << std::endl << std::endl << " 'Combining constraints for electron-density modification.'" << std::endl << " Zhang K. Y. J., Cowtan K., Main P. (1997) Methods in Enzymology, 277, 53-64." << std::endl << std::endl << "$$";
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
  int nncs = 1;
  int ncycles = 3;
  int freeflag = 0;
  int n_refln = 2000;
  int n_param = 10;
  int verbose = 0;

  // command input
  CCP4CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if        ( args[arg] == "-title" ) {
      if ( ++arg < args.size() ) title = args[arg];
    } else if ( args[arg] == "-mtzin-ref" ) {
      if ( ++arg < args.size() ) ipmtz_ref = args[arg];
    } else if ( args[arg] == "-mtzin-wrk" ) {
      if ( ++arg < args.size() ) ipmtz_wrk = args[arg];
    } else if ( args[arg] == "-pdbin-ref" ) {
      if ( ++arg < args.size() ) ippdb_ref = args[arg];
    } else if ( args[arg] == "-seqin-wrk" ) {
      if ( ++arg < args.size() ) ipseq_wrk = args[arg];
    } else if ( args[arg] == "-pdbin-wrk" ) {
      if ( ++arg < args.size() ) ippdb_wrk = args[arg];
    } else if ( args[arg] == "-pdbin-wrk-ha" ) {
      if ( ++arg < args.size() ) ippdb_wrk_ha  = args[arg];
    } else if ( args[arg] == "-pdbin-wrk-mr" ) {
      if ( ++arg < args.size() ) ippdb_wrk_mr  = args[arg];
    } else if ( args[arg] == "-colin-ref-fo" ) {
      if ( ++arg < args.size() ) ipcol_ref_fo = args[arg];
    } else if ( args[arg] == "-colin-ref-hl" ) {
      if ( ++arg < args.size() ) ipcol_ref_hl = args[arg];
    } else if ( args[arg] == "-colin-wrk-fo" ) {
      if ( ++arg < args.size() ) ipcol_wrk_fo = args[arg];
    } else if ( args[arg] == "-colin-wrk-hl" ) {
      if ( ++arg < args.size() ) ipcol_wrk_hl = args[arg];
    } else if ( args[arg] == "-colin-wrk-phifom" ) {
      if ( ++arg < args.size() ) ipcol_wrk_pw = args[arg];
    } else if ( args[arg] == "-colin-wrk-fc" ) {
      if ( ++arg < args.size() ) ipcol_wrk_fc = args[arg];
    } else if ( args[arg] == "-colin-wrk-free" ) {
      if ( ++arg < args.size() ) ipcol_wrk_fr = args[arg];
    } else if ( args[arg] == "-mtzout" ) {
      if ( ++arg < args.size() ) opfile = args[arg];
    } else if ( args[arg] == "-colout" ) {
      if ( ++arg < args.size() ) opcol = args[arg];
    } else if ( args[arg] == "-colout-hl" ) {
      if ( ++arg < args.size() ) opcol_hl = args[arg];
    } else if ( args[arg] == "-colout-fc" ) {
      if ( ++arg < args.size() ) opcol_fc = args[arg];
    } else if ( args[arg] == "-solvent-flatten" ) {
      dosolv = true;
    } else if ( args[arg] == "-histogram-match" ) {
      dohist = true;
    } else if ( args[arg] == "-ncs-average" ) {
      doncsa = true;
    } else if ( args[arg] == "-force-solvent-mask-calculation" ) {
      domask = true;
    } else if ( args[arg] == "-rice-probability" ) {
      dorice = true;
    } else if ( args[arg] == "-resolution" ) {
      if ( ++arg < args.size() ) res_in = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-cycles" ) {
      if ( ++arg < args.size() ) ncycles = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-anisotropy-correction" ) {
      doanis = true;
    } else if ( args[arg] == "-solvent-content" ) {
      if ( ++arg < args.size() ) solc   = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-solvent-mask-filter-radius" ) {
      if ( ++arg < args.size() ) solrad = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-ncs-mask-filter-radius" ) {
      if ( ++arg < args.size() ) ncs_radius = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-ncs-asu-fraction" ) {
      if ( ++arg < args.size() ) ncs_asufrc = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-ncs-operator" ) {
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
    } else if ( args[arg] == "-sharpen" ) {
      if ( ++arg < args.size() ) usharp = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-verbose" ) {
      if ( ++arg < args.size() ) verbose = clipper::String(args[arg]).i();
    } else {
      std::cout << "\nUnrecognized:\t" << args[arg] << std::endl;
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "\nUsage: cparrot\n\t-mtzin-ref <filename>\n\t-pdbin-ref <filename>\n\t-mtzin-wrk <filename>\t\tCOMPULSORY\n\t-seqin-wrk <filename>\n\t-pdbin-wrk <filename>\n\t-pdbin-wrk-ha <filename>\n\t-pdbin-wrk-mr <filename>\n\t-colin-ref-fo <colpath>\n\t-colin-ref-hl <colpath>\n\t-colin-wrk-fo <colpath>\t\tCOMPULSORY\n\t-colin-wrk-hl <colpath> or -colin-wrk-phifom <colpath>\tCOMPULSORY\n\t-colin-wrk-fc <colpath>\n\t-colin-wrk-free <colpath>\n\t-mtzout <filename>\n\t-colout <colpath>\n\t-colout-hl <colpath>\n\t-colout-fc <colpath>\n\t-solvent-flatten\n\t-histogram-match\n\t-ncs-average\n\t-rice-probability\n\t-do-anisotropy-correction\n\t-cycles <cycles>\n\t-resolution <resolution/A>\n\t-solvent-content <fraction>\n\t-solvent-mask-filter-radius <radius>\n\t-ncs-mask-filter-radius <radius>\n\t-ncs-asu-fraction <fraction>\n\t-ncs-operator <alpha>,<beta>,<gamma>,<x>,<y>,<z>,<x>,<y>,<z>\nAn input mtz is specified, F's and HL coefficients are required.\n";
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
  std::string msg;
  if ( opcol_hl == "NONE" ) opcol_hl = opcol;
  if ( opcol_fc == "NONE" ) opcol_fc = opcol;
  if ( !( dosolv || dohist || doncsa ) ) dosolv = dohist = doncsa = true;
  if ( ipmtz_ref == "NONE" || ippdb_ref == "NONE" )
    ParrotUtil::set_reference( ipmtz_ref, ippdb_ref );
  std::cout << std::endl << std::endl;

  // Get coordinate models
  clipper::MiniMol mol_ref, mol_wrk, mol_wrk_ha, mol_wrk_mr;
  ParrotUtil::read_model( mol_ref, ippdb_ref );
  ParrotUtil::read_model( mol_wrk, ippdb_wrk );
  ParrotUtil::read_model( mol_wrk_ha, ippdb_wrk_ha );
  ParrotUtil::read_model( mol_wrk_mr, ippdb_wrk_mr );
  if ( mol_wrk.size() == 0 ) domask = true;
  std::cout << std::endl;

  // Get resolution for calculation
  mtzfile.open_read( ipmtz_ref );
  double res_ref = clipper::Util::max( mtzfile.resolution().limit(), res_in );
  mtzfile.close_read();
  mtzfile.open_read( ipmtz_wrk );
  double res_wrk = clipper::Util::max( mtzfile.resolution().limit(), res_in );
  mtzfile.close_read();
  resol = clipper::Resolution( clipper::Util::max( res_ref, res_wrk ) );
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
    for ( int n = 0; n < solcs.size(); n++ )
      if ( solcs[n].first > solcs[nncs].first ) nncs = n;
    double solest = solcs[nncs].second;
    if ( solc <= 0.0 ) solc = solest;
    // OUTPUT
    prog.summary_beg();
    printf("\nSolvent content estimation from sequence:\n");
    printf(" N(NCS)   Solvent_fraction    Probability\n");
    for ( int n = 1; n < solcs.size(); n++ )
      printf("    %2i       %8.3f         %8.3f\n",
	     n, solcs[n].second, solcs[n].first );
    std::cout << std::endl << "Solvent content from sequence: " << solest
	      << "   (assuming " << nncs << "-fold NCS)."
	      << std::endl << std::endl;
    if ( solcs.size() > 2 && solcs[nncs].first < 0.98 )
      std::cout << "$TEXT:Warning: $$ $$" << std::endl
		<< "WARNING: Assuming " << nncs << "-fold NCS. "
		<< "This is only a guess." << std::endl
		<< "Consider other possibilities "
		<< "and set solvent content accordingly" << std::endl
		<< "$$" << std::endl;
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
  wrk_fp.compute( wrk_f, wrk_pw, Compute_fphi_from_fsigf_phifom() );

  // get NCS stats
  map_wrk.fft_from( wrk_fp );
  const clipper::Map_stats ncsstats = NCSaver::ncs_stats( map_wrk, ncs_radius );
  const double ncs_level = 4.0*ncsstats.std_dev();

  // get NCS operators
  if ( doncsa && mol_wrk_ha.size() > 0 ) {
    std::cout << std::endl << "NCS from heavy atoms: " << std::endl;
    const double tol_dst = 3.0;
    const double tol_ang = 0.1;
    NCSfind ncsatom( tol_dst, tol_ang );
    NCSaver ncsaver;
    nxops = ncsatom.find_ncs_candidates( mol_wrk_ha.atom_list(),
					 spgr_wrk, cell_wrk );
    nxops = ncsaver.filter_ncs_candidates( nxops, map_wrk, ncs_radius,
					   ncs_level, ncs_asufrc, 3 );
    for ( int r = 0; r < nxops.size(); r++ ) {
      clipper::Rotation rot( nxops[r].rot() );
      clipper::Polar_ccp4 polar = rot.polar_ccp4();
      clipper::Euler_ccp4 euler = rot.euler_ccp4();
      std::cout << std::endl << "NCS operator found relating heavy atoms: "
		<< std::endl;
      std::cout << " Polar rotation/deg: "
		<< clipper::Util::rad2d(polar.omega()) << ","
		<< clipper::Util::rad2d(polar.phi() ) << ","
		<< clipper::Util::rad2d(polar.kappa()) << std::endl;
      std::cout << " Euler rotation/deg: "
		<< clipper::Util::rad2d(euler.alpha()) << ","
		<< clipper::Util::rad2d(euler.beta() ) << ","
		<< clipper::Util::rad2d(euler.gamma()) << std::endl;
      std::cout << " Source: " << nxops[r].src().format() << std::endl;
      std::cout << " Target: " << nxops[r].tgt().format() << std::endl;
    }
    std::cout << std::endl;
  }
  if ( doncsa && mol_wrk_mr.size() > 0 ) {
    std::cout << std::endl << "NCS from atomic model: " << std::endl;
    for ( int c1 = 0; c1 < mol_wrk_mr.size(); c1++ )
      for ( int c2 = 0; c2 < mol_wrk_mr.size(); c2++ )
	if ( c1 != c2 ) {
	  const double rmsd = 1.0;
	  Local_rtop nxop =
	    NCSfind::local_rtop( mol_wrk_mr[c1], mol_wrk_mr[c2], rmsd, 16 );
	  if ( !nxop.is_null() ) {
	    nxops.push_back( nxop );
	    clipper::Rotation rot( nxop.rot() );
	    clipper::Polar_ccp4 polar = rot.polar_ccp4();
	    clipper::Euler_ccp4 euler = rot.euler_ccp4();
	    std::cout << std::endl << "NCS operator found relating chains "
		      << mol_wrk_mr[c1].id() << " and "
		      <<  mol_wrk_mr[c2].id() << std::endl;
	    std::cout << " Polar rotation/deg: "
		      << clipper::Util::rad2d(polar.omega()) << ","
		      << clipper::Util::rad2d(polar.phi() ) << ","
		      << clipper::Util::rad2d(polar.kappa()) << std::endl;
	    std::cout << " Euler rotation/deg: "
		      << clipper::Util::rad2d(euler.alpha()) << ","
		      << clipper::Util::rad2d(euler.beta() ) << ","
		      << clipper::Util::rad2d(euler.gamma()) << std::endl;
	    std::cout << " Source: " << nxop.src().format() << std::endl;
	    std::cout << " Target: " << nxop.tgt().format() << std::endl;
	  }
	}
    std::cout << std::endl << std::endl;
  }

  // Main program loop
  for ( int cycle = 0; cycle < ncycles; cycle++ ) {
    std::cout << std::endl << "-- Cycle: " << cycle+1
	      << " --------------------------------" << std::endl << std::endl;

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
	Local_rtop nxo0 = nxop;
	clipper::NXmap<float> ncsmsk;
	const double asuvol = map_wrk.cell().volume() /
	  double( map_wrk.spacegroup().num_symops() );
	const double map_radius = pow( asuvol, 0.333 );
	ncsaver.ncs_mask( ncsmsk, map_wrk, nxop,
			  map_radius, ncs_radius, ncs_level, 3 );
	ncsaver.ncs_refine( nxop, map_wrk, ncsmsk, ncsref );
	ncsaver.ncs_average( map_ncs, map_nwt, map_wrk, ncsmsk, nxop );
	
	// store refined operator
	nxops[r] = nxop;
	ncsvols[r] = ncsaver.mask_volumes()[0]/asuvol;
	ncscors[r] = ncsaver.correlation_new();

	// output
	printf( "NCS operator: %3i\n", r );
	printf( " NCS masking: mask volume as fraction of ASU: %8.2f\n",
		ncsvols[r] );
	if ( ncsref ) {
	  clipper::Euler_ccp4 euler1 = nxo0.rot().euler_ccp4();
	  clipper::Euler_ccp4 euler2 = nxop.rot().euler_ccp4();
	  printf( " NXop refinement- correlation before: %6.3f, after: %6.3f\n", ncsaver.correlation_old(), ncsaver.correlation_new() );
	  printf( " NXop old: %6.1f %6.1f %6.1f   %6.1f %6.1f %6.1f   %6.1f %6.1f %6.1f\n", clipper::Util::rad2d(euler1.alpha()), clipper::Util::rad2d(euler1.beta() ), clipper::Util::rad2d(euler1.gamma()), nxo0.src().x(), nxo0.src().y(), nxo0.src().z(), nxo0.tgt().x(), nxo0.tgt().y(), nxo0.tgt().z() );
	  printf( " NXop new: %6.1f %6.1f %6.1f   %6.1f %6.1f %6.1f   %6.1f %6.1f %6.1f\n", clipper::Util::rad2d(euler2.alpha()), clipper::Util::rad2d(euler2.beta() ), clipper::Util::rad2d(euler2.gamma()), nxop.src().x(), nxop.src().y(), nxop.src().z(), nxop.tgt().x(), nxop.tgt().y(), nxop.tgt().z() );
	}
	std::cout << std::endl;
      }

      // output
      prog.summary_beg();
      printf("NCS operator statistics:\n");
      printf(" Operator_number  Mask_volume/ASU  Correlation\n");
      for ( int n = 0; n < ncsvols.size(); n++ )
        printf( "         %4i         %8.3f     %8.3f\n",
		n+1, ncsvols[n], ncscors[n] );
      prog.summary_end(); std::cout << std::endl;
    }

    // density modify
    ParrotUtil::density_modify( map_mod, map_wrk, map_ncs, map_nwt, msk_prt, msk_sol, map_ref, map_sim, msk_ref, dosolv, dohist );

    // make stats
    ParrotUtil::Map_stats stat_ref, stat_sim, stat_prt, stat_sol;
    stat_ref = ParrotUtil::masked_stats( map_ref, msk_ref );
    stat_sim = ParrotUtil::masked_stats( map_sim, msk_ref );
    stat_prt = ParrotUtil::masked_stats( map_wrk, msk_prt );
    stat_sol = ParrotUtil::masked_stats( map_wrk, msk_sol );
    clipper::Generic_ordinal ord_ref, ord_sim, ord_prt, ord_sol, ord_mod_prt, ord_mod_sol;
    ord_ref = ParrotUtil::masked_ordinal( map_ref, msk_ref, stat_ref );
    ord_sim = ParrotUtil::masked_ordinal( map_sim, msk_ref, stat_ref );
    ord_prt = ParrotUtil::masked_ordinal( map_wrk, msk_prt, stat_prt );
    ord_sol = ParrotUtil::masked_ordinal( map_wrk, msk_sol, stat_prt );
    ord_mod_prt = ParrotUtil::masked_ordinal( map_mod, msk_prt, stat_prt );
    ord_mod_sol = ParrotUtil::masked_ordinal( map_mod, msk_sol, stat_prt );

    // gamma correct
    double sclrnd = 0.01 * stat_prt.std_dev();
    for ( MRI ix = msk_prt.first(); !ix.last(); ix.next() )
      map_wrk_gamma[ix] = map_wrk[ix] + sclrnd * ParrotUtil::random();
    ParrotUtil::density_modify( map_mod_gamma, map_wrk_gamma, map_ncs, map_nwt, msk_prt, msk_sol, map_ref, map_sim, msk_ref, dosolv, dohist );
    double s11, s12;
    for ( MRI ix = msk_prt.first(); !ix.last(); ix.next() ) {
      map_wrk_gamma[ix] -= map_wrk[ix];
      map_mod_gamma[ix] -= map_mod[ix];
      s11 += map_wrk_gamma[ix] * map_wrk_gamma[ix];
      s12 += map_mod_gamma[ix] * map_wrk_gamma[ix];
    }
    double gamma = s12 / s11;

    // print stats
    const int ntab = 20;
    const double rng = stat_prt.max()-stat_prt.min();
    const double scl = stat_sim.std_dev() / stat_prt.std_dev();
    printf("$TABLE :Electron density histograms:\n");
    printf("$GRAPHS :Protein:N:1,4,5,6::Solvent:N:1,7,8::Simulation:N:1,3,4: $$\n");
    printf("rho_min rho_max   Simulatn P_init P_trgt P_mod  S_init S_mod $$\n");
    printf("$$\n");
    for ( int i = 0; i < ntab; i++ ) {
      double r1 = stat_prt.min() + rng * double(i  ) / double(ntab);
      double r2 = stat_prt.min() + rng * double(i+1) / double(ntab);
      double rr1 = scl * ( r1 - stat_prt.mean() ) + stat_sim.mean();
      double rr2 = scl * ( r2 - stat_prt.mean() ) + stat_sim.mean();
      double rr = ord_ref.ordinal(rr2) - ord_ref.ordinal(rr1);
      double rs = ord_sim.ordinal(rr2) - ord_sim.ordinal(rr1);
      double wp = ord_prt.ordinal(r2) - ord_prt.ordinal(r1);
      double ws = ord_sol.ordinal(r2) - ord_sol.ordinal(r1);
      double mp = ord_mod_prt.ordinal(r2) - ord_mod_prt.ordinal(r1);
      double ms = ord_mod_sol.ordinal(r2) - ord_mod_sol.ordinal(r1);
      printf( "%7.3f %7.3f      %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f\n",
      	r1, r2, rs, wp, rr, mp, ws, ms );
    }
    printf("$$\n");
    std::cout << std::endl << std::endl;

    // apply gamma
    std::cout << "Gamma " << gamma << std::endl << std::endl;
    for ( MRI ix = msk_prt.first(); !ix.last(); ix.next() )
      map_mod[ix] -= gamma * map_wrk[ix];

    // back transform
    map_mod.fft_to( wrk_fp );

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

    std::cout << "Log likelihood:" << sfw.log_likelihood_free() << std::endl;
    printf("$TABLE :SigmaA statistics:\n");
    printf("$GRAPHS :SigmaA statistics:N:1,2,3: $$\n");
    printf(" 1/resol^2  sigmaA(s)  sigmaA(w) $$\n");
    printf("$$\n");
    int npweight = sfw.params_scale().size();
    clipper::Resolution_ordinal ord;
    ord.init( flagwt, cell_wrk, 1.0 ); ord.invert();
    for ( int i = 0; i < npweight; i++ ) {
      double err = sfw.params_error()[i];
      double s = ord.ordinal( (double(i)+0.5)/double(npweight) );
      double sigmaa1 = sfw.params_scale()[i];
      double sigmaa2 = ( err < 1.0 ) ? sqrt(1.0-err) : 0.0;
      printf( " %8.3f   %8.3f   %8.3f\n", s, sigmaa1, sigmaa2 );
    }
    printf("$$\n");
    std::cout << std::endl << std::endl;

  }  // cycle loop

  // output new results
  Compute_scale_u_aniso_fphi compute_aniso( 1.0, uaniso );
  wrk_fp.compute( wrk_fp, compute_aniso );
  if ( opcol_hl[0] != '/' ) opcol_hl = oppath + opcol_hl;
  if ( opcol_fc[0] != '/' ) opcol_fc = oppath + opcol_fc;
  mtzfile.open_append( ipmtz_wrk, opfile );
  mtzfile.export_hkl_data( mod_hl, opcol_hl );
  mtzfile.export_hkl_data( wrk_fp, opcol_fc );
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
		<< nxop.tgt().y() << "," << nxop.src().z()
		<< std::endl;
    }
  }
  std::cout << std::endl;

  prog.set_termination_message( "Normal termination" );
}
