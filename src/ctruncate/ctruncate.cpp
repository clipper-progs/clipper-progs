//
//     CTRUNCATE
//     Copyright (C) 2006-2017 Norman Stein, Charles Ballard
//
//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Application.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//

#include "clipper/clipper.h"
#include "clipper/clipper-contrib.h"
#include "clipper/clipper-ccp4.h"
#include "clipper/clipper-mmdb.h"
#include "clipper/core/clipper_util.h"
#include "clipper/core/atomsf.h"
#include "clipper/core/coords.h"
#include "ccp4/ccp4_general.h"
#include "ccp4/ccp4_program.h"
#include "clipper/ccp4/ccp4_mtz_io.h"
#include "ccp4/cmtzlib.h"
#include "ccp4/csymlib.h"
#include <iostream>
#include <math.h>
#include <algorithm>
#include <iomanip>
#include <ctime>
#include <fstream>

//#include "ccp4/ccp4_fortran.h"
#include "intensity_target.h"  // contains additions to resol_targetfn.h
#include "intensity_scale.h"   // contains additions to sfscale.cpp, sfscale.h, function_object_bases.h
#include "alt_hkl_datatypes.h"

#include "cpsf_utils.h"
#include "ctruncate_truncate.h"
#include "ctruncate_utils.h"
#include "ctruncate_twin.h"
#include "ctruncate_parity.h"
#include "ctruncate_moments.h"
#include "ctruncate_analyse.h"
#include "ctruncate_aniso.h"
#include "ctruncate_matthews.h"
#include "ctruncate_wilson.h"
#include "best.h"

#include "mmdb2/mmdb_tables.h"

#include "ccp4/ccp4_utils.h"

using namespace clipper;
using namespace ctruncate;

// replacement for Wilson/Truncate


int main(int argc, char **argv)
{
    clipper::String prog_string = "ctruncate";
    clipper::String prog_vers = "1.17.23";
    clipper::String prog_date = "$Date: 2017/03/02";
	ctruncate::CCP4Program prog( prog_string.c_str(), prog_vers.c_str(), prog_date.c_str() );
    
    // defaults
    clipper::String outfile = "ctruncate_out.mtz";
    clipper::String outcol = "";
    clipper::String freecol = "/*/*/[FreeR_flag]";
    clipper::String appendcol = "";
    clipper::String meancol = "NONE";
    clipper::String pluscol = "";
    clipper::String minuscol = "";
    clipper::String anocols = "NONE";
    clipper::String ipfile = "NONE";
    clipper::String twintest = "first_principles";
    clipper::String ipseq = "NONE";
    clipper::String composition = "protein";
    clipper::String prior_select = "auto";
    clipper::String xmlfile = "";
    std::vector<clipper::String> history;
    
    bool doaniso = true;
    bool debug = false;
    bool freein = false;
    bool amplitudes = false;
    bool outImean = false;
    bool anomalous = false;
    bool refl_mean = false;
    bool is_nucl = false;
    bool outxml = false;
    
    int mtzinarg = 0;
    int mtzoutarg = 0;
    int nbins = 60;
    int nresidues = 0;
    int nprm = 60;
    
    enum MODE {AUTO,WILSON,FLAT,SIVIA};
    MODE prior = AUTO;
	
    clipper::Resolution reso_Patt = clipper::Resolution( 4.0 );
    clipper::Resolution reso_trunc;
    
    clipper::Resolution reso_u1, reso_u2, reso_u3, reso_u;
    
    //time information
    time_t now = std::time(0);
    std::string date_time(ctime(&now) );
    date_time.replace(date_time.find('\n'),1,1,' ');
    
    // command input
    prog.summary_beg();
    printf("\nUSER SUPPLIED INPUT:\n");
    CCP4CommandInput args( argc, argv, true ); 
    prog.summary_end();
    int arg = 0;
    while ( ++arg < args.size() ) {
        if ( args[arg] == "-mtzin" || args[arg] == "-hklin") {
            if ( ++arg < args.size() ) {
                ipfile = args[arg];
                mtzinarg = arg;
            }
        } else if ( args[arg] == "-mtzout" || args[arg] == "-hklout") {
            if ( ++arg < args.size() ) {
                outfile = args[arg];
                mtzoutarg = arg;
            }
        } else if ( args[arg] == "-colin" ) {
            if ( ++arg < args.size() ) meancol = args[arg];
			refl_mean = true;
        } else if ( args[arg] == "-colplus" ) {
            if ( ++arg < args.size() ) pluscol = args[arg];
            anomalous = true;
            printf("obsolete argument - use -colano instead\n");
            printf("e.g. -colano /*/*/[I(+),SIGI(+),I(-),SIGI(-)]\n");
            return(0);
        } else if ( args[arg] == "-colminus" ) {
            if ( ++arg < args.size() ) minuscol = args[arg];
            anomalous = true;
            printf("obsolete argument - use -colano instead\n");
            printf("e.g. -colano /*/*/[I(+),SIGI(+),I(-),SIGI(-)]\n");
            return(0);
        } else if ( args[arg] == "-colano" ) {
            if ( ++arg < args.size() ) anocols = args[arg];
            anomalous = true;
        } else if ( args[arg] == "-freein" ) {
            freein = true;
            if ( ++arg < args.size() ) freecol = args[arg];
        } else if ( args[arg] == "-colout" ) {
            if ( ++arg < args.size() ) appendcol = args[arg];
        } else if ( args[arg] == "-nbins" ) {
            if ( ++arg < args.size() ) nbins = clipper::String(args[arg]).i();
        } else if ( args[arg] == "-nres" ) {
            if ( ++arg < args.size() ) nresidues = clipper::String(args[arg]).i();
        } else if ( args[arg] == "-twintest" ) {
            if ( ++arg < args.size() ) twintest = args[arg];
        } else if ( args[arg] == "-seqin" ) {
            if ( ++arg < args.size() ) ipseq = args[arg];
        } else if ( args[arg] == "-tNCSres" ) {
            if ( ++arg < args.size() ) reso_u1 = clipper::Resolution( clipper::String(args[arg]).f() );
        } else if ( args[arg] == "-twinres" ) {
            if ( ++arg < args.size() ) reso_u2 = clipper::Resolution( clipper::String(args[arg]).f() );
        } else if ( args[arg] == "-anisores" ) {
            if ( ++arg < args.size() ) reso_u3 = clipper::Resolution( clipper::String(args[arg]).f() );
        } else if ( args[arg] == "-reso" ) {
            if ( ++arg < args.size() ) reso_u = clipper::Resolution( clipper::String(args[arg]).f() );
        } else if ( args[arg] == "-no-aniso" ) {
            doaniso = false;
        } else if ( args[arg] == "-amplitudes" ) {
            amplitudes = true;
        } else if ( args[arg] == "-Imean" ) {
            outImean = true;
        } else if ( args[arg] == "-comp" ) {
			if ( ++arg < args.size() ) composition = args[arg];
		} else if ( args[arg] == "-prior" ) {
			if ( ++arg < args.size() ) prior_select = args[arg];
        } else if ( args[arg] == "-debug" ) {
            debug = true;
        } else if ( args[arg] == "-history" ) {
            if ( ++arg < args.size() ) history.push_back( args[arg] );
        } else if ( args[arg] == "-xmlout" ) {
            if ( ++arg < args.size() ) {
                xmlfile = args[arg];
                outxml = true;
            }
        } else if ( args[arg] == "-i" ) {
            CCP4::ccp4_prog_info();
            return(0);
        } else if ( args[arg] == "-help") {
            std::cout << std::endl << "Important input" << std::endl;
            std::cout << " -mtzin <filename> or -hklin <filename>" << std::endl;
            std::cout << "    input mtz file name " << std::endl;
            std::cout << " -mtzout <filename> or -hklout <filename>" << std::endl;
            std::cout << "    output mtz file name [default ctruncate_out.mtz]. " << std::endl;
            std::cout << " -colin <column labels>" << std::endl;
            std::cout << "    input column labels eg /*/*/[IMEAN,SIGIMEAN] or" << std::endl;
            std::cout << "    /*/*/[I(+),SIGI(+),I(-),SIGI(-)]" << std::endl;
            std::cout << " -colano <anomalous column labels>" << std::endl;
            std::cout << "    alternative to -colin for anomalous data" << std::endl;
            std::cout << " -amplitudes" << std::endl;
            std::cout << "    input data is in form of amplitudes, only run tests" << std::endl;
            std::cout << std::endl;
            std::cout << "Lesser input:" << std::endl;
            std::cout << " -colout <output column prelabel>" << std::endl;
            std::cout << "    [default: no prelabel]" << std::endl;
            std::cout << " -freein <free R column label>" << std::endl;
            std::cout << "    copy free R column to output file" << std::endl;
            std::cout << " -Imean" << std::endl;
            std::cout << "    output Imean as well as anomalous" << std::endl;
            std::cout << " -no-aniso" << std::endl;
            std::cout << "    do not apply anomalous correction" << std::endl;
            std::cout << " -history <history>" << std::endl;
            std::cout << "    add to mtz file history" << std::endl;
            std::cout << " -xmlout <filename>" << std::endl;
            std::cout << "    output ccp4i2 xml file [default program.xml]" << std::endl;
            return(0);
        } else {
            printf("Unrecognised argument\n");
            return(1);
        }
        
    }
	
	
	//is colin likeley to be anomalous or not (using , seperator)
	if (refl_mean) {
		int i(0),pos(0);
		while ( (pos = meancol.find(',',++pos)) != std::string::npos ) {
			++i;
		}
		anomalous = anomalous || (i == 3);
		refl_mean = (i == 1);
	}
	
    if (anomalous) {
        clipper::CCP4MTZ_type_registry::add_group( "G_sigG_ano", "FANO" );
        clipper::CCP4MTZ_type_registry::add_group( "J_sigJ_ano", "IANO" );
        clipper::CCP4MTZ_type_registry::add_type( "ISym", "Y", 1.0);
        clipper::CCP4MTZ_type_registry::add_group( "ISym", "ISYM" );
		if (anocols == "NONE") anocols = meancol;
    }
    if ( args.size() <= 1 || ( !refl_mean && !anomalous )) {
        CCP4::ccperror(1,"Usage: ctruncate -mtzin <filename>  -mtzout <filename>  -colin <colpath> -colano <colpath> ");
    }
	
	if ( composition == "nucleic" ) is_nucl = true;
    
	if ( prior_select == "wilson" ) prior = WILSON;
	else if ( prior_select == "flat" ) prior = FLAT;
	else if ( prior_select == "sivia" ) prior = SIVIA;
	
    if (mtzinarg == 0) CCP4::ccperror(1, "No input mtz file");
    
    typedef clipper::HKL_data_base::HKL_reference_index HRI;
    
	CCP4MTZfile mtzfile;
	HKL_info hklinf;
    mtzfile.open_read( ipfile );

    //mtzfile.import_hkl_info( hklinf,false ); // include missing reflections
    mtzfile.import_hkl_info( hklinf,true);
    // allocate memory to isig by reading in hklinf before declaring isig
	HKL_data<data32::I_sigI> isig_import(hklinf); // raw I and sigma
    HKL_data<data32::I_sigI> isig(hklinf);   // working I and sigma
    HKL_data<data32::I_sigI> jsig(hklinf);   // post-truncate I and sigma
    HKL_data<data32::F_sigF> fsig(hklinf);   // post-truncate F and sigma 
    HKL_data<data32::I_sigI_ano> isig_ano_import(hklinf);   // raw I(+) and sigma and I(-) and sigma
	HKL_data<data32::I_sigI_ano> isig_ano(hklinf); //working I(+) and sigma and I(-) and sigma
    HKL_data<data32::I_sigI_ano> jsig_ano(hklinf);   // post-truncate anomalous I and sigma
    HKL_data<data32::F_sigF_ano> fsig_ano(hklinf);   // post-truncate anomalous F and sigma
    HKL_data<data32::Flag> free(hklinf);
    
    clipper::HKL_data<clipper::data32::F_sigF> faniso( hklinf );
    
    try {
    if (amplitudes ) {
        if ( refl_mean )
			mtzfile.import_hkl_data( fsig, meancol );
		if (anomalous)
            mtzfile.import_hkl_data( fsig_ano, anocols );
    } else {
        if ( refl_mean )
            mtzfile.import_hkl_data( isig_import, meancol );
        if (anomalous)
            mtzfile.import_hkl_data( isig_ano_import, anocols );
    }
    } catch (...) {
        std::cout << std::endl;
        std::cout << "Error: ";
        if ( refl_mean ) {
            std::cout << meancol << " could not be loaded.";
        }
        if (anomalous) {
            std::cout << anocols << " could not be loaded.";
        }
        std::cout << std::endl;
        return(1);
    }
 
    ReflectionFile reflnfile(mtzfile,ipfile);
    try {
        if (freein) mtzfile.import_hkl_data( free, freecol );
    } catch (...) {
        std::cout << std::endl;
        std::cout << "Error: " << freecol << " could not be loaded." << std::endl;
        return(1);
    }
    
	clipper::String mcol = ( refl_mean ) ? meancol : anocols;
		
    // need this mumbo-jumbo in order to write to output file    
    if ( outcol[0] != '/' ) outcol = mtzfile.assigned_paths()[0].notail()+"/"+outcol;
    
    // hkl_list contains only those (h,k,l) for which at least data column is not NaN.
    // hkl_info contains all (h,k,l) out to the resolution limit regardless of whether there is any measured data.
    // will need hkl_list when writing to output file.
    clipper::Spacegroup spgr = mtzfile.spacegroup();
    char spgr_confidence = mtzfile.spacegroup_confidence();
    clipper::Cell      cell1 = mtzfile.cell();
    clipper::Resolution reso = mtzfile.resolution();
    
    // limit resolution for truncation, and hence output
    reso_trunc = clipper::Resolution(  reso.limit()  );
    
    HKL_info hkl_list;
    hkl_list.init( spgr, cell1, reso );
    mtzfile.import_hkl_list(hkl_list);
    
    MTZdataset cset;
    MTZcrystal cxtl; 
    std::vector<clipper::String> histin;
    mtzfile.import_crystal ( cxtl, mcol );
    mtzfile.import_dataset ( cset, mcol );
    histin = mtzfile.history();
    
    mtzfile.close_read();
	
	reflnfile.output();
	
	// reconstruct Imean if it was not given
    if ( refl_mean ) {
		if (amplitudes) {
			for ( HRI ih = fsig.first(); !ih.last(); ih.next() ) {
				if ( !fsig[ih].missing() )
					isig[ih] = clipper::data32::I_sigI( fsig[ih].f()*fsig[ih].f(), 2.0*fsig[ih].f()*fsig[ih].sigf() );
			}
		} else {
			for ( HRI ih = isig_import.first(); !ih.last(); ih.next() ) {
				isig[ih] = isig_import[ih];	
				if ( isig[ih].sigI() <= 0.0 ) isig[ih].I() = isig[ih].sigI() = clipper::Util::nan();
			}
		}
	}
    if (anomalous) {
        for ( HRI ih = isig_ano_import.first(); !ih.last(); ih.next() ) {
            isig_ano[ih] = isig_ano_import[ih];
            if (isig_ano[ih].sigI_mi() <= 0.0 )  isig_ano[ih].I_mi() = isig_ano[ih].sigI_mi() = clipper::Util::nan();
            if (isig_ano[ih].sigI_pl() <= 0.0 )  isig_ano[ih].I_pl() = isig_ano[ih].sigI_pl() = clipper::Util::nan();
        }
        if (!refl_mean) {
            if ( amplitudes) {
                for ( HRI ih = fsig.first(); !ih.last(); ih.next() ) {
                    //this is none weighted mean
                    isig[ih].I() = clipper::Util::mean(std::pow(fsig_ano[ih].f_pl(),2),std::pow(fsig_ano[ih].f_mi(),2));
                    isig[ih].sigI() = clipper::Util::sig_mean(2.0f*fsig_ano[ih].f_pl()*fsig_ano[ih].sigf_pl(),2.0f*fsig_ano[ih].f_mi()*fsig_ano[ih].sigf_mi(), 0.0f  );
                }
            } else {
                for ( HRI ih = isig_ano.first(); !ih.last(); ih.next() ) {
                    //this is none weighted mean
                    //isig[ih].I() = clipper::Util::mean(isig_ano[ih].I_pl(),isig_ano[ih].I_mi());
                    //isig[ih].sigI() = clipper::Util::sig_mean(isig_ano[ih].sigI_pl(),isig_ano[ih].sigI_mi(), 0.0f );
                    if ( ( clipper::Util::is_nan(isig_ano[ih].I_pl() ) || isig_ano[ih].sigI_pl() <= 0.0 ) &&
                        ( clipper::Util::is_nan(isig_ano[ih].I_mi() ) || isig_ano[ih].sigI_mi() <= 0.0 ) ) {
                        isig[ih].I() = isig[ih].sigI() = clipper::Util::nan();
                    } else if ( clipper::Util::is_nan(isig_ano[ih].I_pl() ) || isig_ano[ih].sigI_pl() <= 0.0 ) {
                        isig[ih].I() = isig_ano[ih].I_mi();
                        isig[ih].sigI() = isig_ano[ih].sigI_mi();
                    } else if ( clipper::Util::is_nan(isig_ano[ih].I_mi() ) || isig_ano[ih].sigI_mi() <= 0.0 ) {
                        isig[ih].I() = isig_ano[ih].I_pl();
                        isig[ih].sigI() = isig_ano[ih].sigI_pl();
                    } else {
                        isig[ih].I() = 0.5*(isig_ano[ih].I_pl()+isig_ano[ih].I_mi() );
                        isig[ih].sigI() = std::sqrt(0.5*(std::pow(isig_ano[ih].sigI_pl(),2)+std::pow(isig_ano[ih].sigI_mi(),2)) );
                    }
                }
			}
		}
	}
	
    
    int Ncentric = 0;
    int Nreflections = 0;
	
	ReflectionData reflndata;
    if ( refl_mean ) {
		if (amplitudes) {
			reflndata(fsig);
		} else {
			reflndata(isig);
		}
	} else {
		if ( amplitudes) {
			reflndata(fsig_ano);
		} else {
			reflndata(isig_ano_import);
		}
	}
	reflndata.output();
	
    //Cell contents analysis
	{
		if ( ipseq != "NONE" ) {
			prog.summary_beg(); printf("CELL CONTENTS:\n\n");
			
			clipper::SEQfile seqf;
			seqf.read_file( ipseq );
			
			ctruncate::Matthews cmath(true,false);
			int nmol = cmath(cell1, spgr, seqf, 1.0/sqrt(isig.invresolsq_range().max() ) );
			std::cout << "Expected number of molecules in ASU : " << nmol << std::endl;
			prog.summary_end();
			cmath.summary();
		} else if (nresidues > 0) {		
			prog.summary_beg(); printf("CELL CONTENTS:\n\n");
			ctruncate::Matthews cmath(true,false);
			int nmol = cmath(cell1, spgr, nresidues, 1.0/sqrt(isig.invresolsq_range().max() ) );
			std::cout << "Expected number of molecules in ASU : " << nmol << std::endl;
			prog.summary_end();
			cmath.summary();
		}
    }
	
	    clipper::Range<double> active_range;

	std::stringstream xml_comp;
	ctruncate::Rings icerings;
    clipper::ftype scalef=1.0;
	{
		HKLAnalysis hklanalysis(isig);
		hklanalysis.output();
		active_range = hklanalysis.active_range();
		icerings = hklanalysis.ice_rings();
        hklanalysis.xml_output(xml_comp);
        // if something went wrong with Wilson scaling, B could be negative, giving exponentially large scaled SF's
        // so only scale if B positive
        if ( hklanalysis.wilson_intercept() > 0 ) scalef = sqrt(hklanalysis.wilson_intercept() );
		//std::cout << (hklanalysis.xml_output(xml_comp)).str() << std::endl;
	}
	
    // anomalous signal
	std::stringstream xml_anomstats;
	if (anomalous ) {
		if (amplitudes) {
			AnomStats anomstats(fsig_ano);
			anomstats.output();
			/*if (outxml)*/ anomstats.xml_output(xml_anomstats);
		} else {
			AnomStats anomstats(isig_ano);
			anomstats.output();
			/*if (outxml)*/ anomstats.xml_output(xml_anomstats);
		}
	}
    
	std::stringstream xml_tncs;
	std::vector<clipper::Symop> ncs_list;
	bool hastncs(false);
	{
		// limit resolution of Patterson calculation for tNCS (default 4 A), or set to
		// limit from completeness analysis
		float invopt = active_range.max();
		float resopt = 1.0/std::sqrt(invopt);
		reso_Patt = clipper::Resolution( std::max( double(resopt), reso_Patt.limit() ) );
		//user override of defaults
		if (!reso_u1.is_null() ) {
			reso_Patt = clipper::Resolution( clipper::Util::max( reso.limit(), reso_u1.limit() ) );
		}
		
        
		// check for pseudo translation (taken from cpatterson)
		ctruncate::tNCS tncs;
		ncs_list = tncs(isig, reso_Patt );
		hastncs = tncs.hasNCS();
		
		prog.summary_beg();
		tncs.output();
		prog.summary_end();
		printf("\n");
		tncs.xml_output(xml_tncs);
		//std::cout << xml_tncs.str() << std::endl;
	}
	
	//setup aniso copy of isig
	HKL_data<data32::I_sigI> ianiso(hklinf);
	for ( HRI ih = ianiso.first(); !ih.last(); ih.next() ) { 
		float I = isig[ih.hkl()].I();
		float sigI = isig[ih.hkl()].sigI();
		ianiso[ih] = clipper::data32::I_sigI( I, sigI );
	}
	
    // anisotropy estimation
	std::stringstream xml_aniso;
    //clipper::U_aniso_orth uao;
    clipper::U_aniso_orth uaoc(0,0,0,0,0,0);
	bool anisobysymm(false);
	bool anisodemo(false);
    {
		clipper::Range<clipper::ftype> range_Aniso(active_range.min(), (!reso_u3.is_null() ) ?  std::min(active_range.max(),1.0/std::pow(reso_u3.limit(),2) ) : active_range.max() );
		AnisoAnalysis aa(isig,range_Aniso);
		anisobysymm=aa.allowed_by_symmetry();
		anisodemo=aa.is_anisotropic();
		aa.output();
		aa.xml_output(xml_aniso);
		//std::cout << xml_aniso.str() << std::endl;
		
		//set ianiso
        if (doaniso && anisobysymm) {
		    uaoc=aa.u_aniso_orth_corr_F();
		    clipper::datatypes::Compute_scale_u_aniso<clipper::data32::I_sigI > compute_s(1.0,aa.u_aniso_orth_corr_F() );
		    ianiso.compute(ianiso, compute_s);
        }
	}
	
	//want to use anisotropy correction and resolution truncation for twinning tests
	
	bool hastwin(false);
    std::stringstream xml_twin;
	{
		clipper::Range<clipper::ftype> range_Twin(active_range.min(), (!reso_u2.is_null() ) ?  std::min(active_range.max(),1.0/std::pow(reso_u2.limit(),2) ) : active_range.max() );
		std::vector<clipper::Symop> tt(ncs_list.size() );
		for (int i= 0; i != tt.size() ; ++i ) tt[i] = ncs_list[i];
		TwinAnalysis twins(ianiso,tt,range_Twin);
		twins.output();
		//std::stringstream xt;
		//std::cout << (twins.xml_output(xml_twin) ).str() << std::endl;
		twins.xml_output(xml_twin);
		//Parity group analysis
		ctruncate::parity(ianiso, active_range.max(), nbins);
    }
		
	HKL_data<data32::I_sigI> xsig(hklinf);
	//normal calculation
    if (!reso_u.is_null() ) {
        reso_trunc =
        clipper::Resolution( clipper::Util::max( reso_trunc.limit(), reso_u.limit() ) );
    }
    std::stringstream xml_trunc;
	if (!amplitudes) {
        MODE prior_user = prior;
        //user override of truncate procedure and output
        if ( prior == AUTO && (hastncs || hastwin ) ) prior = FLAT;
        
        int nrej_ice(0), nrej_norm(0), nrej_pre(0);
        bool ierror(false);
        double invresolsq(reso_trunc.invresolsq_limit() );
        
        for ( HRI ih = xsig.first(); !ih.last(); ih.next() ) {
            double reso = ih.invresolsq();
            xsig[ih] = clipper::data32::I_sigI( (isig[ih.hkl()].I()), isig[ih.hkl()].sigI() );
            if ( icerings.InRing(reso) != -1 )
                if ( icerings.Reject( icerings.InRing(reso) ) ) {
                    xsig[ih].I() = xsig[ih].sigI() = clipper::Util::nan(); // loose ice rings
                    ++nrej_ice;
                }
            if (reso > invresolsq )
                xsig[ih].I() = xsig[ih].sigI() = clipper::Util::nan();
        }
        
        int nprm2=100;
        int nreflns=200;
        
        int Nreflections = 0;
        {
            for ( HRI ih = xsig.first(); !ih.last(); ih.next() ) {
                clipper::ftype reso = ih.invresolsq();
                if ( !xsig[ih].missing() ) ++Nreflections;
            }
            if ( nprm2 == 0 && nreflns != 0 ) {
                nprm2 = std::max( Nreflections/nreflns , 1);
                //} else if ( nreflns == 0 && nprm2 != 0 ) {
                //nprm = nbins;
            } else {
                //nprm2 = std::max( Nreflections/nreflns , nprm2);
                double np1(nprm2+0.499);
                double np2(Nreflections/nreflns);
                double np(std::sqrt(np1*np1*np2*np2/(np1*np1+np2*np2) ) );
                nprm2 = std::max( int(np), 1 );
            }
        }
        
        HKL_data<data32::I_sigI> tr1(hklinf);
        
        // calc scale
        for ( HRI ih = tr1.first(); !ih.last(); ih.next() ) {
            double reso = ih.invresolsq();
            if ( reso > ctruncate::Best::invresolsq_max() ) reso =  ctruncate::Best::invresolsq_max();
            if ( reso < ctruncate::Best::invresolsq_min() ) reso =  ctruncate::Best::invresolsq_min();
            tr1[ih] = clipper::data32::I_sigI( xsig[ih.hkl()].I()/ctruncate::Best::value(reso), 0.0);
        }
        
        std::vector<double> params(nprm2,1.0);
        //precondition the ML calc using least squares fit.  This should give an excellent start
        clipper::BasisFn_binner basis_pre( tr1, nprm2, 1.0 );
        TargetFn_meanInth<clipper::data32::I_sigI> target_pre(tr1,1.0);
        clipper::ResolutionFn pre( hklinf, basis_pre, target_pre, params);
        params = pre.params();
        
        for (int i=0; i != params.size() ; ++i) {
            if (params[i] < 0.0) ierror = true;
        }
        
        //outlier rejection from Read (1999)
        //double rlimit(0.0);
        for ( HRI ih = tr1.first(); !ih.last(); ih.next() ) {
            if ( !xsig[ih.hkl()].missing() ) {
                double I = tr1[ih].I()/(ih.hkl_class().epsilon()*pre.f(ih));
                double rlimit = (ih.hkl_class().centric() ) ? 6.40*6.40 : 4.55*4.55 ;
                if (I > rlimit )  {
                    ++nrej_pre;
                    xsig[ih.hkl()].I() = xsig[ih.hkl()].sigI() = clipper::Util::nan();
                }
            }
        }
        
        if (ierror) {
            //back to least squares fit, but using spline
            std::vector<bool> mask(nprm2,false);
            clipper::BasisFn_spline basis_fo( tr1, nprm2, 1.0 );
            TargetFn_meanInth<clipper::data32::I_sigI> target_fo(tr1,1.0);
            clipper::ResolutionFn Sigma( hklinf, basis_fo, target_fo, params);
            params = Sigma.params();
            
            for ( HRI ih = xsig.first(); !ih.last(); ih.next() ) {
                double reso = ih.invresolsq();
                if ( reso > ctruncate::Best::invresolsq_max() ) reso =  ctruncate::Best::invresolsq_max();
                if ( reso < ctruncate::Best::invresolsq_min() ) reso =  ctruncate::Best::invresolsq_min();
                xsig[ih].I() = Sigma.f(ih) * ctruncate::Best::value(reso);
                xsig[ih].sigI() = 1.0;
            }
            
        } else {
            //reset bins as are scaling exponential
            int nprm2=20;
            int nreflns=500;
            
            //int Nreflections = 0;
            {
                //for ( HRI ih = xsig.first(); !ih.last(); ih.next() ) {
                //    clipper::ftype reso = ih.invresolsq();
                //    if ( !xsig[ih].missing() ) ++Nreflections;
                //}
                if ( nprm2 == 0 && nreflns != 0 ) {
                    nprm2 = std::max( Nreflections/nreflns , 1);
                    //} else if ( nreflns == 0 && nprm2 != 0 ) {
                    //nprm = nbins;
                } else {
                    //nprm2 = std::max( Nreflections/nreflns , nprm2);
                    double np1(nprm2+0.499);
                    double np2(Nreflections/nreflns);
                    double np(std::sqrt(np1*np1*np2*np2/(np1*np1+np2*np2) ) );
                    nprm2 = std::max( int(np), 1 );
                }
            }
            std::vector<double> params(nprm2,1.0);
            
            //reset tr1
            double reso(0.0);
            for ( HRI ih = tr1.first(); !ih.last(); ih.next() ) {
                reso = ih.invresolsq();
                if ( reso > ctruncate::Best::invresolsq_max() ) reso =  ctruncate::Best::invresolsq_max();
                if ( reso < ctruncate::Best::invresolsq_min() ) reso =  ctruncate::Best::invresolsq_min();
                tr1[ih] = clipper::data32::I_sigI( ctruncate::Best::value(reso), 0.0);
            }
            
            // proeceed with ML calc
            std::vector<bool> mask(nprm2,false);
            clipper::BasisFn_spline basis_fo( xsig, nprm2, 1.0 );
            TargetFn_scaleLogLikeI1I2<clipper::data32::I_sigI,clipper::data32::I_sigI> target_fo(tr1,xsig);
            ctruncate::ResolutionFn_nonlinear Sigma( hklinf, basis_fo, target_fo, params, mask, 1.0, false);
            params = Sigma.params();
            
            //generate norm
            for ( HRI ih = xsig.first(); !ih.last(); ih.next() ) {
                double reso = ih.invresolsq();
                if ( reso > ctruncate::Best::invresolsq_max() ) reso =  ctruncate::Best::invresolsq_max();
                if ( reso < ctruncate::Best::invresolsq_min() ) reso =  ctruncate::Best::invresolsq_min();
                xsig[ih] = clipper::data32::I_sigI(Sigma.f(ih) * ctruncate::Best::value(reso),1.0f);
            }
        }

		// scale the norm for the anisotropy
		if (doaniso) {
			if (anisobysymm && anisodemo) {
				clipper::datatypes::Compute_scale_u_aniso<clipper::data32::I_sigI > compute_s(1.0,-uaoc);
				xsig.compute(xsig, compute_s);
			}
		}
		
		if (!ierror) {
			//outlier rejection from Read (1999) using new norm
			double rlimit(0.0);
			for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
                double reso = ih.invresolsq();
                double I(0.0);
                if (reso <= invresolsq ) {
                    if ( !isig[ih].missing() ) {
                        int ring = icerings.InRing(reso);
                        if ( ring == -1 ) I = isig[ih].I()/(ih.hkl_class().epsilon()*xsig[ih.hkl()].I() );
                        else if ( !icerings.Reject( ring) ) I = isig[ih].I()/(ih.hkl_class().epsilon()*xsig[ih.hkl()].I() );
					    rlimit = (ih.hkl_class().centric() ) ? 6.40*6.40 : 4.55*4.55 ;
					    if (I > rlimit )  {
						    ++nrej_norm;
                        }
					}
				}
			}	
        }
        
        int nrej(0);
        if ( refl_mean ) {
            if (prior == FLAT ) truncate( isig, jsig, fsig, scalef, reso_trunc, nrej, debug );
            else if (prior == SIVIA) truncate_sivia(isig, jsig, fsig, scalef, reso_trunc, nrej, debug );
            else truncate( isig, jsig, fsig, xsig, scalef, reso_trunc, nrej, debug );
        }
        if (anomalous) {
            if (prior == FLAT ) truncate( isig_ano, jsig_ano, fsig_ano, scalef, reso_trunc, nrej, debug );
            else if (prior == SIVIA) truncate_sivia( isig_ano, jsig_ano, fsig_ano, scalef, reso_trunc, nrej, debug );
            else truncate( isig_ano, jsig_ano, fsig_ano, xsig, scalef, reso_trunc, nrej, debug );
            int iwarn = 0;
            // use for phil plot
            if (!refl_mean ) {
                for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
                    jsig[ih].I() = jsig_ano[ih].I();
                    jsig[ih].sigI() = jsig_ano[ih].sigI();
                    fsig[ih].f() = fsig_ano[ih].f();
                    fsig[ih].sigf() = fsig_ano[ih].sigf();
                }
            }
        }
        
        std::cout << std::endl;
        prog.summary_beg();
        std::cout << std::endl << "INTENSITY TO AMPLITUDE CONVERSION:" << std::endl << std::endl;
        if (prior == prior_user) {
            if ( prior == FLAT ) std::cout <<  "Calculation using FLAT prior." << std::endl;
            else if ( prior == SIVIA ) std::cout << "Calculation using SIVIA method." << std::endl;
            else std::cout << "Calculation using Wilson prior." << std::endl;
        } else {
            std::cout << "WARNING: FLAT prior in use due to either tNCS or twinning.\nTo override force --prior WILSON." << std::endl;
        }
        printf("During the truncate procedure %d intensities have been flagged as unphysical.\n\n", nrej);
        prog.summary_end();
        std::cout << std::endl;
        std::cout << "Norm calculation summary:" << std::endl << std::endl;
        if (ierror) std::cout << "      WARNING: negative mean I in bins, resorted to least squares norm." << std::endl;
        std::cout << "      Number of outliers and ice ring reflections not used in norm calculation (Read (1999) ): " << nrej_pre+nrej_ice << std::endl;
        if (!ierror) {
            std::cout << "      Number of outliers in detected in final norm (Read (1999) ): " << nrej_norm << std::endl;
            if (nrej_norm > nrej_pre && ( prior != FLAT || prior != SIVIA ) ) std::cout << "      WARNING: prior may be unstable, change in rejected reflections" << std::endl ;
        }
        if (doaniso) std::cout << "      Anisotropy correction applied to norm." << std::endl;
        std::cout << std::endl << std::endl;
        
        xml_trunc << "<Truncation>" << std::endl;
        xml_trunc << "  <prior>";
        if ( prior == FLAT ) xml_trunc <<  "FLAT";
        else if ( prior == SIVIA ) xml_trunc << "SIVIA";
        else xml_trunc << "WILSON";
        xml_trunc << "</prior>" << std::endl;
        xml_trunc << "  <mode>";
        if ( prior == AUTO ) xml_trunc << "AUTO";
            else if ( prior == FLAT ) xml_trunc <<  "FLAT";
            else if ( prior == SIVIA ) xml_trunc << "SIVIA";
            else xml_trunc << "WILSON";
        xml_trunc << "</mode>"  << std::endl;
        xml_trunc << "  <AnisoCorrection>" << ((doaniso) ? "yes" : "no" ) << "</AnisoCorrection>" << std::endl;
        xml_trunc << "  <ResolutionRange id=\"Truncation\" unit=\"Angstrom\" >" << std::endl;
        xml_trunc << std::fixed << std::setprecision(2) << "    <min>" << 1.0/std::sqrt(isig.invresolsq_range().min() ) << "</min>\n    <max>"
        << std::fixed << std::setprecision(2) <<  1.0/std::sqrt(invresolsq ) << "</max>" << std::endl;
        xml_trunc << "  </ResolutionRange>" << std::endl;
        if (prior != prior_user) xml_trunc << "  <Warning id=\"flat\">FLAT prior in use due to either tNCS or twinning. To override force --prior WILSON.</Warning>" << std::endl;
        if (ierror) xml_trunc << "  <Warning id=\"negative\">Negative mean I in bins, resorted to least squares norm.</Warning>" << std::endl;
        xml_trunc << "</Truncation>" << std::endl;
	}

	{
		ctruncate::PattPeak patt_peak(std::sqrt(isig.invresolsq_range().max() ));
		
		float opt_res = patt_peak(isig);
		
		float width_patt = 2.0f*patt_peak.sigma();
		
		float b_patt = 4.0f*clipper::Util::twopi2()*std::pow(width_patt/2.0f,2.0f);
		
        std::cout << std::endl;
		prog.summary_beg();
		std::cout << "Estimated Optical Resolution: " << opt_res << std::endl;
		prog.summary_end();
        std::cout << std::endl << std::endl;
		
		
    }
    
    // I/Sigma and F/sqrt(Sigma) plots
    if (!amplitudes) {
        //could use a clipper::Histogram
        clipper::Range<double> range(-5.0,10.0);
        /*for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
         float eps = ih.hkl_class().epsilon();
         if (!isig[ih].missing() ) range.include(isig[ih].I()/(xsig[ih].I()));
         //if (!fsig[ih].missing() ) range.include(std::sqrt(eps)*fsig[ih].f()/std::sqrt(xsig[ih].I()) );
         }*/
        clipper::Histogram Icount(range,200);
        clipper::Histogram Jcount(range,200);
        clipper::Histogram Fcount(range,200);
        clipper::Histogram IScount(range,200);
        clipper::Histogram JScount(range,200);
        clipper::Histogram FScount(range,200);
        
        for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
            //float eps = ih.hkl_class().epsilon();
            if (!isig[ih].missing() ) {
                Icount.accumulate(isig[ih].I()/xsig[ih].I());
            }
        }
        
        for ( HRI ih = jsig.first(); !ih.last(); ih.next() ) {
            //float eps = ih.hkl_class().epsilon();
            if (!jsig[ih].missing() ) {
                Jcount.accumulate(jsig[ih].I()/xsig[ih].I());
            }
        }	
        
        for ( HRI ih = fsig.first(); !ih.last(); ih.next() ) {
            //float eps = ih.hkl_class().epsilon();
            if (!fsig[ih].missing() ) {
                Fcount.accumulate(fsig[ih].f()/(scalef*std::sqrt(xsig[ih].I()) ) );
            }
        }	
        
        
        for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
            //float eps = ih.hkl_class().epsilon();
            if (!isig[ih].missing() ) {
                IScount.accumulate(isig[ih].I()/isig[ih].sigI() );
            }
        }
        
        for ( HRI ih = jsig.first(); !ih.last(); ih.next() ) {
            //float eps = ih.hkl_class().epsilon();
            if (!jsig[ih].missing() ) {
                JScount.accumulate(jsig[ih].I()/jsig[ih].sigI() );
            }
        }	
        
        for ( HRI ih = fsig.first(); !ih.last(); ih.next() ) {
            //float eps = ih.hkl_class().epsilon();
            if (!fsig[ih].missing() ) {
                FScount.accumulate(fsig[ih].f()/fsig[ih].sigf() );
            }
        }	
        
        printf("$TABLE: Phil plot:\n");
        printf("$GRAPHS");
        printf(": Phil plot - normalised values:A:1,2,3,4:\n");  
        printf(": Phil plot - vs sigma:A:1,5,6,7:\n$$");  
        
        printf(" Value Io/Sigma I/Sigma F/Sigma**0.5 Io/sigIo I/sigI F/sigF$$\n$$\n");
        
        for ( int i=0; i!=200; ++i ) {
            float res = range.min()+float(i)*(range.max()-range.min())/200.0; 
            printf("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", 
                   res,
                   Icount.y(res),
                   Jcount.y(res),
                   Fcount.y(res),
                   IScount.y(res),
                   JScount.y(res),
                   FScount.y(res)
                   );
        }
        
        printf("$$\n\n");
    }		
    
    // output data
    if (!amplitudes) {
		CCP4MTZfile mtzout;
		HKL_data<data32::D_sigD> Dano(hklinf);   // anomalous difference and sigma 
		HKL_data<data32::ISym> freidal_sym(hklinf);
		HKL_data<data32::J_sigJ_ano> isig_ano_export(hklinf); // do not want to output cov term, so use 4 term ano
		HKL_data<data32::G_sigG_ano> fsig_ano_export(hklinf); // do not want to output cov term
		if (anomalous) {
			for ( HRI ih = freidal_sym.first(); !ih.last(); ih.next() ) {
				freidal_sym[ih].isym() = ( !Util::is_nan(fsig_ano[ih.hkl()].f_pl() )  &&  !Util::is_nan(fsig_ano[ih.hkl()].f_mi() ) ) ? 0 :
				( !Util::is_nan(fsig_ano[ih.hkl()].f_pl() ) ) ? 1 :
				( !Util::is_nan(fsig_ano[ih.hkl()].f_mi() ) ) ? 2 : 0;
            }
            for ( HRI ih = Dano.first(); !ih.last(); ih.next() ) {
                Dano[ih].d() = ( !Util::is_nan(fsig_ano[ih.hkl()].f_pl() )  &&  !Util::is_nan(fsig_ano[ih.hkl()].f_mi() ) ) ? (fsig_ano[ih.hkl()].f_pl() - fsig_ano[ih.hkl()].f_mi()) : clipper::Util::nan();
				Dano[ih].sigd() = ( !Util::is_nan(fsig_ano[ih.hkl()].f_pl() )  &&  !Util::is_nan(fsig_ano[ih.hkl()].f_mi() ) ) ? std::sqrt(fsig_ano[ih.hkl()].sigf_pl()*fsig_ano[ih.hkl()].sigf_pl()+fsig_ano[ih.hkl()].sigf_mi()*fsig_ano[ih.hkl()].sigf_mi() ) : clipper::Util::nan();
				if ( ih.hkl_class().centric() ) {
					Dano[ih].d() = ( !Util::is_nan(fsig_ano[ih.hkl()].f_pl() )  ||  !Util::is_nan(fsig_ano[ih.hkl()].f_mi() ) ) ? 0.0 : clipper::Util::nan();
					Dano[ih].sigd() = ( !Util::is_nan(fsig_ano[ih.hkl()].f_pl() )  ||  !Util::is_nan(fsig_ano[ih.hkl()].f_mi() ) ) ? 0.0 : clipper::Util::nan();
				}
			}
			for ( HRI ih = isig_ano_export.first(); !ih.last(); ih.next() ) {
				isig_ano_export[ih] = clipper::data32::J_sigJ_ano(isig_ano_import[ih.hkl()].I_pl(), isig_ano_import[ih.hkl()].I_mi(), isig_ano_import[ih.hkl()].sigI_pl(), isig_ano_import[ih.hkl()].sigI_mi() );
			}
			for ( HRI ih = fsig_ano_export.first(); !ih.last(); ih.next() ) {
				fsig_ano_export[ih] = clipper::data32::G_sigG_ano(fsig_ano[ih.hkl()].f_pl(), fsig_ano[ih.hkl()].f_mi(), fsig_ano[ih.hkl()].sigf_pl(), fsig_ano[ih.hkl()].sigf_mi() );
			}
		}			
		
        //mtzout.open_append( args[mtzinarg], outfile );
        mtzout.open_write( outfile );
        mtzout.export_crystal ( cxtl, outcol );
        mtzout.export_dataset ( cset, outcol );
        //mtzout.export_hkl_info( hklinf );
        mtzout.export_hkl_info( hkl_list );
        //mtzout.export_hkl_data( jsig, outcol );
        //clipper::String labels;
		if ( refl_mean ) {
            clipper::String labels;
			if (appendcol == "") labels = outcol + "[F,SIGF]";
			else labels = outcol + "[F_" + appendcol + ",SIGF_" + appendcol + "]";
			mtzout.export_hkl_data( fsig, labels );
        }
        if (freein) {
            if (appendcol != "") {
                String::size_type loc = freecol.find("]",0);
                freecol.insert(loc,"_"+appendcol);
            }
            mtzout.export_hkl_data( free, outcol + freecol.tail() );
        }
		
        if (anomalous) {
            clipper::String labels;
            if ( !refl_mean ) {
                if (appendcol == "") labels = "[FMEAN,SIGFMEAN]";
                else labels = "[FMEAN_" + appendcol + ",SIGFMEAN_" + appendcol + "]";
                mtzout.export_hkl_data( fsig, outcol+labels.tail() );
                labels.clear();
            }
            if (appendcol == "") labels = "[DANO,SIGDANO]";
            else labels = "[DANO_" + appendcol + ",SIGDANO_" + appendcol + "]";
            mtzout.export_hkl_data( Dano, outcol+labels.tail() );
            labels.clear();
			if (appendcol == "") labels = "[F(+),SIGF(+),F(-),SIGF(-)]";
            else labels = "[F_" + appendcol + "(+),SIGF_" + appendcol + "(+),F_" + appendcol + "(-),SIGF_" + appendcol + "(-)]";
            mtzout.export_hkl_data( fsig_ano_export, outcol+labels.tail() );
            labels.clear();
            if (appendcol == "") labels = "[ISYM]";
            else labels = "[ISYM_" + appendcol + "]";
            mtzout.export_hkl_data( freidal_sym, outcol+labels.tail() );
        }
		
		//output original input
		if ( refl_mean ) {
            clipper::String labels;
            labels = meancol;
			if (appendcol != "") {
				String::size_type loc = labels.find(",",0);
				labels.insert(loc,"_"+appendcol);
				loc = labels.find("]",0);
				labels.insert(loc,"_"+appendcol);
			}
			mtzout.export_hkl_data( isig_import, outcol + labels.tail() );
            labels.clear();
        }

		// KDC hack to output Imean  if explicitly requested (should be combined with above)
		if ( outImean ) mtzout.export_hkl_data( isig, outcol + appendcol + "_MEAN" );

        if (anomalous) {
            clipper::String labels;
            labels = anocols;
            if (appendcol != "") {
                String::size_type loc = labels.find("-",0);
                labels.insert(loc-1,"_"+appendcol);
                loc = labels.find(",",loc);
                loc = labels.find("-",loc);
                //loc = labels.find("-",loc+appendcol.size()+2 );
                labels.insert(loc-1,"_"+appendcol);
                loc = labels.find("+",0);
                labels.insert(loc-1,"_"+appendcol);
                loc = labels.find(",",loc);
                loc = labels.find("+",loc);
                //loc = labels.find("+",loc+appendcol.size()+2 );
                labels.insert(loc-1,"_"+appendcol);
            }
            mtzout.export_hkl_data( isig_ano_export, outcol + labels.tail() );
        }
		
        //copy old history and say something about ctruncate run
        if (history.size() != 0 ) {
			for (int i = 0 ; i != history.size() ; ++i ) histin.push_back(history[i]);
        }
        //char run_date[10];
        //CCP4::ccp4_utils_date(run_date);
        //char run_time[8];
        //CCP4::ccp4_utils_time(run_time);
        clipper::String run_type = ( prior == FLAT ) ? " flat " : " french-wilson ";
        
        //clipper::String chistory = prog_string + " " + prog_vers + run_type +  "run on " + run_date + " " + run_time;
        clipper::String chistory = prog_string + " " + prog_vers + run_type + "run on " + date_time;
        histin.push_back(chistory);
    
        mtzout.set_history(histin);
		
        mtzout.set_spacegroup_confidence(spgr_confidence);
        
        //mtzout.close_append();
        mtzout.close_write();
        
        // Clipper will change H3 to R3, so change it back
        if ((spgr.symbol_hm())[0] == 'R') {			
			CMtz::MTZ *mtz=NULL;
			int read_refs=1;  // need to read in reflections, otherwise they won't be written out
			mtz = CMtz::MtzGet(outfile.c_str(), read_refs);
			// write title to output file
			char title[72];
			CMtz::ccp4_lrtitl(mtz, title);
			strncpy( mtz->title, title, 71 );
			//reset spacegroup
			char spacegroup[20];
			CSym::CCP4SPG *spg = CSym::ccp4spg_load_by_ccp4_num(CMtz::MtzSpacegroupNumber(mtz));
			strcpy(spacegroup,spg->symbol_old);
			strcpy(mtz->mtzsymm.spcgrpname,spacegroup);
			CMtz::MtzPut( mtz, outfile.c_str() );
			CMtz::MtzFree( mtz );
		}
    }
    
    if (outxml) {
        time_t now = std::time(0);
		std::stringstream ss1,ss2,ss3;
        std::ofstream xf;
        xf.open(xmlfile.c_str() );
        xf << prog.xml_start(ss1).str();
        xf << xml_trunc.str();
		reflndata.xml_output(ss3);
		xf << reflnfile.xml_output(ss3).str();
		xf << xml_comp.str();
        xf << xml_anomstats.str();
        xf << xml_tncs.str();
		xf << xml_aniso.str();
        xf << xml_twin.str();
        xf << prog.xml_end(ss2).str();
    }
    
    prog.set_termination_message( "Normal termination" );
    
    return(0);
}





