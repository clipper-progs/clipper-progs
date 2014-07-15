//
//     CTRUNCATE
//     Copyright (C) 2006-2013 Norman Stein, Charles Ballard
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
#include "ctruncate_matthews.h"
#include "ctruncate_wilson.h"
#include "best.h"

#include "mmdb/mmdb_tables.h"

#include "ccp4/ccp4_utils.h"

using namespace clipper;
using namespace ctruncate;

// replacement for Wilson/Truncate


int main(int argc, char **argv)
{
    clipper::String prog_string = "ctruncate";
    clipper::String prog_vers = "1.15.9";
    clipper::String prog_date = "$Date: 2014/07/15";
    CCP4Program prog( prog_string.c_str(), prog_vers.c_str(), prog_date.c_str() );
    
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
    std::vector<clipper::String> history;
    
    bool aniso = true;
    bool debug = false;
    bool freein = false;
    bool amplitudes = false;
    bool anomalous = false;
    bool refl_mean = false;
    bool is_nucl = false;
    
    int mtzinarg = 0;
    int mtzoutarg = 0;
    int nbins = 60;
    int nresidues = 0;
    int nprm = 60;
    
    enum MODE {AUTO,WILSON,FLAT,SIVIA};
    MODE prior = AUTO;
	
    clipper::Resolution reso_Patt = clipper::Resolution( 4.0 );
    clipper::Resolution reso_trunc;
    
    clipper::Resolution reso_u1, reso_u2, reso_u3;
    
    // clipper seems to use its own column labels, then append yours
    
    CCP4MTZfile mtzfile, mtzout;
    HKL_info hklinf;
    
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
        } else if ( args[arg] == "-reso" ) {
            if ( ++arg < args.size() ) reso_u3 = clipper::Resolution( clipper::String(args[arg]).f() );
        } else if ( args[arg] == "-no-aniso" ) {
            aniso = false;
        } else if ( args[arg] == "-amplitudes" ) {
            amplitudes = true;
        } else if ( args[arg] == "-comp" ) {
			if ( ++arg < args.size() ) composition = args[arg];
		} else if ( args[arg] == "-prior" ) {
			if ( ++arg < args.size() ) prior_select = args[arg];
        } else if ( args[arg] == "-debug" ) {
            debug = true;
        } else if ( args[arg] == "-history" ) {
            if ( ++arg < args.size() ) history.push_back( args[arg] );
        } else if ( args[arg] == "-i" ) {
            CCP4::ccp4_prog_info();
            return(0);
        } else {
            printf("Unrecognised argument\n");
            return(0);
        }
        
    }
    if (anomalous) {
        clipper::CCP4MTZ_type_registry::add_group( "G_sigG_ano", "FANO" );
        clipper::CCP4MTZ_type_registry::add_group( "J_sigJ_ano", "IANO" );
        clipper::CCP4MTZ_type_registry::add_type( "ISym", "Y", 1.0);
        clipper::CCP4MTZ_type_registry::add_group( "ISym", "ISYM" );
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
    
    //mtzfile.open_read( args[1] );
    mtzfile.open_read( ipfile );
    mtzfile.import_hkl_info( hklinf );  
    // allocate memory to isig by reading in hklinf before declaring isig
    HKL_data<data32::I_sigI> isig(hklinf);   // raw I and sigma
    HKL_data<data32::I_sigI> jsig(hklinf);   // post-truncate I and sigma
    HKL_data<data32::F_sigF> fsig(hklinf);   // post-truncate F and sigma 
    HKL_data<data32::J_sigJ_ano> isig_ano_import(hklinf);   // raw I(+) and sigma and I(-) and sigma
    HKL_data<data32::J_sigJ_ano> jsig_ano(hklinf);   // post-truncate anomalous I and sigma
    HKL_data<data32::G_sigG_ano> fsig_ano(hklinf);   // post-truncate anomalous F and sigma 
    HKL_data<data32::D_sigD> Dano(hklinf);   // anomalous difference and sigma 
    HKL_data<data32::ISym> freidal_sym(hklinf);
    HKL_data<data32::Flag> free(hklinf);
    
    clipper::HKL_data<clipper::data32::F_sigF> faniso( hklinf );
    
    if (amplitudes ) {
        if ( refl_mean ) {
			mtzfile.import_hkl_data( fsig, meancol );
		}
		if (anomalous) {
            mtzfile.import_hkl_data( fsig_ano, anocols );
        }		
    } else {
        if ( refl_mean ) mtzfile.import_hkl_data( isig, meancol );
        
        if (anomalous) {
            mtzfile.import_hkl_data( isig_ano_import, anocols );
        }
    }
    if (freein) mtzfile.import_hkl_data( free, freecol );
    
	clipper::String mcol = ( refl_mean ) ? meancol : anocols;
	
    prog.summary_beg();
    printf("\nCRYSTAL INFO:\n\n");
    std::cout << "Crystal/dataset names: " << mtzfile.assigned_paths()[0].notail() << "\n"; 
    /*std::cout << mtzfile.assigned_paths()[0].notail().notail().tail() << "\n";  // crystal name
     std::cout << mtzfile.assigned_paths()[0].notail().tail() << "\n";  //dataset name
     String xtlname = mtzfile.assigned_paths()[0].notail().notail().tail();
     String setname = mtzfile.assigned_paths()[0].notail().tail();*/
    prog.summary_end();
    printf("\n");
    
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
	
	// reconstruct Imean if it was not given
    if ( refl_mean ) {
		if (amplitudes) {
			for ( HRI ih = fsig.first(); !ih.last(); ih.next() ) {
				if ( !fsig[ih].missing() )
					isig[ih] = clipper::data32::I_sigI( fsig[ih].f()*fsig[ih].f(), 2.0*fsig[ih].f()*fsig[ih].sigf() );
			}
		}
	} else {
		if ( amplitudes) {
			for ( HRI ih = fsig.first(); !ih.last(); ih.next() ) {
				isig[ih].I() = clipper::Util::mean(std::pow(fsig_ano[ih].f_pl(),2),std::pow(fsig_ano[ih].f_mi(),2));
				isig[ih].sigI() = clipper::Util::sig_mean(2.0f*fsig_ano[ih].f_pl()*fsig_ano[ih].sigf_pl(),2.0f*fsig_ano[ih].f_mi()*fsig_ano[ih].sigf_mi(), 0.0f  );
				//isig[ih].I() = ctruncate::Utils::mean(std::pow(fsig_ano[ih].f_pl(),2),std::pow(fsig_ano[ih].f_mi(),2),2.0f*fsig_ano[ih].f_pl()*fsig_ano[ih].sigf_pl(),2.0f*fsig_ano[ih].f_mi()*fsig_ano[ih].sigf_mi());
				//isig[ih].sigI() = ctruncate::Utils::sig_mean(std::pow(fsig_ano[ih].f_pl(),2),std::pow(fsig_ano[ih].f_mi(),2),2.0f*fsig_ano[ih].f_pl()*fsig_ano[ih].sigf_pl(),2.0f*fsig_ano[ih].f_mi()*fsig_ano[ih].sigf_mi(), 0.0f  );
			}
		} else {
			for ( HRI ih = isig_ano_import.first(); !ih.last(); ih.next() ) {
				//isig[ih].I() = clipper::Util::mean(isig_ano_import[ih].I_pl(),isig_ano_import[ih].I_mi());
				//isig[ih].sigI() = clipper::Util::sig_mean(isig_ano_import[ih].sigI_pl(),isig_ano_import[ih].sigI_mi(), 0.0f );
				//isig[ih].I() = ctruncate::Utils::mean(isig_ano_import[ih].I_pl(),isig_ano_import[ih].I_mi(),isig_ano_import[ih].sigI_pl(),isig_ano_import[ih].sigI_mi());
				//isig[ih].sigI() = ctruncate::Utils::sig_mean(isig_ano_import[ih].I_pl(),isig_ano_import[ih].I_mi(),isig_ano_import[ih].sigI_pl(),isig_ano_import[ih].sigI_mi(),0.0f );
				isig[ih] = clipper::data32::I_sigI(isig_ano_import[ih].I(),isig_ano_import[ih].sigI());
			}				
		}
	}
	
    
    int Ncentric = 0;
    int Nreflections = 0;
    for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
        if ( !isig[ih].missing() ) Nreflections++;
        if ( ih.hkl_class().centric() && !isig[ih].missing()) Ncentric++;
    }
    
    prog.summary_beg();
    
    printf("Cell parameters: %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n", cell1.a(), cell1.b(), cell1.c(), 
           Util::rad2d( cell1.alpha() ), Util::rad2d( cell1.beta() ), Util::rad2d( cell1.gamma() ) );
    printf("\nNumber of reflections: %d\n", Nreflections);
	printf("\nNumber of centric reflections %d\n", Ncentric);
    clipper::Grid_sampling grid;
    
    // can't seem to get max resolution from clipper, so use CCP4 methods
    CMtz::MTZ *mtz1=NULL;
    int read_refs=1;  // not sure what read_refs actually does - reads reflections presumably
    float minres,maxres;
    mtz1 = CMtz::MtzGet(args[mtzinarg].c_str(), read_refs);
    
    // read title
    char title[72];
    CMtz::ccp4_lrtitl(mtz1, title);
    
    CMtz::MtzResLimits(mtz1,&minres,&maxres);
    float invopt = maxres;
    float resopt = 1.0/sqrt(invopt);
    printf("\nMinimum resolution = %7.3f A\nMaximum resolution = %7.3f A\n",1.0/sqrt(minres),1.0/sqrt(maxres));
    prog.summary_end();
    CSym::CCP4SPG *spg1 = CSym::ccp4spg_load_by_ccp4_num(CMtz::MtzSpacegroupNumber(mtz1));
    prog.summary_beg();
    CMtz::MtzFree( mtz1 );
	
    
    char spacegroup[20];
    strcpy(spacegroup,spg1->symbol_old);
    printf("\nSpacegroup: %s (number %4d)\n", spg1->symbol_old, spg1->spg_ccp4_num);
    
    char pointgroup[20];
    strcpy(pointgroup,spg1->point_group);
    printf("Pointgroup: %s\n\n",pointgroup);
    prog.summary_end();
    
    if ( ipseq != "NONE" ) {
        prog.summary_beg(); printf("CELL CONTENTS:\n\n");
        
        clipper::SEQfile seqf;
        seqf.read_file( ipseq );
        
        ctruncate::Matthews cmath(true,false);
        int nmol = cmath(cell1, spgr, seqf, resopt);
        std::cout << "Expected number of molecules in ASU : " << nmol << std::endl;
        prog.summary_end();
        cmath.summary();
    } else if (nresidues > 0) {		
        prog.summary_beg(); printf("CELL CONTENTS:\n\n");
        ctruncate::Matthews cmath(true,false);
        int nmol = cmath(cell1, spgr, nresidues, resopt);
        std::cout << "Expected number of molecules in ASU : " << nmol << std::endl;
        prog.summary_end();
        cmath.summary();
    }
    
    //Completeness information
    //what is our working resolution (85% of I/sigI > 3.0)
    clipper::Range<double> reso_range;
    int NBINS = 60;
    double ACCEPTABLE = 0.85;
    ctruncate::Completeness<data32::I_sigI> compt(NBINS);
    compt(isig);
    compt.plot();
    {
        clipper::Range<double> range(minres,hklinf.invresolsq_range().max() );
        int i = 0;
        for ( ; i != NBINS-1 ; ++i) {
            if ( compt.completeness3(compt.bin2invresolsq(i)) > ACCEPTABLE && compt.completeness3(compt.bin2invresolsq(i+1)) > ACCEPTABLE ) break;
        }
        if ( i != (NBINS-1) ) {
            int j = NBINS-1;
            for ( ; j != 1 ; --j) {
                if ( compt.completeness3(compt.bin2invresolsq(j)) > ACCEPTABLE && compt.completeness3(compt.bin2invresolsq(j-1)) > ACCEPTABLE ) break;
            }
            if (j != 0 )
				if (i != 0) {
					float d = (compt.bin2invresolsq(i)+compt.bin2invresolsq(i-1))/2.0;
					reso_range.include(d);
				} else {
					reso_range.include(minres );
				}
            if (j != NBINS-1 ) {
                float d = (compt.bin2invresolsq(j)+compt.bin2invresolsq(j+1))/2.0;
                reso_range.include(d);	
            } else {
                reso_range.include(range.max() );
            }
        } else {
			reso_range = range;
        }
    }
    
    prog.summary_beg();
    double rr = 0.0;
    printf("\nCOMPLETENESS ANALYSIS (using intensities):\n");
    if ( reso_range.max() == -999999999 && reso_range.min() == 999999999 ) {
        printf("WARNING: The resolution range with I/sigI > 3  and completeness above %4.2f could not be\n",ACCEPTABLE);
        printf("         determined.  This data is of poor quality.\n\n");
    } else {
        printf("\nUsing I/sigI > 3 with completeness above %4.2f, the estimated useful\nResolution Range ",ACCEPTABLE);
        printf("of this data is %7.3fA to %7.3fA\n",1.0/std::sqrt(reso_range.min() ), 1.0/std::sqrt(reso_range.max() ) );
        rr = 1.0/std::sqrt(reso_range.min() ) - 1.0/std::sqrt(reso_range.max() );
    }
    prog.summary_end();
    printf("\n");
	
    // how big is reso_range
    {
        float rmax = hklinf.resolution().limit();
        float rmin = 1.0/std::sqrt(minres);
        float amax = 1.0/std::sqrt(reso_range.max() );
        if ( rr < 4.0 || amax > std::max(6.0,rmax+2.0)  ) {
            reso_range = clipper::Range<double>();
            printf("WARNING: The resolution range with I/sigI > 3 with completeness above 0.85 is small\n");
            do {
                printf("         Recalculating using I/sigI > 2 with completeness aboue %4.2f for\n",ACCEPTABLE);
                printf("         use in the following statistics, which must be treated with extreme caution\n");
                int i = 0;
                clipper::Range<double> range(minres,hklinf.invresolsq_range().max() );
                for ( ; i != NBINS-1 ; ++i) {
                    if ( compt.completeness2(compt.bin2invresolsq(i)) > ACCEPTABLE && compt.completeness2(compt.bin2invresolsq(i+1)) > ACCEPTABLE ) break;
                }
                if ( i != (NBINS-1) ) {
                    int j = NBINS-1;
                    for ( ; j != 1 ; --j) {
                        if ( compt.completeness2(compt.bin2invresolsq(j)) > ACCEPTABLE && compt.completeness2(compt.bin2invresolsq(j-1)) > ACCEPTABLE ) break;
                    }
                    if (j != 0 )
                        if (i != 0) {
                            float d = (compt.bin2invresolsq(i)+compt.bin2invresolsq(i-1))/2.0;
                            reso_range.include(d);
                        } else {
                            reso_range.include(range.min() );
                        }
                    if (j != NBINS-1 ) {
                        float d = (compt.bin2invresolsq(j)+compt.bin2invresolsq(j+1))/2.0;
                        reso_range.include(d);
                    } else {
                        reso_range.include(range.max() );
                    }
                }
                amax = 1.0/std::sqrt(reso_range.max() );
                rr = ( reso_range.max() == -999999999 && reso_range.min() == 999999999 ) ? 0.0 : 1.0/std::sqrt(reso_range.min() ) - amax;
                ACCEPTABLE -= 0.1;
                printf("         Resolution Range of this data for these values is %7.3fA to %7.3fA\n",1.0/std::sqrt(reso_range.min() ), amax );
                if (ACCEPTABLE <= 0.40) break;
            } while ( rr < 4.0 || amax > std::max(6.8,rmax+2.0) );
            if (ACCEPTABLE <= 0.40) reso_range = clipper::Range<double>();
            
            // try on Istandard
            if (reso_range.max() == -999999999 && reso_range.min() == 999999999 ) {
                printf("         Attempt resolution range estimate using Istandard < 1.0 (Ideally would use 0.2) \n");
                int i = 0;
                clipper::Range<double> range(minres,hklinf.invresolsq_range().max() );
                for ( ; i != NBINS-1 ; ++i) {
                    if ( compt.standard(compt.bin2invresolsq(i)) < 1.0 && compt.standard(compt.bin2invresolsq(i+1)) < 1.0 ) break;
                }
                if ( i != (NBINS-1) ) {
                    int j = NBINS-1;
                    for ( ; j != 1 ; --j) {
                        if ( compt.standard(compt.bin2invresolsq(j)) < 1.0 && compt.standard(compt.bin2invresolsq(j-1)) < 1.0 ) break;
                    }
                    if (j != 0 )
                        if (i != 0) {
                            float d = (compt.bin2invresolsq(i)+compt.bin2invresolsq(i-1))/2.0;
                            reso_range.include(d);
                        } else {
                            reso_range.include(minres );
                        }
                    if (j != NBINS-1 ) {
                        float d = (compt.bin2invresolsq(j)+compt.bin2invresolsq(j+1))/2.0;
                        reso_range.include(d);
                    } else {
                        reso_range.include(range.max() );
                    }
                } else {
                    reso_range = range;
                }
                if ( reso_range.max() != -999999999 && reso_range.min() != 999999999 ) {
                    printf("         Resolution Range of this data is %7.3fA to %7.3fA\n",1.0/std::sqrt(reso_range.min() ), 1.0/std::sqrt(reso_range.max() ) );
                }
            }
            
            if (reso_range.max() == -999999999 && reso_range.min() == 999999999 ) {
                reso_range.include(minres);
                reso_range.include(1.0/std::pow(rmax+2.0,2.0));
                printf("WARNING: Arbitary resolution range of %7.3fA to %7.3fA\n",1.0/std::sqrt(reso_range.min() ), 1.0/std::sqrt(reso_range.max() ) );
            }
        }
    }
    printf("\nThe high resolution cut-off will be used in gathering the statistics for the dataset, however the full dataset will be output\n\n");
	
    // limit resolution of Patterson calculation for tNCS (default 4 A), or set to
    // limit from completeness analysis
    invopt = reso_range.max();
    resopt = 1.0/std::sqrt(invopt);
    reso_Patt = clipper::Resolution( std::max( double(resopt), reso_Patt.limit() ) );
    //user override of defaults
    if (!reso_u1.is_null() ) {
        reso_Patt = clipper::Resolution( clipper::Util::max( reso.limit(), reso_u1.limit() ) );
    }
	
    // check for pseudo translation (taken from cpatterson)
    ctruncate::tNCS<float> tncs;
    const std::vector<clipper::Symop>& cf = tncs(isig, reso_Patt );
	
    prog.summary_beg();
    tncs.summary();
    prog.summary_end();
    printf("\n");
	
	//setup aniso copy of isig
	HKL_data<data32::I_sigI> ianiso(hklinf);
	for ( HRI ih = ianiso.first(); !ih.last(); ih.next() ) {  
		if (reso_range.contains(ih.invresolsq() ) ) {
			float I = isig[ih.hkl()].I();
			float sigI = isig[ih.hkl()].sigI();
			ianiso[ih] = clipper::data32::I_sigI( I, sigI );
		}
	}
	
    // anisotropy estimation
    clipper::U_aniso_orth uao;
    clipper::U_aniso_orth uaoc(0,0,0,0,0,0);
	bool anisobysymm(false);
	bool anisodemo(false);
    {
		prog.summary_beg();
		printf("\nANISOTROPY ANALYSIS (using intensities):\n");
		
		
		
		
		clipper::U_aniso_orth uao_sum(1,3,5,7,11,13);
		for (int i = 1; i != spgr.num_symops() ; ++i ) {
			clipper::U_aniso_orth uao(1,3,5,7,11,13);
			clipper::U_aniso_orth tmp = uao.transform(spgr.symop(i).rtop_orth(cell1) );
			uao_sum = uao_sum + tmp;
		}
		
		if (std::fabs(uao_sum.mat00() - uao_sum.mat11() ) > 0.5 || std::fabs(uao_sum.mat00() - uao_sum.mat22()) > 0.5 ) {
			anisobysymm = true;
			try { 
				AnisoCorr<Iscale_logLikeAniso<float>, clipper::datatypes::I_sigI<float>, float > llscl(ianiso, false, false, reso_range);
				uao = -(llscl.u_aniso_orth(Scaling::I) );
				uaoc = -(llscl.u_aniso_orth(Scaling::F) );
			} catch (clipper::Message_fatal) {
				CCP4::ccperror(1, "Anisotropy anlysis failed.");
			}
			
			// Eigenvalue calculation
			AnisoDirection<float> direction(uao);
			
			std::vector<float> v = direction.eigenValues();
			float max = std::exp(direction.max() );
			
			printf("\nEigenvalues: %8.4f %8.4f %8.4f\n", v[0],v[1],v[2]);
			printf("Eigenvalue ratios: %8.4f %8.4f %8.4f\n", std::exp(v[0])/max, std::exp(v[1])/max, std::exp(v[2])/max);
			float ratio = std::exp(std::min(v[0],std::min(v[1],v[2]) ) );
			if ( ratio < 0.5 ) printf("\nWARNING! WARNING! WARNING! Your data is severely anisotropic\n");
			prog.summary_end();
			printf("\n");
			
			printf("\nAnisotropic U (orthogonal coords):\n\n");
			printf("| %8.4f %8.4f %8.4f |\n", uao(0,0) ,  uao(0,1) ,  uao(0,2)  );
			printf("| %8.4f %8.4f %8.4f |\n", uao(1,0) ,  uao(1,1) ,  uao(1,2)  );
			printf("| %8.4f %8.4f %8.4f |\n", uao(2,0) ,  uao(2,1) ,  uao(2,2)  );
			
			clipper::U_aniso_frac uaf = uao.u_aniso_frac( cell1 );
			
			printf("\nAnisotropic U scaling (fractional coords):\n\n"); 
			
			printf("| %11.3e %11.3e %11.3e |\n", uaf(0,0) ,  uaf(0,1) ,  uaf(0,2)  );
			printf("| %11.3e %11.3e %11.3e |\n", uaf(1,0) ,  uaf(1,1) ,  uaf(1,2)  );
			printf("| %11.3e %11.3e %11.3e |\n", uaf(2,0) ,  uaf(2,1) ,  uaf(2,2)  );
			
			printf("\nAnisotropic B scaling (fractional coords):\n\n"); 
			
			printf("| %11.3e %11.3e %11.3e |\n",clipper::Util::u2b( uaf(0,0) ), clipper::Util::u2b( uaf(0,1) ), clipper::Util::u2b( uaf(0,2) ) );
			printf("| %11.3e %11.3e %11.3e |\n",clipper::Util::u2b( uaf(1,0) ), clipper::Util::u2b( uaf(1,1) ), clipper::Util::u2b( uaf(1,2) ) );
			printf("| %11.3e %11.3e %11.3e |\n",clipper::Util::u2b( uaf(2,0) ), clipper::Util::u2b( uaf(2,1) ), clipper::Util::u2b( uaf(2,2) ) );
			// demonstratable anisotropy
			{
				clipper::ftype v1(v[0]/max), v2(v[1]/max), v3(v[2]/max);
				if (std::abs(v1-v2) > 0.01 || std::abs(v1-v3) > 0.01) anisodemo = true;
			}
		} else { 
			printf("\nNo anisotropy by symmetry. \n");
			prog.summary_end();
		}
		
		
		// falloff calculation (Yorgo Modis)
		YorgoModis<data32::I_sigI> ym(resopt,60,uao);
		ym(isig);
		ym.plot();		
	}
	
	//want to use anisotropy correction and resolution truncation for twinning tests
	
	float lval(0.0);
    {
		printf("\nTWINNING ANALYSIS:\n\n");
		
		clipper::Range<clipper::ftype> range_Twin(reso_range.min(),
												  (!reso_u2.is_null() ) ? 
												  std::min(reso_range.max(),1.0/std::pow(reso_u2.limit(),2) ) :
												  reso_range.max() );
		
		printf("\nData has been truncated at %6.2f - %6.2f A resolution\n",1.0/std::sqrt(range_Twin.min()), 1.0/std::sqrt(range_Twin.max()));
		printf("Anisotropy correction has been applied before calculating twinning tests\n\n");
		
		
		
		TwinSymops ts1(cell1,spgr);
		
		L_test ltest;
		lval=ltest(ianiso,const_cast<std::vector<clipper::Symop> & >(cf), range_Twin);
		ltest.summary();
		ltest.loggraph();
		
		Moments<data32::I_sigI> m(isig,range_Twin);
		m.loggraph();
		
		printf("\nMean acentric moments I from input data:\n\n");
		printf("  <I^2>/<I>^2 = %6.3f (Expected = %6.3f, Perfect Twin = %6.3f)\n", m.acentric_second(), m.theo_untwinned_acentric_second(), m.theo_perfect_acentric_second() );
		printf("  <I^3>/<I>^3 = %6.3f (Expected value = %6.3f, Perfect Twin = %6.3f)\n", m.acentric_third(), m.theo_untwinned_acentric_third(), m.theo_perfect_acentric_third() );
		printf("  <I^4>/<I>^4 = %6.3f (Expected value = %6.3f, Perfect Twin = %6.3f)\n", m.acentric_fourth(), m.theo_untwinned_acentric_fourth(), m.theo_perfect_acentric_fourth());
		float m_fraction = m.fraction();
		
		if (anisobysymm && anisodemo) {
			Moments<data32::I_sigI> mc(ianiso,range_Twin);
			printf("\n\nMean acentric moments I from anisotropically corrected data:\n\n");
			printf("  <I^2>/<I>^2 = %6.3f (Expected = %6.3f, Perfect Twin = %6.3f)\n", mc.acentric_second(), mc.theo_untwinned_acentric_second(), mc.theo_perfect_acentric_second() );
			printf("  <I^3>/<I>^3 = %6.3f (Expected value = %6.3f, Perfect Twin = %6.3f)\n", mc.acentric_third(), mc.theo_untwinned_acentric_third(), mc.theo_perfect_acentric_third() );
			printf("  <I^4>/<I>^4 = %6.3f (Expected value = %6.3f, Perfect Twin = %6.3f)\n", mc.acentric_fourth(), mc.theo_untwinned_acentric_fourth(), mc.theo_perfect_acentric_fourth());	
			m_fraction = mc.fraction();
		}
		std::cout << std::endl << std::endl;
		
		std::vector<clipper::ftype> hval(ts1.size() );
		std::vector<clipper::ftype> bval(ts1.size() );
		std::vector<clipper::ftype> mval(ts1.size() );
		std::vector<clipper::ftype> mrval(ts1.size() );
		
		std::vector<H_test> htests(ts1.size() );
		for (int i = 0; i != ts1.size() ; ++i ) {
			hval[i]=htests[i](ianiso,ts1[i],range_Twin);
			//htests[i].summary();
			htests[i].loggraph();
		}
		
		std::vector<Britton_test> btests(ts1.size() );
		for (int i = 0; i != ts1.size() ; ++i ) {
			bval[i]=btests[i](ianiso,ts1[i],range_Twin);
			//btests[i].summary();
			btests[i].loggraph();
		}
		
		std::vector<MLBritton_test> mdtests(ts1.size() );
		for (int i = 0; i != ts1.size() ; ++i ) {
			clipper::ftype product(0.0);
			int jp;
			if ( tncs.hasNCS() ) {
				for (int j=0; j != tncs.numOps() ; ++j) {
					clipper::Vec3<clipper::ftype> vect = tncs[i].rtop_orth(cell1).trn();
					clipper::Mat33<int> tmp = ts1[i].rot(); 
					clipper::Rotation rot(clipper::Mat33<clipper::ftype>(
																		 clipper::ftype(tmp(0,0) )/12.0,
																		 clipper::ftype(tmp(0,1) )/12.0,
																		 clipper::ftype(tmp(0,2) )/12.0,
																		 clipper::ftype(tmp(1,0) )/12.0,
																		 clipper::ftype(tmp(1,1) )/12.0,
																		 clipper::ftype(tmp(1,2) )/12.0,
																		 clipper::ftype(tmp(2,0) )/12.0,
																		 clipper::ftype(tmp(2,1) )/12.0,
																		 clipper::ftype(tmp(2,2) )/12.0));
					clipper::Vec3<clipper::ftype> vecr(rot.x(),rot.y(),rot.z() );
					clipper::ftype p = vecr.unit()*vect.unit();
					if (p > product) {
						product = p;
						jp = j;
					}
				}
			}
			if (product > 0.95 ) mval[i]=mdtests[i](ianiso,ts1[i],cf[jp],range_Twin);
			else mval[i]=mdtests[i](ianiso,ts1[i],range_Twin);
			mrval[i]=mdtests[i].deltaR();
			//mdtests[i].summary();
			//mdtests[i].loggraph();
		}
		
		std::cout << "Twin fraction estimates excluding operators" << std::endl;
		std::cout << "  Twin fraction estimate from L-test:  " << std::setw(4) << std::setprecision(2) << ltest.fraction() << std::endl;
		std::cout << "  Twin fraction estimate from moments: " << std::setw(4) << std::setprecision(2) << m_fraction << std::endl << std::endl;
		
		std::cout << "Twin fraction estimates by operator" << std::endl << std::endl;
		if ( ts1.size() > 0 ) {
			std::cout << "---------------------------------------------------------------------------------------" << std::endl;
			std::cout << "| " << std::setw(40) << "operator" <<         " | L-test | H-test | Murray | ML Britton    |" << std::endl;
			std::cout << "---------------------------------------------------------------------------------------" << std::endl;
			for (int i = 0; i != ts1.size() ; ++i ) {
				std::cout << "| " << std::setw(40) << htests[i].description() << " |  " << ((0.1 <= lval && lval < 0.440) ? " Yes " : " No  " ) 
				<< " |  " << std::setw(4) << std::setprecision(2) << hval[i] << "  |  " << bval[i] << "  |  " << mval[i] << " (";
				if (mrval[i] == 100.0 ) std::cout << " N/A ";
				else std::cout << std::setw(5) << mrval[i];
				std::cout << ") |" << std::endl;
			}
			std::cout << "---------------------------------------------------------------------------------------" << std::endl;
		} else {
			std::cout << "  No operators found" << std::endl << std::endl;
		}
		
		prog.summary_beg();
		if (ts1.size() == 0 ) twin_summary(0.0,lval);
		else twin_summary((*(std::max_element(hval.begin(),hval.end()))),lval);
		prog.summary_end(); 
		printf("\n");
		
		
		//printf("Starting parity group analysis:\n");
		
		//Parity group analysis
		
		ctruncate::parity(ianiso, invopt, nbins);	
    }
    // Ice rings
    ctruncate::Rings icerings;
    icerings.DefaultIceRings();
    icerings.ClearSums();
	
	{
		ctruncate::IceRings_analyse ice;
		bool icer=ice(isig,icerings);
		
        printf("\n");
        prog.summary_beg();
        printf("\nICE RINGS:\n\n");
        if (icer) printf("Possible Ice Rings\n\n");
        else printf("No Ice Rings detected\n\n");
        prog.summary_end();
        printf("\n");
		
		std::cout << ice.format() << std::endl;
		
	}
	
	// if something went wrong with Wilson scaling, B could be negative, giving exponentially large scaled SF's
	// so only scale if B positive
	float scalef = 1.0;
	{
		//Wilson plot
		//std::vector<float> wilson(2,0.0f);
		//nprm = 60*std::max(int(sqrt(float(Nreflections))),nbins );
		clipper::MMoleculeSequence seq;
		
		WilsonB::MODE wilson_flag;
		if (is_nucl) {
			wilson_flag = WilsonB::RNA;
		} else { 
			wilson_flag = WilsonB::BEST;
		}
		WilsonB wilson( wilson_flag);
		if ( ipseq != "NONE" ) {
			clipper::SEQfile seqf;
			seqf.read_file( ipseq );
			seqf.import_molecule_sequence( seq );
			MPolymerSequence poly = seq[0];
			wilson(isig,poly,&reso_range, &icerings);
		} else if (nresidues > 0) {
			wilson(isig,nresidues,&reso_range, &icerings);
		} else {
			wilson(isig,&reso_range,&icerings);
		}		
		
		clipper::String comment("Smooth");
		printf("\n");
		prog.summary_beg();
		wilson.summary();
		prog.summary_end();
		printf("\n");
		wilson.plot();
		
		if ( wilson.intercept() > 0 ) scalef = sqrt(wilson.intercept() );
	}
	
		
	HKL_data<data32::I_sigI> xsig(hklinf);
	//normal calculation
	{		
		for ( HRI ih = xsig.first(); !ih.last(); ih.next() ) {
			double reso = ih.invresolsq();
			xsig[ih] = clipper::data32::I_sigI( (isig[ih.hkl()].I()), isig[ih.hkl()].sigI() );
			if ( icerings.InRing(reso) != -1 ) 
				if ( icerings.Reject( icerings.InRing(reso) ) ) xsig[ih].I() = xsig[ih].sigI() = clipper::Util::nan(); // loose ice rings
		}

		
		
		int nprm2=12;
		int nreflns=500;
		
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
			tr1[ih] = clipper::data32::I_sigI( xsig[ih.hkl()].I()/ctruncate::BEST(reso), 0.0);
		}		
		
		int nrej_pre(0);
		std::vector<double> params(nprm2,1.0);
		{
			//precondition the ML calc using least squares fit.  This should give an excellent start
			clipper::BasisFn_binner basis_pre( tr1, nprm2, 1.0 );
			TargetFn_meanInth<clipper::data32::I_sigI> target_pre(tr1,1.0);
			clipper::ResolutionFn pre( hklinf, basis_pre, target_pre, params);
			params = pre.params();
			
			//outlier rejection from Read (1999)			
			double rlimit(0.0);
			for ( HRI ih = tr1.first(); !ih.last(); ih.next() ) {
				if ( !xsig[ih.hkl()].missing() ) {
				double I = tr1[ih].I()/(ih.hkl_class().epsilon()*/*tr1[ih].I()*/pre.f(ih));
				rlimit = (ih.hkl_class().centric() ) ? 6.40*6.40 : 4.55*4.55 ;
					if (I > rlimit )  {
						++nrej_pre;
						}
				}
			}				
			std::cout << "Number of outliers not used in norm calculation (Read (1999) ): " << nrej_pre << std::endl;
			//reset tr1
			double reso(0.0);
			for ( HRI ih = tr1.first(); !ih.last(); ih.next() ) {
				reso = ih.invresolsq();
				tr1[ih] = clipper::data32::I_sigI( ctruncate::BEST(reso), 0.0);
			}					
		}
		std::vector<bool> mask(nprm2,false);
		clipper::BasisFn_spline basis_fo( xsig, nprm2, 1.0 );
		//TargetFn_meanInth<clipper::data32::I_sigI> target_fo( xsig, 1 );
		TargetFn_scaleLogLikeI1I2<clipper::data32::I_sigI,clipper::data32::I_sigI> target_fo(tr1,xsig);
		ctruncate::ResolutionFn_nonlinear Sigma( hklinf, basis_fo, target_fo, params, mask, 1.0, false);
		params = Sigma.params();
		
		for ( HRI ih = xsig.first(); !ih.last(); ih.next() ) {
			double reso = ih.invresolsq();
			xsig[ih] = clipper::data32::I_sigI(Sigma.f(ih) * tr1[ih].I(),1.0f);
		}
		
		// scale the norm for the anisotropy
		if (aniso) {            
			if (anisobysymm && anisodemo) {
				clipper::datatypes::Compute_scale_u_aniso<clipper::data32::I_sigI > compute_s(1.0,-uaoc);
				xsig.compute(xsig, compute_s);
				//printf("\nAnisotropic U (orthogonal coords):\n\n");
				//printf("| %8.4f %8.4f %8.4f |\n", uaoc(0,0) ,  uaoc(0,1) ,  uaoc(0,2)  );
				//printf("| %8.4f %8.4f %8.4f |\n", uaoc(1,0) ,  uaoc(1,1) ,  uaoc(1,2)  );
				//printf("| %8.4f %8.4f %8.4f |\n", uaoc(2,0) ,  uaoc(2,1) ,  uaoc(2,2)  );
				
			}
		}		
		
		{
			int nrej(0);
			//outlier rejection from Read (1999) using new norm		
			double rlimit(0.0);
			for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
				if ( !isig[ih].missing() ) {
					double I = isig[ih].I()/(ih.hkl_class().epsilon()*xsig[ih.hkl()].I() );
					rlimit = (ih.hkl_class().centric() ) ? 6.40*6.40 : 4.55*4.55 ;
					if (I > rlimit )  {
						++nrej;
					}
				}
			}	
			std::cout << "Number of outliers (Read (1999) ): " << nrej << std::endl << std::endl;
			if (nrej > nrej_pre && ( prior != FLAT || prior != SIVIA ) ) std::cout << std::endl << "WARNING: prior may be unstable" << std::endl << std::endl;
		}
	}
		
	int nrej = 0; 
	
	if (!amplitudes) {
		
		//user override of truncate procedure and output
		if (!reso_u3.is_null() ) {
			reso_trunc = 
			clipper::Resolution( clipper::Util::max( reso_trunc.limit(), reso_u3.limit() ) );
		}
		if ( prior == AUTO && (tncs.hasNCS() || 0.440 > lval) ) {
			printf("\nWARNING: FLAT prior in use due to either tNCS or twinning.\nTo override force --prior WILSON\n\n");
			prior = FLAT;
		}
		
		if ( refl_mean ) {
			if (prior == FLAT ) truncate( isig, jsig, fsig, scalef, spg1, reso_trunc, nrej, debug );
			else if (prior == SIVIA) truncate_sivia(isig, jsig, fsig, scalef, spg1, reso_trunc, nrej, debug );
			else truncate( isig, jsig, fsig, xsig, scalef, spg1, reso_trunc, nrej, debug );
		}
		if (anomalous) {
			if (prior == FLAT ) truncate( isig_ano_import, jsig_ano, fsig_ano, scalef, spg1, reso_trunc, nrej, debug );
			else if (prior == SIVIA) truncate_sivia( isig_ano_import, jsig_ano, fsig_ano, scalef, spg1, reso_trunc, nrej, debug );
            else truncate( isig_ano_import, jsig_ano, fsig_ano, xsig, scalef, spg1, reso_trunc, nrej, debug );
			int iwarn = 0;
			for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
				freidal_sym[ih].isym() = ( !Util::is_nan(fsig_ano[ih].f_pl() )  &&  !Util::is_nan(fsig_ano[ih].f_mi() ) ) ? 0 :
				( !Util::is_nan(fsig_ano[ih].f_pl() ) ) ? 1 :
				( !Util::is_nan(fsig_ano[ih].f_mi() ) ) ? 2 : 0;
				Dano[ih].d() = fsig_ano[ih].d();
				Dano[ih].sigd() = fsig_ano[ih].sigd();
				if ( ih.hkl_class().centric() ) {
					Dano[ih].d() = 0.0;
					Dano[ih].sigd() = 0.0;
				}
			}
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
	}
	printf("\n");
	prog.summary_beg();
	printf("\nINTENSITY TO AMPLITUDE CONVERSION:\n\n");
	if ( prior == FLAT ) printf("Calculation using flat prior\n");
	else printf("Calculation using Wilson prior\n");
	
	printf("%d intensities have been rejected as unphysical\n", nrej);
	prog.summary_end();
	printf("\n");
	
	
	
	// following code is for when truncate calc switched off - do not delete it
	// usually already have F's in this case; really just need to skip truncate calc
	
	/*
	 for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
	 if ( !isig[ih].missing() ) {
	 float I = isig[ih].I();
	 float sigma = isig[ih].sigI();
	 HKL hkl = ih.hkl();
	 float weight = (float) CSym::ccp4spg_get_multiplicity( spg1, hkl.h(), hkl.k(), hkl.l() );
	 float sqwt = sqrt(weight);
	 if (I < 0.0) {
	 fsig[ih].f() = 0.0;
	 fsig[ih].sigf() = 0.0;
	 }
	 else {
	 fsig[ih].f() = sqrt(I)*scalef*sqwt;
	 fsig[ih].sigf() = 0.5*(sigma/sqrt(I))*scalef*sqwt; //check this
	 }
	 }
	 }*/
	
	
	// moments of E using clipper binning
	// moments_Z(ianiso,resopt,nbins,prog);
	
	
	// construct cumulative distribution function for intensity (using Z rather than E)
	int ntw = cumulative_plot(isig, xsig);
	
	if (ntw > 2) {
		prog.summary_beg();
		printf("\nWARNING: ****  Cumulative Distribution shows Possible Twinning ****\n\n");
		prog.summary_end();
	}
	
	
	// falloff calculation (Yorgo Modis)
	//yorgo_modis_plot(fsig,maxres,60,prog,uao);
	
	// anomalous signal
	if (anomalous ) {
		if (amplitudes) AnomStats<float> anomstats(fsig_ano);
        else AnomStats<float> anomstats(isig_ano_import);
	}
	
	
	{
		ctruncate::PattPeak patt_peak(std::sqrt(maxres));
		
		float opt_res = patt_peak(xsig);
		
		float width_patt = 2.0f*patt_peak.sigma();
		
		float b_patt = 4.0f*clipper::Util::twopi2()*std::pow(width_patt/2.0f,2.0f);
		
		prog.summary_beg();
		std::cout << "Estimated Optical Resolution: " << opt_res << std::endl;
		prog.summary_end();
		
		
		
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
        //mtzout.open_append( args[mtzinarg], outfile );
        mtzout.open_write( outfile );
        mtzout.export_crystal ( cxtl, outcol );
        mtzout.export_dataset ( cset, outcol );
        mtzout.export_hkl_info( hkl_list );
        //mtzout.export_hkl_data( jsig, outcol );
        clipper::String labels;
		if ( refl_mean ) {
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
            if ( !refl_mean ) {
                if (appendcol == "") labels = outcol + "[FMEAN,SIGFMEAN]";
                else labels = outcol + "[FMEAN_" + appendcol + ",SIGFMEAN_" + appendcol + "]";
                mtzout.export_hkl_data( fsig, labels );
            }
            if (appendcol == "") labels = outcol + "[DANO,SIGDANO]";
            else labels = outcol + "[DANO_" + appendcol + ",SIGDANO_" + appendcol + "]";
            mtzout.export_hkl_data( Dano, labels );
            if (appendcol == "") labels = outcol + "[F(+),SIGF(+),F(-),SIGF(-)]";
            else labels = outcol + "[F_" + appendcol + "(+),SIGF_" + appendcol + "(+),F_" + appendcol + "(-),SIGF_" + appendcol + "(-)]";
            mtzout.export_hkl_data( fsig_ano, labels );
            if (appendcol == "") labels = outcol + "[ISYM]";
            else labels = outcol + "[ISYM_" + appendcol + "]";
            mtzout.export_hkl_data( freidal_sym, labels );
        }
		
		//output original input
		if ( refl_mean ) {
			if (appendcol != "") {
				String::size_type loc = meancol.find(",",0);
				meancol.insert(loc,"_"+appendcol);
				loc = meancol.find("]",0);
				meancol.insert(loc,"_"+appendcol);
			}
			mtzout.export_hkl_data( isig, outcol + meancol.tail() );
        }
		
        if (anomalous) {
            if (appendcol != "") {
                String::size_type loc = anocols.find("+",0);
                anocols.insert(loc-1,"_"+appendcol);
                loc = anocols.find(",",0);
                loc = anocols.find("+",loc+1);
                anocols.insert(loc-1,"_"+appendcol);
                loc = anocols.find("-",0);
                anocols.insert(loc-1,"_"+appendcol);
                loc = anocols.find(",",loc);
                loc = anocols.find("-",loc+1);
                anocols.insert(loc-1,"_"+appendcol);
            }
            mtzout.export_hkl_data( isig_ano_import, outcol + anocols.tail() );
        }
		
        //copy old history and say something about ctruncate run
        if (history.size() != 0 ) {
			for (int i = 0 ; i != history.size() ; ++i ) histin.push_back(history[i]);
        } else {
            char run_date[10];
            CCP4::ccp4_utils_date(run_date);
            char run_time[8];
            CCP4::ccp4_utils_time(run_time);
            clipper::String run_type = ( prior == FLAT ) ? " flat " : " french-wilson ";
			
            clipper::String history = prog_string + " " + prog_vers + run_type +  "run on " + run_date + " " + run_time;
            histin.push_back(history);
        }
        mtzout.set_history(histin);
		
        mtzout.set_spacegroup_confidence(spgr_confidence);
        
        //mtzout.close_append();
        mtzout.close_write();
        
        // Clipper will change H3 to R3, so change it back
        if (spacegroup[0] == 'H') {
			CMtz::MTZ *mtz2=NULL;
			read_refs=1;  // need to read in reflections, otherwise they won't be written out
			mtz2 = CMtz::MtzGet(outfile.c_str(), read_refs);
			// write title to output file
			strncpy( mtz2->title, title, 71 );
			strcpy(mtz2->mtzsymm.spcgrpname,spacegroup);
			CMtz::MtzPut( mtz2, outfile.c_str() );
			CMtz::MtzFree( mtz2 );
		}
    }
    prog.set_termination_message( "Normal termination" );
    
    return(0);
}





