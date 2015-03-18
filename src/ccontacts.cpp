// Clipper app to calculate interesting contact distances
// // See below for a complete description of the calculations
// // 2013 Jon Agirre & Kevin Cowtan @ University of York
// // mailto: jon.agirre@york.ac.uk
// // mailto: kevin.cowtan@york.ac.uk
// //

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <iostream>
#include <clipper/clipper.h>
#include <clipper/clipper-mmdb.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-minimol.h>


struct salt_bridge;

std::vector < salt_bridge > calculate_salt_bridges ( clipper::MiniMol&, clipper::MAtomNonBond& );
void print_salt_bridges ( std::vector < salt_bridge >, clipper::String, bool );

void insert_coot_prologue_scheme ( std::fstream& output );
void insert_coot_prologue_python ( std::fstream& output );
void insert_coot_files_loadup_scheme ( std::fstream& output, const clipper::String& pdb, bool mode );
void insert_coot_files_loadup_python ( std::fstream& output, const clipper::String& pdb, bool mode );
void insert_coot_go_to_res_scheme ( std::fstream& output, const clipper::Coord_orth& res_centre, const clipper::String& text );
void insert_coot_go_to_res_python ( std::fstream& output, const clipper::Coord_orth& res_centre, const clipper::String& text );
void insert_coot_epilogue_scheme ( std::fstream& output );
void insert_coot_epilogue_python ( std::fstream& output );


int main(int argc, char** argv)
{

    clipper::String program_version = "0.1";
    CCP4Program prog( "ccontacts", program_version.c_str(), "$Date: 2015/02/01" );

    prog.set_termination_message( "Failed" );
    std::cout << std::endl << "Copyright 2015 Jon Agirre, Kevin Cowtan and The University of York." << std::endl;
    std::cout << std::endl << "References:" << std::endl << std::endl << "- Salt bridge detection:\tKumar and Nussinov, 1999, JMB 293:1241-1255." << std::endl;
    std::cout << std::endl << "Please send any comments, bugs or feature requests to jon.agirre@york.ac.uk" << std::endl;
    clipper::String ippdb    = "NONE";
    clipper::String ipmode   = "long";
    clipper::String title    = "generic title";
    bool i2mode = false;

    // command input
    CCP4CommandInput args( argc, argv, true );
    int arg = 0;
    
    while ( ++arg < args.size() )
    {
        if ( args[arg] == "-title" )
        {
            if ( ++arg < args.size() )
               title = args[arg];
        }
        if ( args[arg] == "-pdbin" )
        {
            if ( ++arg < args.size() )
               ippdb = args[arg];
        }
        if ( args[arg] == "-typein" )
        {
            if ( ++arg < args.size() )
               ipmode = args[arg];
        }
        if ( args[arg] == "-mode" )
        {
            if ( ++arg < args.size() )
                if ( args[arg] == "ccp4i2" )
                   i2mode = true;
        } 
    }
   
    if (ippdb == "NONE") 
    {
        std::cout << "Usage: ccontacts" << std::endl << std::endl;
        std::cout << "\t-pdbin <.pdb>\t\tCOMPULSORY: input PDB file to examine" << std::endl;
        std::cout << "\t-typein\t\t\tlong,short,all" << std::endl;
        std::cout << "\t-mode\t\t\t{ccp4i2,normal}" << std::endl << std::endl << std::endl;
        std::cout << "\tThe program will also produce a visual tour of the detected contacts in the form" << std::endl;
        std::cout << "\tof Scheme and Python scripts for use with Coot" << std::endl;
        std::cout << "\n\tTo use them: 'coot --script ccontacts-results.scm' or 'coot --script ccontacts-results.py'" << std::endl;
        prog.set_termination_message( "Failed" );
        return 1;
    }

    std::cout << std::endl << "Reading " << ippdb.trim().c_str() << "... ";
    fflush(0);

    clipper::MMDBfile mfile;
    clipper::MiniMol mmol;

    const int mmdbflags = mmdb::MMDBF_IgnoreBlankLines | mmdb::MMDBF_IgnoreDuplSeqNum | mmdb::MMDBF_IgnoreNonCoorPDBErrors | mmdb::MMDBF_IgnoreRemarks | mmdb::MMDBF_EnforceUniqueChainID;
    mfile.SetFlag( mmdbflags );

    try
    {
        mfile.read_file( ippdb.trim() );
        mfile.import_minimol( mmol );
    }
    catch (...)
    {
        std::cout << std::endl << "There has been an unexpected error reading coordinates" << std::endl ;
    }
    
    std::cout << "done." << std::endl;

    clipper::MAtomNonBond manb = clipper::MAtomNonBond( mmol, 1.0 );

    std::vector < salt_bridge > results = calculate_salt_bridges ( mmol, manb );
    print_salt_bridges ( results, ippdb, i2mode );

    prog.set_termination_message( "Normal termination" );
    return 0;

}

struct salt_bridge 
{
    clipper::MMonomer positive, negative;
    clipper::MPolymer positive_chain, negative_chain;
    clipper::Coord_orth centroid_positive, centroid_negative;
    clipper::ftype distance;
    int symop_positive, symop_negative;
};

std::vector < salt_bridge > calculate_salt_bridges ( clipper::MiniMol& mmol, clipper::MAtomNonBond& manb )
{

    std::vector < salt_bridge > results;
    int natom_1, natom_2 = 0;
    clipper::MAtom atom_1, atom_2;
    clipper::Coord_orth centroid_positive, centroid_negative;

    for ( int poly = 0; poly < mmol.size(); poly++)
        for ( int mon = 0; mon < mmol[poly].size(); mon++)
            if ( mmol[poly][mon].type().trim() == "LYS" )
            {
                natom_1 = mmol[poly][mon].lookup("NZ",clipper::MM::ANY);
                    
                if ( natom_1 != -1 )
                {
                    atom_1 = mmol[poly][mon][natom_1];
                    centroid_positive = atom_1.coord_orth();
                    const std::vector<clipper::MAtomIndexSymmetry> neighbourhood = manb.atoms_near ( atom_1.coord_orth(), 2.8 );
                    
                    for ( int i = 0; i < neighbourhood.size() ; i++ )
                    {
                        clipper::MMonomer mon_tmp = mmol[neighbourhood[i].polymer()][neighbourhood[i].monomer()];

                        if ( mon_tmp.type().trim() == "ASP" )
                        {
                            int index_1 = mon_tmp.lookup("OD1",clipper::MM::ANY);
                            int index_2 = mon_tmp.lookup("OD2",clipper::MM::ANY);
                            
                            if ( ( index_1 != -1 ) && ( index_2 != -1 ) ) // everything alright, calculate centroid and add salt bridge if distance checks
                            {
                                const clipper::Coord_orth sum = mon_tmp[index_1].coord_orth() + mon_tmp[index_2].coord_orth();
                                centroid_negative = 0.5 * sum ;
                                
                                if ( clipper::Coord_orth::length (centroid_positive, centroid_negative ) < 4.0 )
                                {
                                    salt_bridge sb_tmp;
                                    sb_tmp.positive = mmol[poly][mon];
                                    sb_tmp.positive_chain = mmol[poly];
                                    sb_tmp.negative_chain = mmol[neighbourhood[i].polymer()];
                                    sb_tmp.negative = mon_tmp;
                                    sb_tmp.centroid_positive = centroid_positive;
                                    sb_tmp.centroid_negative = centroid_negative;
                                    sb_tmp.distance = clipper::Coord_orth::length (centroid_positive, centroid_negative );
                                    sb_tmp.symop_positive = 0;
                                    sb_tmp.symop_negative = neighbourhood[i].symmetry();
                                    results.push_back ( sb_tmp );
                                    break; // do not check the rest of the neighbours because we have what we were looking for
                                } 
                            }
                        }
                        else if ( mon_tmp.type().trim() == "GLU" )
                        {
                            int index_1 = mon_tmp.lookup("OE1",clipper::MM::ANY);
                            int index_2 = mon_tmp.lookup("OE2",clipper::MM::ANY);
                            
                            if ( ( index_1 != -1 ) && ( index_2 != -1 ) ) // everything alright, calculate centroid and add salt bridge if distance checks
                            {
                                const clipper::Coord_orth sum = mon_tmp[index_1].coord_orth() + mon_tmp[index_2].coord_orth();
                                centroid_negative = 0.5 * sum ;
                                
                                if ( clipper::Coord_orth::length (centroid_positive, centroid_negative ) < 4.0 )
                                {
                                    salt_bridge sb_tmp;
                                    sb_tmp.positive = mmol[poly][mon];
                                    sb_tmp.positive_chain = mmol[poly];
                                    sb_tmp.negative_chain = mmol[neighbourhood[i].polymer()];
                                    sb_tmp.negative = mon_tmp;
                                    sb_tmp.centroid_positive = centroid_positive;
                                    sb_tmp.centroid_negative = centroid_negative;
                                    sb_tmp.distance = clipper::Coord_orth::length (centroid_positive, centroid_negative );
                                    sb_tmp.symop_positive = 0;
                                    sb_tmp.symop_negative = neighbourhood[i].symmetry();
                                    results.push_back ( sb_tmp );
                                    break; // do not check the rest of the neighbours because we have what we were looking for
                                } 
                            }
                        }
                    }
                }
                else 
                    break;
            }
            else if ( mmol[poly][mon].type().trim() == "ARG" )
            {
                natom_1 = mmol[poly][mon].lookup("NH1",clipper::MM::ANY);
                natom_2 = mmol[poly][mon].lookup("NH2",clipper::MM::ANY);

                if (( natom_1 != -1 ) && ( natom_2 != -1 ))
                {
                    atom_1 = mmol[poly][mon][natom_1];
                    atom_2 = mmol[poly][mon][natom_2];
                    
                    const clipper::Coord_orth sum_positive = atom_1.coord_orth() + atom_2.coord_orth();
                    centroid_positive = 0.5 * sum_positive;
                    
                    const std::vector<clipper::MAtomIndexSymmetry> neighbourhood = manb.atoms_near ( atom_1.coord_orth(), 2.8 );
                    
                    for ( int i = 0; i < neighbourhood.size() ; i++ )
                    {
                        clipper::MMonomer mon_tmp = mmol[neighbourhood[i].polymer()][neighbourhood[i].monomer()];

                        if ( mon_tmp.type().trim() == "ASP" )
                        {
                            int index_1 = mon_tmp.lookup("OD1",clipper::MM::ANY);
                            int index_2 = mon_tmp.lookup("OD2",clipper::MM::ANY);
                            
                            if ( ( index_1 != -1 ) && ( index_2 != -1 ) ) // everything alright, calculate centroid and add salt bridge if distance checks
                            {
                                const clipper::Coord_orth sum = mon_tmp[index_1].coord_orth() + mon_tmp[index_2].coord_orth();
                                centroid_negative = 0.5 * sum ;
                                
                                if ( clipper::Coord_orth::length (centroid_positive, centroid_negative ) < 4.0 )
                                {
                                    salt_bridge sb_tmp;
                                    sb_tmp.positive = mmol[poly][mon];
                                    sb_tmp.positive_chain = mmol[poly];
                                    sb_tmp.negative_chain = mmol[neighbourhood[i].polymer()];
                                    sb_tmp.negative = mon_tmp;
                                    sb_tmp.centroid_positive = centroid_positive;
                                    sb_tmp.centroid_negative = centroid_negative;
                                    sb_tmp.distance = clipper::Coord_orth::length (centroid_positive, centroid_negative );
                                    sb_tmp.symop_positive = 0;
                                    sb_tmp.symop_negative = neighbourhood[i].symmetry();
                                    results.push_back ( sb_tmp );
                                    break; // do not check the rest of the neighbours because we have what we were looking for
                                } 
                            }
                        }
                        else if ( mon_tmp.type().trim() == "GLU" )
                        {
                            int index_1 = mon_tmp.lookup("OE1",clipper::MM::ANY);
                            int index_2 = mon_tmp.lookup("OE2",clipper::MM::ANY);
                            
                            if ( ( index_1 != -1 ) && ( index_2 != -1 ) ) // everything alright, calculate centroid and add salt bridge if distance checks
                            {
                                const clipper::Coord_orth sum = mon_tmp[index_1].coord_orth() + mon_tmp[index_2].coord_orth();
                                centroid_negative = 0.5 * sum ;
                                
                                if ( clipper::Coord_orth::length (centroid_positive, centroid_negative ) < 4.0 )
                                {
                                    salt_bridge sb_tmp;
                                    sb_tmp.positive = mmol[poly][mon];
                                    sb_tmp.positive_chain = mmol[poly];
                                    sb_tmp.negative_chain = mmol[neighbourhood[i].polymer()];
                                    sb_tmp.negative = mon_tmp;
                                    sb_tmp.centroid_positive = centroid_positive;
                                    sb_tmp.centroid_negative = centroid_negative;
                                    sb_tmp.distance = clipper::Coord_orth::length (centroid_positive, centroid_negative );
                                    sb_tmp.symop_positive = 0;
                                    sb_tmp.symop_negative = neighbourhood[i].symmetry();
                                    results.push_back ( sb_tmp );
                                    break; // do not check the rest of the neighbours because we have what we were looking for
                                } 
                            }
                        }
                    }
                }
                else 
                    break;
            }

    return results;
}


void print_salt_bridges ( std::vector < salt_bridge > salt_bridges, clipper::String ippdb, bool i2mode ) 
{

    std::fstream of_scm; std::fstream of_py;
    of_scm.open("ccontacts-results.scm", std::fstream::out);
    of_py.open ("ccontacts-results.py" , std::fstream::out);
    
    insert_coot_prologue_scheme ( of_scm );
    insert_coot_prologue_python ( of_py );
    insert_coot_files_loadup_scheme (of_scm, ippdb, i2mode );
    insert_coot_files_loadup_python (of_py , ippdb, i2mode );
    
    std::cout << std::endl << "List of detected salt bridges: " << std::endl ;
    std::cout << std::endl << "\tPositive    Distance    Negative";
    std::cout << std::endl << "\t--------------------------------" << std::endl;
    
    int inter, intra;
    inter = 0; intra = 0;

    for ( int i = 0; i < salt_bridges.size() ; i++ )
    {
        clipper::String message = "";

        std::cout << "\t" << salt_bridges[i].positive_chain.id() << "/" << std::fixed << std::setw( 3 ) << salt_bridges[i].positive.id().trim() << " " << salt_bridges[i].positive.type();
        std::cout << std::fixed << std::setw( 4 ) << std::setprecision( 2 ) << " -- [" << salt_bridges[i].distance << "] -- " ;
        std::cout << salt_bridges[i].negative_chain.id() << "/" << std::fixed << std::setw (3) << salt_bridges[i].negative.id().trim() << " " << salt_bridges[i].negative.type();
        
        message = salt_bridges[i].positive_chain.id() + "/" + salt_bridges[i].positive.id().trim() + " " 
                  + salt_bridges[i].positive.type().trim() + " - " + salt_bridges[i].negative_chain.id() 
                  + "/" + salt_bridges[i].negative.id().trim() + " " + salt_bridges[i].negative.type().trim() ;
        
        if ( salt_bridges[i].symop_negative != 0 ) std::cout << "*" ;
        std::cout << std::endl ;
        salt_bridges[i].positive_chain.id() == salt_bridges[i].negative_chain.id() ? intra++ : inter++;

        insert_coot_go_to_res_scheme ( of_scm, salt_bridges[i].centroid_positive, message );
        insert_coot_go_to_res_python ( of_py, salt_bridges[i].centroid_positive, message );
    }

    std::cout << std::endl << "Total: " << salt_bridges.size() << "\tIntramolecular: " << intra << "\tInterfacing: " << inter << std::endl; 

    insert_coot_epilogue_scheme ( of_scm );
    insert_coot_epilogue_python ( of_py );

    return;
}

void insert_coot_prologue_scheme ( std::fstream& output )
{

    output  << "; This script has been created by ccontacts (Jon Agirre, University of York)\n"        
            << "(set-graphics-window-size 1873 968)\n"
            << "(set-graphics-window-position 0 0)\n"
            << "(set-go-to-atom-window-position 0 19)\n"
            << "(set-display-control-dialog-position 366 20)\n"
            << "(vt-surface 2)\n"
            << "(set-clipping-front  0.00)\n"
            << "(set-clipping-back  0.00)\n"
            << "(set-map-radius 10.00)\n"
            << "(set-iso-level-increment  0.0500)\n"
            << "(set-diff-map-iso-level-increment  0.0050)\n"
            << "(set-colour-map-rotation-on-read-pdb 21.00)\n"
            << "(set-colour-map-rotation-on-read-pdb-flag 1)\n"
            << "(set-colour-map-rotation-on-read-pdb-c-only-flag 1)\n"
            << "(set-swap-difference-map-colours 0)\n"
            << "(set-background-colour  0.00  0.00  0.00)\n"
            << "(set-symmetry-size 13.00)\n"
            << "(set-symmetry-colour-merge  0.50)\n"
            << "(set-symmetry-colour  0.10  0.20  0.80)\n"
            << "(set-symmetry-atom-labels-expanded 0)\n"
            << "(set-active-map-drag-flag 1)\n"
            << "(set-show-aniso 0)\n"
            << "(set-aniso-probability 50.00)\n"
            << "(set-smooth-scroll-steps 40)\n"
            << "(set-smooth-scroll-limit 10.00)\n"
            << "(set-font-size 2)\n"
            << "(set-rotation-centre-size  0.10)\n"
            << "(set-default-bond-thickness 5)\n"
            << "(scale-zoom  0.20)\n"
            << "(set-nomenclature-errors-on-read \"auto-correct\")"
            << "(set-run-state-file-status 0)\n";

}


void insert_coot_prologue_python ( std::fstream& output )
{

 output  << "# This script has been created by ccontacts (Jon Agirre, University of York)\n"           
            << "set_graphics_window_size (1873, 968)\n"
            << "set_graphics_window_position (0, 0)\n"
            << "set_go_to_atom_window_position (0, 19)\n"
            << "vt_surface (2)\n"
            << "set_clipping_front  (0.00)\n"
            << "set_clipping_back  (0.00)\n"
            << "set_map_radius (10.00)\n"
            << "set_iso_level_increment  (0.0500)\n"
            << "set_diff_map_iso_level_increment  (0.0050)\n"
            << "set_colour_map_rotation_on_read_pdb (21.00)\n"
            << "set_colour_map_rotation_on_read_pdb_flag (1)\n"
            << "set_colour_map_rotation_on_read_pdb_c_only_flag (1)\n"
            << "set_swap_difference_map_colours (0)\n"
            << "set_background_colour  (0.00,  0.00,  0.00)\n"
            << "set_symmetry_size (13.00)\n"
            << "set_symmetry_colour_merge  (0.50)\n"
            << "set_symmetry_colour  (0.10,  0.20,  0.80)\n"
            << "set_symmetry_atom_labels_expanded (0)\n"
            << "set_active_map_drag_flag (1)\n"
            << "set_show_aniso (0)\n"
            << "set_aniso_probability (50.00)\n"
            << "set_smooth_scroll_steps (40)\n"
            << "set_smooth_scroll_limit (10.00)\n"
            << "set_font_size (2)\n"
            << "set_rotation_centre_size (0.10)\n"
            << "set_default_bond_thickness (4)\n"
            << "scale_zoom (0.20)\n"
            << "set_nomenclature_errors_on_read (\"auto-correct\")\n"
            << "set_run_state_file_status (0)\n"
            << "toggle_idle_spin_function\n";   

}


void insert_coot_files_loadup_scheme ( std::fstream& output, const clipper::String& pdb, bool mode )
{
    if (!mode) 
        output << "(handle-read-draw-molecule \"" << pdb << "\")\n";

    output << "(interesting-things-gui \"Contact list\"\n\t(list\n\t\t";
}


void insert_coot_files_loadup_python ( std::fstream& output, const clipper::String& pdb, bool mode )
{
    if (!mode) 
        output  << "handle_read_draw_molecule (\"" << pdb << "\")\n";
        
    output << "interesting_things_gui (\"Contact list\",[\n";
}


void insert_coot_go_to_res_scheme ( std::fstream& output, const clipper::Coord_orth& res_centre, const clipper::String& text )
{
    output  << "\t(list\t\"" << text << "\"\t" << res_centre.x() << "\t" << res_centre.y() << "\t" << res_centre.z() << ")\n";
}


void insert_coot_go_to_res_python ( std::fstream& output, const clipper::Coord_orth& res_centre, const clipper::String& text )
{
    output  << "\t[\"" << text << "\",\t" << res_centre.x() << ",\t" << res_centre.y() << ",\t" << res_centre.z() << "],\n";
}


void insert_coot_epilogue_scheme ( std::fstream& output )
{
    output  << "\n\n))\n(set-scroll-wheel-map 3)\n"
            << "(set-matrix 60.00)\n"
            << "(set-show-symmetry-master 0)\n";
}


void insert_coot_epilogue_python ( std::fstream& output )
{
    output  << "\n\n])\nset_scroll_wheel_map (3)\n"
            << "set_matrix (60.00)\n"
            << "set_show_symmetry_master (0)\n";
}
