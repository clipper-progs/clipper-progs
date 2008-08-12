// clipper CNS->MTZ utility
/* (C) 2007 Kevin Cowtan */

#include <iostream>
#include <fstream>
#include <set>
#include <stdlib.h>

#include "convert2mtz.h"

int main( int argc, char** argv )
{
  CCP4Program prog( "convert2mtz", "0.4", "$Date: 2007/03/30" );

  // defaults
  clipper::String title;
  clipper::String ipfile = "NONE";
  clipper::String opfile = "NONE";
  clipper::String ccell = "NONE";
  clipper::String cspgr = "NONE";
  clipper::String ipfilepdb = "NONE";
  clipper::String colin = "NONE";
  bool complete = true;
  bool anom = false;
  int seed = 54321;
  int verbose = 0;

  // command input
  CCP4CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if ( args[arg] == "-title" ) {
      if ( ++arg < args.size() ) title = args[arg];
    } else if ( args[arg] == "-hklin" || args[arg] == "-cnsin" ) {
      if ( ++arg < args.size() ) ipfile = args[arg];
    } else if ( args[arg] == "-mtzout" ) {
      if ( ++arg < args.size() ) opfile = args[arg];
    } else if ( args[arg] == "-cell" ) {
      if ( ++arg < args.size() ) ccell = args[arg];
    } else if ( args[arg] == "-spacegroup" ) {
      if ( ++arg < args.size() ) cspgr = args[arg];
    } else if ( args[arg] == "-pdbin" ) {
      if ( ++arg < args.size() ) ipfilepdb = args[arg];
    } else if ( args[arg] == "-colin" ) {
      if ( ++arg < args.size() ) colin = args[arg];
    } else if ( args[arg] == "-anomalous" ) {
      anom = true;
    } else if ( args[arg] == "-no-complete" ) {
      complete = false;
    } else if ( args[arg] == "-seed" ) {
      if ( ++arg < args.size() ) seed = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-verbose" ) {
      if ( ++arg < args.size() ) verbose = clipper::String(args[arg]).i();
    } else {
      std::cout << "Unrecognized:\t" << args[arg] << "\n";
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "Usage: convert2mtz\n\t-hklin <filename> [COMPULSORY]\n\t-mtzout <filename>\n\t-cell a,b,c,a,b,g\n\t-spacegroup <spacegroup>\n\t-pdbin <filename>\n\t-anomalous\n\t-no-complete\n\t-colin <column names>\n\t-seed <seed>\nConvert formatted reflection file to MTZ.\nCell and spacegroup may be specified directly or by providing a PDB file.\nFor CNS files, column labels and anomalous are detected automatically.\nFor other files (e.g. Shelx/XtalView), give -colin.\n";
    exit(1);
  }

  if ( ipfilepdb == "NONE" && ccell == "NONE" )
    { std::cout << "Missing cell.";       exit(1); }
  if ( ipfilepdb == "NONE" && cspgr == "NONE" )
    { std::cout << "Missing spacegroup."; exit(1); }

  // set defaults
  srand( seed );
  if ( opfile == "NONE" ) {
    opfile = ipfile;
    int ldot = int(opfile.rfind("."));
    if ( ldot != int(std::string::npos) ) opfile = opfile.substr(0,ldot);
    opfile = opfile + ".mtz";
  }

  // set crystal info
  clipper::Cell cell;
  clipper::Spacegroup spgr;
  clipper::Resolution reso;
  if ( ipfilepdb != "NONE" ) {
    clipper::MMDBManager mmdb;
    mmdb.ReadCoorFile( (char*)ipfilepdb.c_str() );
    spgr = mmdb.spacegroup();
    cell = mmdb.cell();
  } else {
    double cd[] = { 100.0, 100.0, 100.0, 90.0, 90.0, 90.0 };
    std::vector<clipper::String> cellx = ccell.split( ", " );
    for ( int i = 0; i < cellx.size(); i++ ) cd[i] = cellx[i].f();
    clipper::Cell_descr celld( cd[0], cd[1], cd[2], cd[3], cd[4], cd[5] );
    clipper::Spgr_descr spgrd( cspgr );
    cell = clipper::Cell( celld );
    spgr = clipper::Spacegroup( spgrd );
  }  

  // read input file
  std::vector<clipper::String> words;
  std::ifstream file( ipfile.c_str() );
  int lvl = 0;
  std::string word;
  while ( file.good() ) {
    char c;
    file.get( c );
    if ( c == '{' ) lvl++;
    if ( lvl == 0 ) {
      if ( isspace(c) || c == '=' ) {
	if ( word.length() != 0 ) {
	  words.push_back( word );
	  word = "";
	}
      } else {
	word = word + c;
      }
    }
    if ( c == '}' ) lvl--;
  }
  file.close();

  // file info
  int nref = 0;
  std::vector<std::string> cols, typs;
  std::vector<std::pair<std::string,std::vector<std::string> > > grps;
  cols.push_back( "H" ); cols.push_back( "K" ); cols.push_back( "L" );
  typs.push_back( "H" ); typs.push_back( "H" ); typs.push_back( "H" );
  std::vector<float> vals;

  // loop over words and read all info from file
  int pos = 0;
  // header section
  while ( pos < words.size() ) {
    // CNS headers (irrelevent for non-CNS files)
    if        ( words[pos].substr(0,4) == "NREF" ) {
      if ( ++pos < words.size() ) nref = clipper::String(words[pos]).i();
    } else if ( words[pos].substr(0,4) == "ANOM" ) {
      if ( ++pos < words.size() ) anom = ( words[pos][0] == 'T' );
    } else if ( words[pos].substr(0,4) == "HERM" ) {
      if ( ++pos < words.size() ) anom = ( words[pos][0] == 'F' );
    } else if ( words[pos].substr(0,4) == "DECL" ) {
      // COLUMN DECLARATIONS
      ++pos;
      std::string col, dom, typ;
      while ( pos < words.size() ) {
	if ( words[pos] == "END" ) {
	  if ( dom.substr(0,4) == "RECI" ) {
	    if ( typ.substr(0,4) == "INTE" ) {
	      cols.push_back( col );
	      typs.push_back( "I" );
	    } else if ( typ.substr(0,4) == "REAL" ) {
	      cols.push_back( col );
	      typs.push_back( "R" );
	    } else if ( typ.substr(0,4) == "COMP" ) {
	      cols.push_back( col );
	      typs.push_back( "F" );
	      cols.push_back( col+".phase" );
	      typs.push_back( "P" );
	    }
	  }
	  break;
	} else if ( words[pos].substr(0,4) == "NAME" ) {
	  if ( ++pos < words.size() ) col = words[pos];
	} else if ( words[pos].substr(0,4) == "DOMA" ) {
	  if ( ++pos < words.size() ) dom = words[pos];
	} else if ( words[pos].substr(0,4) == "TYPE" ) {
	  if ( ++pos < words.size() ) typ = words[pos];
	}
	++pos;
      }
    } else if ( words[pos].substr(0,4) == "GROU" ) {
      // GROUP DECLARATIONS
      ++pos;
      std::string gtyp;
      std::vector<std::string> gcol;
      while ( pos < words.size() ) {
	if ( words[pos] == "END" ) {
	  grps.push_back( std::pair<std::string,std::vector<std::string> >( gtyp, gcol ) );
	  break;
	} else if ( words[pos].substr(0,4) == "TYPE" ) {
	  if ( ++pos < words.size() ) gtyp = words[pos];
	} else if ( words[pos].substr(0,4) == "OBJE" ) {
	  if ( ++pos < words.size() ) gcol.push_back( words[pos] );
	}
	++pos;
      }
    } else {
      if ( words[pos].find_first_not_of("0123456789.-") == std::string::npos ) {
	vals.push_back( clipper::String(words[pos]).f() );
      }
    }
    ++pos;
  }

  // Now do group assignments
  for ( int g = 0; g < grps.size(); g++ ) {
    if ( grps[g].first == "HL" ) {
      for ( int c = 0; c < grps[g].second.size(); c++ ) {
	for ( int i = 0; i < cols.size(); i++ ) {
	  if ( cols[i] == grps[g].second[c] )
	    typs[i] = "A";
	}
      }
    }
  }

  // if columns were input, they override the headers
  if ( colin != "NONE" ) {
    cols.clear(); typs.clear();
    cols.push_back( "H" ); cols.push_back( "K" ); cols.push_back( "L" );
    typs.push_back( "H" ); typs.push_back( "H" ); typs.push_back( "H" );
    std::vector<clipper::String> cin = colin.split(" ,");
    for ( int i = 0; i < cin.size(); i++ ) {
      cols.push_back( cin[i] );
      typs.push_back( "R" );
    }
    nref = vals.size() / cols.size();
  }

  // check data length
  int ncol = cols.size();
  if ( ncol*nref != vals.size() ) {
    std::cout << "Error: data length does not match number of reflections.\n";
    exit(1);
  }

  // column assignment heuristics
  for ( int i = 3; i < ncol; i++ ) {
    std::string c;
    for ( int j = 0; j < cols[i].length(); j++ ) c += toupper( cols[i][j] );
    if ( typs[i] == "R" ) {
      if      ( cols[i].find("SIG") != std::string::npos )
	typs[i] = "Q";
      else if ( cols[i].find("FOM") != std::string::npos )
	typs[i] = "W";
      else if ( cols[i].find("HLA") != std::string::npos )
	typs[i] = "A";
      else if ( cols[i].find("FREE")!= std::string::npos )
	typs[i] = "I";
      else if ( cols[i][0] == 'F' )
	typs[i] = "F";
      else if ( cols[i][0] == 'E' )
	typs[i] = "E";
      else if ( cols[i][0] == 'I' )
	typs[i] = "J";
      else if ( cols[i][0] == 'P' )
	typs[i] = "P";
      else if ( cols[i][0] == 'S' )
	typs[i] = "Q";
      else if ( cols[i][0] == 'W' )	   
	typs[i] = "W";
      else if ( cols[i][0] == 'M' )	   
	typs[i] = "W";
      else if ( cols[i][0] == 'A' )	   
	typs[i] = "A";
      else if ( cols[i][0] == 'B' )	   
	typs[i] = "A";
      else if ( cols[i][0] == 'C' )	   
	typs[i] = "A";
      else if ( cols[i][0] == 'D' )	   
	typs[i] = "A";
      else if ( cols[i][0] == 'H' )	   
	typs[i] = "A";
    }
  }

  // post-assignment diagnostics
  std::cout << std::endl << "Columns found:" << std::endl;
  for ( int i = 0; i < cols.size(); i++ ) {
    std::cout << " " << i << " " << cols[i] << " " << typs[i] << std::endl;
  }
  std::cout << "Note: these do not represent the final MTZ column types" << std::endl << std::endl;

  // get reflections
  std::set<clipper::HKL,HKLlessthan> hkl_set;
  const clipper::HKL zero( 0, 0, 0 );
  for ( int i = 0; i < nref; i++ ) {
    clipper::HKL hkl( Util::intr(vals[i*ncol  ]),
		      Util::intr(vals[i*ncol+1]),
		      Util::intr(vals[i*ncol+2]) );
    clipper::HKL asu;
    if ( hkl != zero ) {
      for ( int s = 0; s < spgr.num_primops(); s++ ) {
	asu = hkl.transform( spgr.symop(s) );
	if ( spgr.recip_asu( asu ) ) break;
	asu = -asu;
	if ( spgr.recip_asu( asu ) ) break;
      }
      hkl_set.insert( asu );
    }
  }
  std::vector<clipper::HKL> hkl_list( hkl_set.begin(), hkl_set.end() );
  // get resolution
  double smax = 0.0;
  for ( int i = 0; i < hkl_list.size(); i++ ) {
    double s = hkl_list[i].invresolsq( cell );
    if ( s > smax ) smax = s;
  }
  reso = clipper::Resolution( 0.999999/sqrt(smax) );

  // prepare data arrays
  clipper::HKL_info hkls( spgr, cell, reso );
  if ( complete ) {
    hkl_list.clear();
    clipper::HKL_info hkl0( spgr, cell, reso, true );
    for ( int i = 0; i < hkl0.num_reflections(); i++ )
      if ( hkl0.hkl_of(i) != zero )
	hkl_list.push_back( hkl0.hkl_of(i) );
  }
  hkls.add_hkl_list( hkl_list );
  std::vector<clipper::String> opcols;
  std::vector<clipper::HKL_data<dataI> > datai;
  std::vector<clipper::HKL_data<dataF> > dataf;
  std::vector<clipper::HKL_data<dataIano> > dataiano;
  std::vector<clipper::HKL_data<dataFano> > datafano;
  std::vector<clipper::HKL_data<dataSigIano> > datasigiano;
  std::vector<clipper::HKL_data<dataSigFano> > datasigfano;
  std::vector<clipper::HKL_data<dataE> > datae;
  std::vector<clipper::HKL_data<dataSig> > datasig;
  std::vector<clipper::HKL_data<dataPhi> > dataphi;
  std::vector<clipper::HKL_data<dataFom> > datafom;
  std::vector<clipper::HKL_data<dataABCD> > dataabcd;
  std::vector<clipper::HKL_data<dataFlag> > dataflag;
  clipper::HKL_data<clipper::data32::Flag> datafree(hkls);

  // fill data arrays
  // loop over the columns
  bool friedel;
  int isym;
  for ( int col = 3; col < ncol; col++ ) {
    // I or Iano
    if ( typs[col] == "J" ) {
      if ( !anom ) {
	clipper::HKL_data<dataI> data( hkls );
	for ( int ref = 0; ref < nref; ref++ ) {
	  clipper::HKL hkl( Util::intr(vals[ref*ncol  ]),
			    Util::intr(vals[ref*ncol+1]),
			    Util::intr(vals[ref*ncol+2]) );
	  dataI dat( vals[ref*ncol+col] );
	  data.set_data( hkl, dat );
	}
	datai.push_back( data );
	opcols.push_back("i"+cols[col]);
      } else {
	clipper::HKL_data<dataIano> data( hkls );
	for ( int ref = 0; ref < nref; ref++ ) {
	  clipper::HKL hkl( Util::intr(vals[ref*ncol  ]),
			    Util::intr(vals[ref*ncol+1]),
			    Util::intr(vals[ref*ncol+2]) );
	  clipper::HKL hkl1 = hkls.find_sym( hkl, isym, friedel );
	  int index = hkls.index_of( hkl1 );
	  float v = vals[ref*ncol+col];
	  if ( friedel ) data[index] = dataIano( data[index].I_pl(), v );
          else           data[index] = dataIano( v, data[index].I_mi() );
	}
	dataiano.push_back( data );
	opcols.push_back("I"+cols[col]+"+,"+cols[col]+"-");
      }
    }
    // F or Fano
    if ( typs[col] == "F" ) {
      if ( !anom ) {
	clipper::HKL_data<dataF> data( hkls );
	for ( int ref = 0; ref < nref; ref++ ) {
	  clipper::HKL hkl( Util::intr(vals[ref*ncol  ]),
			    Util::intr(vals[ref*ncol+1]),
			    Util::intr(vals[ref*ncol+2]) );
	  dataF dat( vals[ref*ncol+col] );
	  data.set_data( hkl, dat );
	}
	dataf.push_back( data );
	opcols.push_back("f"+cols[col]);
      } else {
	clipper::HKL_data<dataFano> data( hkls );
	for ( int ref = 0; ref < nref; ref++ ) {
	  clipper::HKL hkl( Util::intr(vals[ref*ncol  ]),
			    Util::intr(vals[ref*ncol+1]),
			    Util::intr(vals[ref*ncol+2]) );
	  clipper::HKL hkl1 = hkls.find_sym( hkl, isym, friedel );
	  int index = hkls.index_of( hkl1 );
	  float v = vals[ref*ncol+col];
	  if ( friedel ) data[index] = dataFano( data[index].F_pl(), v );
          else           data[index] = dataFano( v, data[index].F_mi() );
	}
	datafano.push_back( data );
	opcols.push_back("F"+cols[col]+"+,"+cols[col]+"-");
      }
    }
    // sigF or sigFano or sigIano
    if ( typs[col] == "Q" ) {
      if ( !anom ) {
	clipper::HKL_data<dataSig> data( hkls );
	for ( int ref = 0; ref < nref; ref++ ) {
	  clipper::HKL hkl( Util::intr(vals[ref*ncol  ]),
			    Util::intr(vals[ref*ncol+1]),
			    Util::intr(vals[ref*ncol+2]) );
	  dataSig dat( vals[ref*ncol+col] );
	  data.set_data( hkl, dat );
	}
	datasig.push_back( data );
	opcols.push_back("S"+cols[col]);
      } else {
	if ( typs[col-1] == "J" ) {
	  clipper::HKL_data<dataSigIano> data( hkls );
	  for ( int ref = 0; ref < nref; ref++ ) {
	    clipper::HKL hkl( Util::intr(vals[ref*ncol  ]),
			      Util::intr(vals[ref*ncol+1]),
			      Util::intr(vals[ref*ncol+2]) );
	    clipper::HKL hkl1 = hkls.find_sym( hkl, isym, friedel );
	    int index = hkls.index_of( hkl1 );
	    float v = vals[ref*ncol+col];
	    if ( friedel ) data[index] = dataSigIano(data[index].sigI_pl(), v);
	    else           data[index] = dataSigIano(v, data[index].sigI_mi());
	  }
	  datasigiano.push_back( data );
	  opcols.push_back("1"+cols[col]+"+,"+cols[col]+"-");
	} else {
	  clipper::HKL_data<dataSigFano> data( hkls );
	  for ( int ref = 0; ref < nref; ref++ ) {
	    clipper::HKL hkl( Util::intr(vals[ref*ncol  ]),
			      Util::intr(vals[ref*ncol+1]),
			      Util::intr(vals[ref*ncol+2]) );
	    clipper::HKL hkl1 = hkls.find_sym( hkl, isym, friedel );
	    int index = hkls.index_of( hkl1 );
	    float v = vals[ref*ncol+col];
	    if ( friedel ) data[index] = dataSigFano(data[index].sigF_pl(), v);
	    else           data[index] = dataSigFano(v, data[index].sigF_mi());
	  }
	  datasigfano.push_back( data );
	  opcols.push_back("2"+cols[col]+"+,"+cols[col]+"-");
	}
      }
    }
    // E
    if ( typs[col] == "E" ) {
      clipper::HKL_data<dataE> data( hkls );
      for ( int ref = 0; ref < nref; ref++ ) {
	clipper::HKL hkl( Util::intr(vals[ref*ncol  ]),
			  Util::intr(vals[ref*ncol+1]),
			  Util::intr(vals[ref*ncol+2]) );
	dataE dat( vals[ref*ncol+col] );
	data.set_data( hkl, dat );
      }
      datae.push_back( data );
      opcols.push_back("E"+cols[col]);
    }
    // Phase
    if ( typs[col] == "P" ) {
      clipper::HKL_data<dataPhi> data( hkls );
      for ( int ref = 0; ref < nref; ref++ ) {
	clipper::HKL hkl( Util::intr(vals[ref*ncol  ]),
			  Util::intr(vals[ref*ncol+1]),
			  Util::intr(vals[ref*ncol+2]) );
	dataPhi dat( Util::d2rad(vals[ref*ncol+col]) );
	data.set_data( hkl, dat );
      }
      dataphi.push_back( data );
      opcols.push_back("P"+cols[col]);
    }
    // FOM
    if ( typs[col] == "W" ) {
      clipper::HKL_data<dataFom> data( hkls );
      for ( int ref = 0; ref < nref; ref++ ) {
	clipper::HKL hkl( Util::intr(vals[ref*ncol  ]),
			  Util::intr(vals[ref*ncol+1]),
			  Util::intr(vals[ref*ncol+2]) );
	dataFom dat( vals[ref*ncol+col] );
	data.set_data( hkl, dat );
      }
      datafom.push_back( data );
      opcols.push_back("W"+cols[col]);
    }
    // Flag
    if ( typs[col] == "I" ) {
      clipper::HKL_data<dataFlag> data( hkls );
      for ( int ref = 0; ref < nref; ref++ ) {
	clipper::HKL hkl( Util::intr(vals[ref*ncol  ]),
			  Util::intr(vals[ref*ncol+1]),
			  Util::intr(vals[ref*ncol+2]) );
	dataFlag dat( Util::max(Util::intr(vals[ref*ncol+col]), 0 ) );
	data.set_data( hkl, dat );
      }
      dataflag.push_back( data );
      opcols.push_back("X"+cols[col]);
    }
    // HL coeffs
    if ( typs[col] == "A" && col < ncol-3 ) {
      clipper::HKL_data<dataABCD> data( hkls );
      for ( int ref = 0; ref < nref; ref++ ) {
	clipper::HKL hkl( Util::intr(vals[ref*ncol  ]),
			  Util::intr(vals[ref*ncol+1]),
			  Util::intr(vals[ref*ncol+2]) );
	dataABCD dat( vals[ref*ncol+col  ], vals[ref*ncol+col+1],
		      vals[ref*ncol+col+2], vals[ref*ncol+col+3] );
	data.set_data( hkl, dat );
      }
      dataabcd.push_back( data );
      opcols.push_back("A"+cols[col]+","+cols[col+1]+","+cols[col+2]+","+cols[col+3]);
      col += 3;
    }
  }

  // Check for Fano/Iano without sigFano/sigIano
  if ( dataiano.size() > 0 && datasigiano.size() == 0 ) {
    clipper::HKL_data<dataSigIano> data( hkls );
    dataSigIano dat( 1.0, 1.0 );
    data = dat;
    datasigiano.push_back( data );
    opcols.push_back("1SIGI+,SIGI-");
  }
  if ( datafano.size() > 0 && datasigfano.size() == 0 ) {
    clipper::HKL_data<dataSigFano> data( hkls );
    dataSigFano dat( 1.0, 1.0 );
    data = dat;
    datasigfano.push_back( data );
    opcols.push_back("2SIGF+,SIGF-");
  }

  // Free-R flag conversion
  typedef clipper::HKL_info::HKL_reference_index HRI;
  if ( dataflag.size() > 0 ) {
    int n, n0, n1, fl, nfl;
    n = n0 = n1 = 0;
    for ( HRI ih = dataflag[0].first(); !ih.last(); ih.next() ) {
      if ( dataflag[0][ih].flag() >= 0 ) n++;
      if ( dataflag[0][ih].flag() == 0 ) n0++;
      if ( dataflag[0][ih].flag() == 1 ) n1++;
    }
    fl = 1; nfl = n1;  // conventional CNS flag definition
    if ( n0 < n1/2 ) { fl = 0; nfl = n0; }  // SHELX flag definition
    // generate CCP4 free flags
    int nfree = Util::intr( double(n) / double(nfl) );
    for ( HRI ih = dataflag[0].first(); !ih.last(); ih.next() )
      if ( dataflag[0][ih].flag() == fl )  // free: add to free set
	datafree[ih] = dataFlag( 0 );
      else                           // missing or work: pick a working set
	datafree[ih] = dataFlag( rand()%(nfree-1) + 1 );
  } else {
    int nfree = nref / 1000;
    if ( nfree < 10 ) nfree = 10;
    if ( nfree > 20 ) nfree = 20;
    for ( HRI ih = datafree.first(); !ih.last(); ih.next() )
      datafree[ih] = dataFlag( rand()%nfree );
  }

  // tidy phases
  for ( int i = 0; i < dataphi.size(); i++ )
    for ( HRI ih = dataphi[i].first(); !ih.last(); ih.next() )
      dataphi[i][ih].phi() = clipper::Util::mod( dataphi[i][ih].phi(),
                                               clipper::Util::twopi() );

  // export
  clipper::CCP4MTZfile mtzout;
  mtzout.open_write( opfile );
  mtzout.export_hkl_info( hkls );
  for ( int i = 0; i < opcols.size(); i++ ) {
    // game name
    clipper::String ty = opcols[i].substr(0,1);
    clipper::String nm = "/*/*/["+opcols[i].substr(1)+"]";
    // count which column to export
    int pos = 0;
    for ( int j = 0; j < i; j++ )
      if ( opcols[i][0] == opcols[j][0] ) pos++;
    // and export
    if ( ty == "I" ) mtzout.export_hkl_data( dataiano[pos],    nm );
    if ( ty == "1" ) mtzout.export_hkl_data( datasigiano[pos], nm );
    if ( ty == "F" ) mtzout.export_hkl_data( datafano[pos],    nm );
    if ( ty == "2" ) mtzout.export_hkl_data( datasigfano[pos], nm );
    if ( ty == "i" ) mtzout.export_hkl_data( datai[pos],       nm );
    if ( ty == "f" ) mtzout.export_hkl_data( dataf[pos],       nm );
    if ( ty == "E" ) mtzout.export_hkl_data( datae[pos],       nm );
    if ( ty == "S" ) mtzout.export_hkl_data( datasig[pos],     nm );
    if ( ty == "P" ) mtzout.export_hkl_data( dataphi[pos],     nm );
    if ( ty == "W" ) mtzout.export_hkl_data( datafom[pos],     nm );
    if ( ty == "X" ) mtzout.export_hkl_data( dataflag[pos],    nm );
    if ( ty == "A" ) mtzout.export_hkl_data( dataabcd[pos],    nm );
  }
  mtzout.export_hkl_data( datafree, "/*/*/[FreeF_flag]" );
  mtzout.close_write();

  // output
  std::cout << ( anom ? "Collecting HKL and -HKL" : "Not collecting HKL and -HKL - anomalous data will be lost." ) << std::endl;
  std::cout << ( complete ? "Inserting missing reflections" : "Not inserting missing reflections" ) << std::endl;
  std::cout << "Number of reflections input : " << nref << std::endl;
  std::cout << "Number of reflections output: " << hkls.num_reflections() << std::endl;
}
