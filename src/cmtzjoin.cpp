// Clipper app to join mtzs
/* Copyright 2012 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>

#include <set>

struct Compare_HKL{ bool operator() ( const clipper::HKL& h1, const clipper::HKL& h2 ) const { return ( h1.h()<h2.h() || ( h1.h()==h2.h() && ( h1.k()<h2.k() || ( h1.k()==h2.k() && h1.l()<h2.l() ) ) ) ); } };

    typedef float dtype;
    using clipper::Util;
    using clipper::Datatype_base; using clipper::String;
    using clipper::ftype; using clipper::xtype;

    class F_sigF_anom : private Datatype_base
    { public:
      F_sigF_anom() { set_null(); }
      void set_null() { Util::set_null(f_pl_); Util::set_null(f_mi_); Util::set_null(sigf_pl_); Util::set_null(sigf_mi_); }
      static String type() { return "F_sigF_anom"; }
      void friedel() { dtype f=f_pl_; f_pl_=f_mi_; f_mi_=f;
                       f=sigf_pl_; sigf_pl_=sigf_mi_; sigf_mi_=f; }
      void shift_phase(const ftype&) {}
      bool missing() const { return (Util::is_nan(f_pl_) && Util::is_nan(f_mi_)); }
      static int data_size() { return 4; }
      static String data_names() { return "F+ sigF+ F- sigF-"; }
      void data_export( xtype a[] ) const { a[0] = f_pl(); a[1] = sigf_pl(); a[2] = f_mi(); a[3] = sigf_mi(); }
      void data_import( const xtype a[] ) { f_pl() = a[0]; sigf_pl() = a[1]; f_mi() = a[2]; sigf_mi() = a[3]; }
      //! this type is scalable - apply magnitude scale factor
      void scale(const ftype& s) { f_pl_ *= s; sigf_pl_ *= s; f_mi_ *= s; sigf_mi_ *= s; }
      // accessors
      const dtype& f_pl() const { return f_pl_; }  //<! read access
      const dtype& sigf_pl() const { return sigf_pl_; }  //<! read access
      const dtype& f_mi() const { return f_mi_; }  //<! read access
      const dtype& sigf_mi() const { return sigf_mi_; }  //<! read access
      dtype& f_pl() { return f_pl_; }  //<! write access
      dtype& sigf_pl() { return sigf_pl_; }  //<! write access
      dtype& f_mi() { return f_mi_; }  //<! write access
      dtype& sigf_mi() { return sigf_mi_; }  //<! write access
    private:
      dtype f_pl_, f_mi_, sigf_pl_, sigf_mi_;
    };

    class I_sigI_anom : private Datatype_base
    { public:
      I_sigI_anom() { set_null(); }
      void set_null() { Util::set_null(I_pl_); Util::set_null(I_mi_); Util::set_null(sigI_pl_); Util::set_null(sigI_mi_); }
      static String type() { return "I_sigI_anom"; }
      void friedel() { dtype I=I_pl_; I_pl_=I_mi_; I_mi_=I;
                       I=sigI_pl_; sigI_pl_=sigI_mi_; sigI_mi_=I; }
      void shift_phase(const ftype&) {}
      bool missing() const { return (Util::is_nan(I_pl_) && Util::is_nan(I_mi_)); }
      static int data_size() { return 4; }
      static String data_names() { return "I+ sigI+ I- sigI-"; }
      void data_export( xtype a[] ) const { a[0] = I_pl(); a[1] = sigI_pl(); a[2] = I_mi(); a[3] = sigI_mi(); }
      void data_import( const xtype a[] ) { I_pl() = a[0]; sigI_pl() = a[1]; I_mi() = a[2]; sigI_mi() = a[3]; }
      //! this type is scalable - apply magnitude scale factor
      void scale(const ftype& s) { I_pl_ *= (s*s); sigI_pl_ *= (s*s); I_mi_ *= (s*s); sigI_mi_ *= (s*s); }
      // accessors
      const dtype& I_pl() const { return I_pl_; }  //<! read access
      const dtype& sigI_pl() const { return sigI_pl_; }  //<! read access
      const dtype& I_mi() const { return I_mi_; }  //<! read access
      const dtype& sigI_mi() const { return sigI_mi_; }  //<! read access
      dtype& I_pl() { return I_pl_; }  //<! write access
      dtype& sigI_pl() { return sigI_pl_; }  //<! write access
      dtype& I_mi() { return I_mi_; }  //<! write access
      dtype& sigI_mi() { return sigI_mi_; }  //<! write access
    private:
      dtype I_pl_, I_mi_, sigI_pl_, sigI_mi_;
    };

int main( int argc, char** argv )
{
  CCP4Program prog( "cmtzjoin", "0.1", "$Date: 2012/05/01" );

  // defaults
  clipper::String title;
  clipper::String opfile;
  std::vector<clipper::String> ipfiles, ipcols, opcols;

  // command input
  CCP4CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if ( args[arg] == "-title" ) {
      if ( ++arg < args.size() ) title = args[arg];
    } else if ( args[arg] == "-mtzout" ) {
      if ( ++arg < args.size() ) opfile = args[arg];
    } else if ( args[arg] == "-mtzin" ) {
      if ( ++arg < args.size() ) ipfiles.push_back( args[arg] );
    } else if ( args[arg] == "-colin" ) {
      if ( ++arg < args.size() ) ipcols.push_back( args[arg] );
    } else if ( args[arg] == "-colout" ) {
      if ( ++arg < args.size() ) opcols.push_back( args[arg] );
    } else {
      std::cout << "Unrecognized:\t" << args[arg] << "\n";
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "Usage: cmtzjoin\n\t-mtzin <filename>\n\n";
    return 1;
  }

  // get spacegroup/cell/HKLs
  std::set<clipper::HKL,Compare_HKL> hkl_set;
  clipper::Spacegroup spgr = clipper::Spacegroup::null();
  clipper::Cell cell;
  double reslim = 1.0e6;

  using clipper::HKL_info;        using clipper::HKL_data;
  using clipper::data32::F_sigF;  using clipper::data32::I_sigI;
  using clipper::data32::Phi_fom; using clipper::data32::ABCD;
  using clipper::data32::F_phi;   using clipper::data32::Flag;
  clipper::CCP4MTZ_type_registry::add_group( "F_sigF_anom", "FANM" );
  clipper::CCP4MTZ_type_registry::add_group( "I_sigI_anom", "IANM" );

  // first pass read to get crystal data and reflection list
  unsigned int nullhash(12345678);
  unsigned int sgrphash(nullhash), pgrphash(nullhash);
  std::vector<std::vector<clipper::String> > cols_all;
  std::vector<clipper::String> types_all;
  for ( int i = 0; i < ipfiles.size(); i++ ) {
    clipper::CCP4MTZfile mtzin;
    mtzin.set_column_label_mode( clipper::CCP4MTZfile::Legacy );
    mtzin.open_read( ipfiles[i] );
    std::vector<clipper::String> mtzcols = mtzin.column_paths();
    // get column names and types
    std::vector<clipper::String> cols;
    clipper::String types;
    if ( i < ipcols.size() ) {  // column names specified
      cols = ipcols[i].split(",");
      for ( int c1 = 0; c1 < cols.size(); c1++ ) {
	for ( int c2 = 0; c2 < mtzcols.size(); c2++ ) {
	  std::vector<clipper::String> nametype = mtzcols[c2].split( " " );
	  clipper::String name = nametype[0].tail();
	  clipper::String type = nametype[1];
	  if ( name == cols[c1] ) types += type;
	}
      }
    } else {  // column names not specified
      for ( int c = 0; c < mtzcols.size(); c++ ) {
	std::vector<clipper::String> nametype = mtzcols[c].split( " " );
        if ( nametype[1] != "H" ) {
	  cols.push_back( nametype[0].tail() );
	  types += nametype[1];
	}
      }
      clipper::String ipnames = cols[0];
      for ( int c = 1; c < cols.size(); c++ ) ipnames += ','+cols[c];
      ipcols.push_back( ipnames );
    }
    cols_all.push_back( cols );
    types_all.push_back( types );
    // check spacegroup
    bool isphase = (types.find("P")!=std::string::npos) || (types.find("A")!=std::string::npos);
    unsigned int shash = mtzin.spacegroup().hash();
    unsigned int phash = mtzin.spacegroup().generator_ops().pgrp_ops().hash();
    if ( pgrphash != nullhash && pgrphash != phash ) {
      std::cout << "Pointgroup mismatch in " << ipfiles[i] << std::endl;
      std::cerr << "Pointgroup mismatch in " << ipfiles[i] << std::endl;
      return 2;
    }
    if ( isphase && sgrphash != nullhash && sgrphash != shash ) {
      std::cout << "Spacegroup mismatch in " << ipfiles[i] << std::endl;
      std::cerr << "Spacegroup mismatch in " << ipfiles[i] << std::endl;
      return 2;
    }
    if ( pgrphash == nullhash ) pgrphash = phash;
    if ( isphase && sgrphash == nullhash ) sgrphash = shash;
    if ( cell.is_null() ) cell = mtzin.cell();
    if ( spgr.is_null() ) spgr = mtzin.spacegroup();
    if ( isphase )        spgr = mtzin.spacegroup();
    reslim = std::min( reslim, mtzin.resolution().limit() );
    // read reflections
    clipper::HKL_info hkls;
    mtzin.import_hkl_info( hkls );
    for ( int h = 0; h < hkls.num_reflections(); h++ )
      hkl_set.insert( hkls.hkl_of(h) );
    mtzin.close_read();
    //std::cerr << i << " " << hkl_set.size() << std::endl;
  }
  std::vector<clipper::HKL> hkl_list( hkl_set.begin(), hkl_set.end() );
  std::cerr << hkl_list.size() << std::endl;
  clipper::Resolution reso( reslim );
  clipper::HKL_info hklinfo( spgr, cell, reso, false );
  hklinfo.add_hkl_list( hkl_list );

  //for ( int h = 0; h < hklinfo.num_reflections(); h++ )
  //  std::cout << hklinfo.hkl_of(h).format() << std::endl;

  // now read data columns
  std::vector<clipper::MTZcrystal> xtals;
  std::vector<clipper::MTZdataset> dsets;
  std::vector<clipper::HKL_data_base*> hkldata;
  for ( int i = 0; i < ipfiles.size(); i++ ) {
    clipper::MTZcrystal xtal;
    clipper::MTZdataset dset;
    clipper::HKL_data_base* data = NULL;
    clipper::CCP4MTZfile mtzin;
    mtzin.set_column_label_mode( clipper::CCP4MTZfile::Legacy );
    mtzin.open_read( ipfiles[i] );
    // get column names and types
    std::vector<clipper::String> cols = cols_all[i];
    clipper::String types = types_all[i];
    std::cout << ipcols[i] << " " << types << std::endl;
    if ( types.length() != cols.size() ) types = "";
    if      ( types == "FQ" )   data = new HKL_data<F_sigF> ( hklinfo );
    else if ( types == "JQ" )   data = new HKL_data<I_sigI> ( hklinfo );
    else if ( types == "GLGL" ) data = new HKL_data<F_sigF_anom>( hklinfo );
    else if ( types == "KMKM" ) data = new HKL_data<I_sigI_anom>( hklinfo );
    else if ( types == "FP" )   data = new HKL_data<F_phi>  ( hklinfo );
    else if ( types == "PW" )   data = new HKL_data<Phi_fom>( hklinfo );
    else if ( types == "AAAA" ) data = new HKL_data<ABCD>   ( hklinfo );
    else if ( types == "I" )    data = new HKL_data<Flag>   ( hklinfo );
    if ( data != NULL ) {
      mtzin.import_crystal( xtal, ipcols[i] );
      mtzin.import_dataset( dset, ipcols[i] );
      mtzin.import_hkl_data( *data, ipcols[i] );
      xtals.push_back( xtal );
      dsets.push_back( dset );
      hkldata.push_back( data );
    }
    mtzin.close_read();
  }


  // do column type conversions if required
  for ( int i = 0; i < hkldata.size(); i++ ) {
    int s = opcols[i].split(",").size();
    if ( hkldata[i]->type() == "Phi_fom" && s == 4 ) {
      std::cout << hkldata[i]->type() << " " << opcols[i] << std::endl;
      HKL_data<ABCD>* datanew = new HKL_data<ABCD> ( hklinfo );
      datanew->compute( *dynamic_cast<HKL_data<Phi_fom>*>(hkldata[i]),
                        clipper::data32::Compute_abcd_from_phifom() );
      hkldata[i] = datanew;
    }
    if ( hkldata[i]->type() == "ABCD" && s == 2 ) {
      std::cout << hkldata[i]->type() << " " << opcols[i] << std::endl;
      HKL_data<Phi_fom>* datanew = new HKL_data<Phi_fom> ( hklinfo );
      datanew->compute( *dynamic_cast<HKL_data<ABCD>*>(hkldata[i]),
                        clipper::data32::Compute_phifom_from_abcd() );
      hkldata[i] = datanew;
    }
    if ( hkldata[i]->type() == "I_sigI_anom" && s == 2 ) {
      std::cout << hkldata[i]->type() << " " << opcols[i] << std::endl;
      HKL_data<I_sigI>* datanew = new HKL_data<I_sigI> ( hklinfo );
      for ( clipper::HKL_info::HKL_reference_index ih = hklinfo.first(); !ih.last(); ih.next() ) {
        I_sigI_anom isigianom = (*dynamic_cast<HKL_data<I_sigI_anom>*>(hkldata[i]))[ih];
        clipper::datatypes::I_sigI<dtype> isigi;
        if        ( Util::is_nan( isigianom.I_pl() ) ) {
          isigi.I() = isigianom.I_mi();
          isigi.sigI() = isigianom.sigI_mi();
        } else if ( Util::is_nan( isigianom.I_mi() ) ) {
          isigi.I() = isigianom.I_pl();
          isigi.sigI() = isigianom.sigI_pl();
        } else {
          isigi.I() = 0.5 * ( isigianom.I_pl() + isigianom.I_mi() );
          isigi.sigI() = 0.5 * sqrt( isigianom.sigI_pl()*isigianom.sigI_pl() +
				                             isigianom.sigI_mi()*isigianom.sigI_mi() );
        }
        (*datanew)[ih] = isigi;
      }
      hkldata[i] = datanew;
    }
    if ( hkldata[i]->type() == "F_sigF_anom" && s == 2 ) {
      std::cout << hkldata[i]->type() << " " << opcols[i] << std::endl;
      HKL_data<F_sigF>* datanew = new HKL_data<F_sigF> ( hklinfo );
      for ( clipper::HKL_info::HKL_reference_index ih = hklinfo.first(); !ih.last(); ih.next() ) {
        F_sigF_anom fsigfanom = (*dynamic_cast<HKL_data<F_sigF_anom>*>(hkldata[i]))[ih];
        clipper::datatypes::F_sigF<dtype> fsigf;
        if        ( Util::is_nan( fsigfanom.f_pl() ) ) {
          fsigf.f() = fsigfanom.f_mi();
          fsigf.sigf() = fsigfanom.sigf_mi();
        } else if ( Util::is_nan( fsigfanom.f_mi() ) ) {
          fsigf.f() = fsigfanom.f_pl();
          fsigf.sigf() = fsigfanom.sigf_pl();
        } else {
          fsigf.f() = 0.5 * ( fsigfanom.f_pl() + fsigfanom.f_mi() );
          fsigf.sigf() = 0.5 * sqrt( fsigfanom.sigf_pl()*fsigfanom.sigf_pl() +
				                             fsigfanom.sigf_mi()*fsigfanom.sigf_mi() );
        }
        (*datanew)[ih] = fsigf;
      }
      hkldata[i] = datanew;
    }
  }


  // now write out the lot
  clipper::CCP4MTZfile mtzout;
  mtzout.open_write( opfile );
  mtzout.export_hkl_info( hklinfo );
  for ( int i = 0; i < hkldata.size(); i++ ) {
    clipper::String path = ( "/" + xtals[i].crystal_name() + 
			     "/" + dsets[i].dataset_name() +
			     "/["+opcols[i]+"]" );
    std::cout << path << std::endl;
    mtzout.export_crystal( xtals[i], path );
    mtzout.export_dataset( dsets[i], path );
    mtzout.export_hkl_data( *(hkldata[i]), path );
  }
  mtzout.close_write();
}
