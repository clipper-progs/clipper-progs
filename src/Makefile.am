bin_PROGRAMS = cpirate chltofom

cpirate_SOURCES = cpirate.cpp ccp4-extras.cpp piralib.cpp ccp4-extras.h piralib.h
chltofom_SOURCES = chltofom.cpp ccp4-extras.cpp

noinst_SCRIPTS = makereference.csh

EXTRA_DIST = README

INCLUDES = $(FFTW_CXXFLAGS) $(MMDB_CXXFLAGS)  $(CCP4_CXXFLAGS) $(CCTBX_CXXFLAGS) \
	   $(BOOST_CXXFLAGS) -I$(includedir)

LDFLAGS = $(FFTW_LIBS) $(CCP4_LIBS) $(MMDB_LIBS) $(CCTBX_LIBS) 

LDADD = LDADD = ../clipper/core/libclipper-core.la ../clipper/contrib/libclipper-contrib.la ../clipper/cif/libclipper-cif.la ../clipper/mmdb/libclipper-mmdb.la ../clipper/ccp4/libclipper-ccp4.la ../clipper/phs/libclipper-phs.la

