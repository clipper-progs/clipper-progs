bin_PROGRAMS = cpirate chltofom

cpirate_SOURCES = cpirate.cpp ccp4-extras.cpp piralib.cpp ccp4-extras.h piralib.h
chltofom_SOURCES = chltofom.cpp ccp4-extras.cpp

noinst_SCRIPTS = makereference.csh

EXTRA_DIST = README

INCLUDES = $(FFTW_CXXFLAGS) $(MMDB_CXXFLAGS)  $(CCP4_CXXFLAGS) $(CCTBX_CXXFLAGS) \
	   $(BOOST_CXXFLAGS) -I$(includedir)

LDFLAGS = -L$(libdir) -lclipper-contrib -lclipper-cif -lclipper-mmdb  -lclipper-ccp4 -lclipper-phs -lclipper-core $(FFTW_LIBS) $(CCP4_LIBS) $(MMDB_LIBS) $(CCTBX_LIBS) 
