bin_PROGRAMS = cpirate chltofom cmakereference cphasematch csfcalc cfft cinvfft csigmaa cecalc cphasecombine cpatterson

cpirate_SOURCES = cpirate.cpp piralib.cpp piralib.h pirancslib.cpp pirancslib.h

chltofom_SOURCES = chltofom.cpp

cmakereference_SOURCES = cmakereference.cpp ftplib.h ftplib.c compress42.c

cphasematch_SOURCES = cphasematch.cpp

csfcalc_SOURCES = csfcalc.cpp

cfft_SOURCES = cfft.cpp

cinvfft_SOURCES = cinvfft.cpp

csigmaa_SOURCES = csigmaa.cpp

cecalc_SOURCES = cecalc.cpp

cphasecombine_SOURCES = cphasecombine.cpp

cpatterson_SOURCES = cpatterson.cpp

EXTRA_DIST = README ftplib.c compress42.c

INCLUDES = $(FFTW_CXXFLAGS) $(MMDB_CXXFLAGS)  $(CCP4_CXXFLAGS) $(CCTBX_CXXFLAGS) \
	   $(BOOST_CXXFLAGS) -I$(includedir)

LIBS = $(FFTW_LIBS) $(CCP4_LIBS) $(MMDB_LIBS) $(CCTBX_LIBS) @LIBS@

LDADD = ../clipper/contrib/libclipper-contrib.la ../clipper/cif/libclipper-cif.la ../clipper/mmdb/libclipper-mmdb.la ../clipper/ccp4/libclipper-ccp4.la ../clipper/phs/libclipper-phs.la ../clipper/core/libclipper-core.la

