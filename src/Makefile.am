DIST_SUBDIRS = pirate buccaneer
SUBDIRS = pirate buccaneer

bin_PROGRAMS = chltofom cphasematch csfcalc cfft cinvfft csigmaa cecalc cphasecombine cpatterson caniso convert2mtz

chltofom_SOURCES = chltofom.cpp

cphasematch_SOURCES = cphasematch.cpp

csfcalc_SOURCES = csfcalc.cpp

cfft_SOURCES = cfft.cpp

cinvfft_SOURCES = cinvfft.cpp

csigmaa_SOURCES = csigmaa.cpp

cecalc_SOURCES = cecalc.cpp

cphasecombine_SOURCES = cphasecombine.cpp

cpatterson_SOURCES = cpatterson.cpp

caniso_SOURCES = caniso.cpp

convert2mtz_SOURCES = convert2mtz.cpp

EXTRA_DIST = README

INCLUDES = $(FFTW_CXXFLAGS) $(MMDB_CXXFLAGS)  $(CCP4_CXXFLAGS) $(CCTBX_CXXFLAGS) \
	   $(BOOST_CXXFLAGS) -I$(includedir)

LIBS = $(FFTW_LIBS) $(CCP4_LIBS) $(MMDB_LIBS) $(CCTBX_LIBS) @LIBS@

LDADD = ../clipper/contrib/libclipper-contrib.la ../clipper/cif/libclipper-cif.la ../clipper/mmdb/libclipper-mmdb.la ../clipper/ccp4/libclipper-ccp4.la ../clipper/phs/libclipper-phs.la ../clipper/core/libclipper-core.la

