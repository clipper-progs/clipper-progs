DIST_SUBDIRS = pirate buccaneer parrot sloop
SUBDIRS = pirate buccaneer parrot sloop

bin_PROGRAMS = chltofom cphasematch csfcalc cfft cinvfft csigmaa cecalc cphasecombine cpatterson caniso convert2mtz ctruncate csymmatch cmaplocal cncsfrommodel cmodeltoseq cmapcut

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

ctruncate_SOURCES = ctruncate.cpp ctruncate_utils.cpp cpsf_utils.cpp twinlaws.cpp ctruncate_truncate.cpp ctruncate_twin.cpp ctruncate_parity.cpp ctruncate_moments.cpp ctruncate_analyse.cpp ctruncate_matthews.cpp ctruncate_wilson.cpp intensity_scale.cpp best.cpp

csymmatch_SOURCES = csymmatch.cpp

cmaplocal_SOURCES = cmaplocal.cpp

cncsfrommodel_SOURCES = cncsfrommodel.cpp

cmodeltoseq_SOURCES = cmodeltoseq.cpp

cmapcut_SOURCES = cmapcut.cpp

EXTRA_DIST = README

INCLUDES = $(PTHREAD_CFLAGS) $(CLIPPER_CXXFLAGS) $(FFTW_CXXFLAGS) $(MMDB_CXXFLAGS)  $(CCP4_CXXFLAGS) $(CCTBX_CXXFLAGS) \
	   $(BOOST_CXXFLAGS) -I$(includedir)

LIBS = $(PTHREAD_CFLAGS) $(CLIPPER_LIBS) $(FFTW_LIBS) $(CCP4_LIBS) $(MMDB_LIBS) $(CCTBX_LIBS) $(PTHREAD_LIBS) @LIBS@

