bin_PROGRAMS = csloop

csloop_SOURCES = protein_db.cpp        protein_db_utils.h
protein_db.h     csloop.cpp            protein_db_utils.cpp


INCLUDES = $(PTHREAD_CFLAGS) $(CLIPPER_CXXFLAGS) $(FFTW_CXXFLAGS) $(MMDB_CXXFLAGS)  $(CCP4_CXXFLAGS) \
	   -I$(includedir)

LIBS = $(PTHREAD_CFLAGS) $(CLIPPER_LIBS) $(FFTW_LIBS) $(CCP4_LIBS) $(MMDB_LIBS) $(CCTBX_LIBS) $(PTHREAD_LIBS) @LIBS@

pkgdatadir = $(datadir)/sloop

dist_pkgdata_DATA = prot500.db

