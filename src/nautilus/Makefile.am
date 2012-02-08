bin_PROGRAMS = cnautilus

cnautilus_SOURCES = \
cnautilus.cpp			nautilus-sequence.cpp		nautilus-tools.cpp \
nautilus-join.cpp		nautilus-ss-find.cpp		nautilus-util.cpp \
nautilus-rebuild-bases.cpp	nautilus-target.cpp		nucleicacid_db.cpp

noinst_HEADERS = \
nautilus-join.h			nautilus-ss-find.h		nautilus-util.h \
nautilus-rebuild-bases.h	nautilus-target.h		nucleicacid_db.h \
nautilus-sequence.h		nautilus-tools.h

AM_CPPFLAGS = $(CLIPPER_CXXFLAGS) $(FFTW_CXXFLAGS) $(MMDB_CXXFLAGS) \
           $(CCP4_CXXFLAGS) -I$(includedir)

LIBS = $(PTHREAD_CFLAGS) $(CLIPPER_LIBS) $(FFTW_LIBS) $(CCP4_LIBS) $(MMDB_LIBS) @LIBS@ $(PTHREAD_LIBS)


