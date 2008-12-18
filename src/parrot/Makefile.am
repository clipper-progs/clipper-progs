bin_PROGRAMS = cparrot

cparrot_SOURCES = \
cparrot.cpp         parrot-ncsaver.h    parrot-ncs.h      simulate-lib.h \
parrot-lib.cpp      parrot-ncs.cpp      simplex-lib.cpp \
parrot-lib.h        parrot-ncsfind.cpp  simplex-lib.h \
parrot-ncsaver.cpp  parrot-ncsfind.h    simulate-lib.cpp


INCLUDES = $(CLIPPER_CXXFLAGS) $(FFTW_CXXFLAGS) $(MMDB_CXXFLAGS) \
           $(CCP4_CXXFLAGS) -I$(includedir)

LIBS = $(CLIPPER_LIBS) $(FFTW_LIBS) $(CCP4_LIBS) $(MMDB_LIBS) @LIBS@


