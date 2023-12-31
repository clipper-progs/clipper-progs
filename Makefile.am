ACLOCAL_AMFLAGS = -I m4

dist_data_DATA = src/nautilus/data/nautilus_lib.pdb src/sloop/prot500.db

# duplicated (also in documentation/html), not installed
dist_noinst_DATA = html/cpirate.html html/chltofom.html html/cbuccaneer.html \
		 html/csequins.html html/parrot.html html/cmapcut.html \
		 html/ctruncate.html html/comit.html html/cnautilus.html

bin_PROGRAMS = src/chltofom src/cphasematch src/csfcalc src/cfft \
	       src/cinvfft src/csigmaa src/cecalc src/cphasecombine \
	       src/cpatterson src/caniso src/convert2mtz \
	       src/csymmatch src/cmaplocal src/cncsfrommodel \
	       src/cmodeltoseq src/cmapcut src/comit \
	       src/cmapcoeff \
	       src/cmtzsplit src/cmtzjoin \
	       src/clipper_test \
	       src/ctruncate/ctruncate \
	       src/pirate/cpirate src/pirate/cmakereference \
	       src/buccaneer/cbuccaneer src/buccaneer/csequins \
	       src/parrot/cparrot \
	       src/sloop/csloop \
	       src/nautilus/cnautilus

AM_CXXFLAGS = $(FFTW_CXXFLAGS) $(PTHREAD_CFLAGS)

# LDADD is used for programs for which program_LDADD is not set
LDADD = -lclipper-ccp4 -lclipper-contrib -lclipper-minimol -lclipper-mmdb \
	-lclipper-core $(PTHREAD_LIBS) -lccp4c -lmmdb2 $(FFTW_LIBS)

src_pirate_cmakereference_LDADD = -lclipper-cif $(LDADD) $(SOCKET_LIB)


src_chltofom_SOURCES = src/chltofom.cpp

src_cphasematch_SOURCES = src/cphasematch.cpp

src_csfcalc_SOURCES = src/csfcalc.cpp

src_cfft_SOURCES = src/cfft.cpp

src_clipper_test_SOURCES = src/clipper_test.cpp

src_cinvfft_SOURCES = src/cinvfft.cpp

src_csigmaa_SOURCES = src/csigmaa.cpp

src_cecalc_SOURCES = src/cecalc.cpp

src_cphasecombine_SOURCES = src/cphasecombine.cpp

src_cpatterson_SOURCES = src/cpatterson.cpp

src_caniso_SOURCES = src/caniso.cpp

src_convert2mtz_SOURCES = src/convert2mtz.cpp src/convert2mtz.h

src_csymmatch_SOURCES = src/csymmatch.cpp

src_cmaplocal_SOURCES = src/cmaplocal.cpp

src_cncsfrommodel_SOURCES = src/cncsfrommodel.cpp

src_cmodeltoseq_SOURCES = src/cmodeltoseq.cpp

src_cmapcut_SOURCES = src/cmapcut.cpp

src_comit_SOURCES = src/comit.cpp

src_cmapcoeff_SOURCES = src/cmapcoeff.cpp

src_cmtzjoin_SOURCES = src/cmtzjoin.cpp

src_cmtzsplit_SOURCES = src/cmtzsplit.cpp

src_ctruncate_ctruncate_SOURCES = src/ctruncate/ctruncate.cpp \
				  src/ctruncate/ctruncate_utils.cpp \
				  src/ctruncate/cpsf_utils.cpp \
				  src/ctruncate/twinlaws.cpp \
				  src/ctruncate/ctruncate_truncate.cpp \
				  src/ctruncate/ctruncate_twin.cpp \
				  src/ctruncate/ctruncate_parity.cpp \
				  src/ctruncate/ctruncate_moments.cpp \
				  src/ctruncate/ctruncate_analyse.cpp \
				  src/ctruncate/ctruncate_matthews.cpp \
				  src/ctruncate/ctruncate_wilson.cpp \
				  src/ctruncate/intensity_scale.cpp \
				  src/ctruncate/best.cpp \
                                  src/ctruncate/best_rna.cpp \
                                  src/ctruncate/alt_hkl_datatypes.cpp \
                                  src/ctruncate/ctruncate_aniso.cpp \
				  src/ctruncate/ctruncate_utils.h \
				  src/ctruncate/cpsf_utils.h \
				  src/ctruncate/twinlaws.h \
				  src/ctruncate/ctruncate_truncate.h \
				  src/ctruncate/ctruncate_twin.h \
				  src/ctruncate/ctruncate_parity.h \
				  src/ctruncate/ctruncate_moments.h \
				  src/ctruncate/ctruncate_analyse.h \
				  src/ctruncate/ctruncate_matthews.h \
				  src/ctruncate/ctruncate_wilson.h \
				  src/ctruncate/intensity_scale.h \
				  src/ctruncate/best.h \
                                  src/ctruncate/best_rna.h \
				  src/ctruncate/intensity_target.h \
				  src/ctruncate/alt_hkl_datatypes.h \
                                  src/ctruncate/ctruncate_aniso.h

src_pirate_cpirate_SOURCES = src/pirate/cpirate.cpp \
			     src/pirate/simulate-lib.cpp \
			     src/pirate/pirate-lib.cpp \
			     src/pirate/pirate-ncsaver.cpp \
			     src/pirate/pirate-ncsfind.cpp \
			     src/pirate/simplex-lib.cpp \
			     src/pirate/simulate-lib.h \
			     src/pirate/pirate-lib.h \
			     src/pirate/pirate-ncsaver.h \
			     src/pirate/pirate-ncsfind.h \
			     src/pirate/simplex-lib.h

src_pirate_cmakereference_SOURCES = src/pirate/cmakereference.cpp \
				    src/pirate/ftplib.h \
				    src/pirate/ftplib.c \
				    src/pirate/compress42.c

EXTRA_DIST = src/pirate/makereference.csh


src_buccaneer_cbuccaneer_SOURCES = src/buccaneer/buccaneer-build.cpp \
				   src/buccaneer/buccaneer-join.h \
				   src/buccaneer/buccaneer-prune.cpp \
				   src/buccaneer/buccaneer-build.h \
				   src/buccaneer/buccaneer-lib.cpp \
				   src/buccaneer/buccaneer-prune.h \
				   src/buccaneer/buccaneer-correct.cpp \
				   src/buccaneer/buccaneer-lib.h \
				   src/buccaneer/buccaneer-sequence.cpp \
				   src/buccaneer/buccaneer-correct.h \
				   src/buccaneer/buccaneer-link.cpp \
				   src/buccaneer/buccaneer-sequence.h \
				   src/buccaneer/buccaneer-filter.cpp \
				   src/buccaneer/buccaneer-link.h \
				   src/buccaneer/buccaneer-util.cpp \
				   src/buccaneer/buccaneer-filter.h \
				   src/buccaneer/buccaneer-ncsbuild.cpp \
				   src/buccaneer/buccaneer-util.h \
				   src/buccaneer/buccaneer-find.cpp \
				   src/buccaneer/buccaneer-ncsbuild.h \
				   src/buccaneer/cbuccaneer.cpp \
				   src/buccaneer/buccaneer-find.h \
				   src/buccaneer/buccaneer-prep.cpp \
				   src/buccaneer/simplex-lib.cpp \
				   src/buccaneer/buccaneer-grow.cpp \
				   src/buccaneer/buccaneer-prep.h \
				   src/buccaneer/simplex-lib.h \
				   src/buccaneer/buccaneer-grow.h \
				   src/buccaneer/buccaneer-prot.cpp \
				   src/buccaneer/simulate-lib.cpp \
				   src/buccaneer/buccaneer-join.cpp \
				   src/buccaneer/buccaneer-prot.h \
				   src/buccaneer/simulate-lib.h \
				   src/buccaneer/buccaneer-known.cpp \
				   src/buccaneer/buccaneer-known.h \
				   src/buccaneer/buccaneer-merge.cpp \
				   src/buccaneer/buccaneer-merge.h \
				   src/buccaneer/buccaneer-tidy.cpp \
				   src/buccaneer/buccaneer-tidy.h

src_buccaneer_csequins_SOURCES = src/buccaneer/csequins.cpp \
				 src/buccaneer/buccaneer-sequence.cpp \
				 src/buccaneer/simulate-lib.cpp \
				 src/buccaneer/buccaneer-prot.cpp \
				 src/buccaneer/buccaneer-lib.cpp


src_parrot_cparrot_SOURCES = src/parrot/cparrot.cpp \
			     src/parrot/parrot-ncsaver.h \
			     src/parrot/parrot-ncs.h \
			     src/parrot/simulate-lib.h \
			     src/parrot/parrot-lib.cpp \
			     src/parrot/parrot-ncs.cpp \
			     src/parrot/simplex-lib.cpp \
			     src/parrot/parrot-lib.h \
			     src/parrot/parrot-ncsfind.cpp \
			     src/parrot/simplex-lib.h \
			     src/parrot/parrot-ncsaver.cpp \
			     src/parrot/parrot-ncsfind.h \
			     src/parrot/simulate-lib.cpp


src_sloop_csloop_SOURCES = src/sloop/csloop.cpp \
			   src/sloop/protein_db.cpp \
			   src/sloop/protein_db_utils.cpp \
			   src/sloop/protein_db.h \
			   src/sloop/protein_db_utils.h


src_nautilus_cnautilus_SOURCES = src/nautilus/cnautilus.cpp \
				 src/nautilus/nautilus-sequence.cpp \
				 src/nautilus/nautilus-tools.cpp \
				 src/nautilus/nautilus-join.cpp \
				 src/nautilus/nautilus-ss-find.cpp \
				 src/nautilus/nautilus-util.cpp \
				 src/nautilus/nautilus-tidy.cpp \
				 src/nautilus/nautilus-rebuild-bases.cpp \
				 src/nautilus/nautilus-target.cpp \
				 src/nautilus/nucleicacid_db.cpp \
				 src/nautilus/nautilus-join.h \
				 src/nautilus/nautilus-ss-find.h \
				 src/nautilus/nautilus-util.h \
				 src/nautilus/nautilus-tidy.h \
				 src/nautilus/nautilus-rebuild-bases.h \
				 src/nautilus/nautilus-target.h \
				 src/nautilus/nucleicacid_db.h \
				 src/nautilus/nautilus-sequence.h \
				 src/nautilus/nautilus-tools.h

