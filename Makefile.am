ACLOCAL_AMFLAGS = -I m4

dist_html_DATA = html/cpirate.html html/chltofom.html html/cbuccaneer.html \
		 html/csequins.html html/parrot.html html/cmapcut.html \
		 html/ctruncate.html html/comit.html html/cnautilus.html

bin_PROGRAMS = src/chltofom src/cphasematch src/csfcalc src/cfft \
	       src/cinvfft src/csigmaa src/cecalc src/cphasecombine \
	       src/cpatterson src/caniso src/convert2mtz src/ctruncate \
	       src/csymmatch src/cmaplocal src/cncsfrommodel \
	       src/cmodeltoseq src/cmapcut src/comit \
	       src/pirate/cpirate src/pirate/cmakereference \
	       src/buccaneer/cbuccaneer src/buccaneer/csequins \
	       src/parrot/cparrot \
	       src/sloop/csloop \
	       src/nautilus/cnautilus

AM_CXXFLAGS = $(FFTW_CXXFLAGS) $(PTHREAD_CFLAGS)

LDADD = -lclipper-ccp4 -lclipper-contrib -lclipper-minimol -lclipper-mmdb \
	-lclipper-core $(PTHREAD_LIBS) -lccp4c -lmmdb $(FFTW_LIBS)

src_pirate_cmakereference_LDADD = -lclipper-cif $(LDADD)


src_chltofom_SOURCES = src/chltofom.cpp

src_cphasematch_SOURCES = src/cphasematch.cpp

src_csfcalc_SOURCES = src/csfcalc.cpp

src_cfft_SOURCES = src/cfft.cpp

src_cinvfft_SOURCES = src/cinvfft.cpp

src_csigmaa_SOURCES = src/csigmaa.cpp

src_cecalc_SOURCES = src/cecalc.cpp

src_cphasecombine_SOURCES = src/cphasecombine.cpp

src_cpatterson_SOURCES = src/cpatterson.cpp

src_caniso_SOURCES = src/caniso.cpp

src_convert2mtz_SOURCES = src/convert2mtz.cpp src/convert2mtz.h

src_ctruncate_SOURCES = src/ctruncate.cpp \
			src/ctruncate_utils.cpp \
			src/cpsf_utils.cpp \
			src/twinlaws.cpp \
			src/ctruncate_truncate.cpp \
			src/ctruncate_twin.cpp \
			src/ctruncate_parity.cpp \
			src/ctruncate_moments.cpp \
			src/ctruncate_analyse.cpp \
			src/ctruncate_matthews.cpp \
			src/ctruncate_wilson.cpp \
			src/intensity_scale.cpp \
			src/best.cpp \
			src/ctruncate_utils.h \
			src/cpsf_utils.h \
			src/twinlaws.h \
			src/ctruncate_truncate.h \
			src/ctruncate_twin.h \
			src/ctruncate_parity.h \
			src/ctruncate_moments.h \
			src/ctruncate_analyse.h \
			src/ctruncate_matthews.h \
			src/ctruncate_wilson.h \
			src/intensity_scale.h \
			src/best.h \
			src/intensity_target.h \
			src/alt_hkl_datatypes.h

src_csymmatch_SOURCES = src/csymmatch.cpp

src_cmaplocal_SOURCES = src/cmaplocal.cpp

src_cncsfrommodel_SOURCES = src/cncsfrommodel.cpp

src_cmodeltoseq_SOURCES = src/cmodeltoseq.cpp

src_cmapcut_SOURCES = src/cmapcut.cpp

src_comit_SOURCES = src/comit.cpp


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



sloopdatadir = $(datadir)/sloop
dist_sloopdata_DATA = src/sloop/prot500.db


src_nautilus_cnautilus_SOURCES = src/nautilus/cnautilus.cpp \
				 src/nautilus/nautilus-sequence.cpp \
				 src/nautilus/nautilus-tools.cpp \
				 src/nautilus/nautilus-join.cpp \
				 src/nautilus/nautilus-ss-find.cpp \
				 src/nautilus/nautilus-util.cpp \
				 src/nautilus/nautilus-rebuild-bases.cpp \
				 src/nautilus/nautilus-target.cpp \
				 src/nautilus/nucleicacid_db.cpp \
				 src/nautilus/nautilus-join.h \
				 src/nautilus/nautilus-ss-find.h \
				 src/nautilus/nautilus-util.h \
				 src/nautilus/nautilus-rebuild-bases.h \
				 src/nautilus/nautilus-target.h \
				 src/nautilus/nucleicacid_db.h \
				 src/nautilus/nautilus-sequence.h \
				 src/nautilus/nautilus-tools.h

