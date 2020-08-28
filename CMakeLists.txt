cmake_minimum_required(VERSION 2.8)
project(clipper-progs C CXX)

find_package(CCP4 REQUIRED mmdb2 ccp4c clipper-core clipper-ccp4 clipper-mmdb
                           clipper-contrib clipper-minimol clipper-cif)
include_directories(${CCP4_INCLUDE_DIRS})

# clipper-cif is needed only for cmakereference, it will be added explicitely
list(REMOVE_ITEM CCP4_LIBRARIES ${CLIPPER-CIF_LIBRARY})

if (WIN32)
  set(SOCKET_LIB ws2_32)
endif()

add_executable(cextend src/cextend.cpp)
target_link_libraries(cextend ${CCP4_LIBRARIES})

add_executable(ccontacts src/ccontacts.cpp)
target_link_libraries(ccontacts ${CCP4_LIBRARIES})

add_executable(csheetbend src/sheetbend/csheetbend.cpp
	                  src/sheetbend/shiftfield.cpp
	                  src/sheetbend/pseudoregularize.cpp )
target_link_libraries(csheetbend ${CCP4_LIBRARIES})

add_executable(chltofom src/chltofom.cpp)
target_link_libraries(chltofom ${CCP4_LIBRARIES})

add_executable(cphasematch src/cphasematch.cpp)
target_link_libraries(cphasematch ${CCP4_LIBRARIES})

add_executable(csfcalc src/csfcalc.cpp)
target_link_libraries(csfcalc ${CCP4_LIBRARIES})

add_executable(cfft src/cfft.cpp)
target_link_libraries(cfft ${CCP4_LIBRARIES})

add_executable(cinvfft src/cinvfft.cpp)
target_link_libraries(cinvfft ${CCP4_LIBRARIES})

add_executable(csigmaa src/csigmaa.cpp)
target_link_libraries(csigmaa ${CCP4_LIBRARIES})

add_executable(cecalc src/cecalc.cpp)
target_link_libraries(cecalc ${CCP4_LIBRARIES})

add_executable(cphasecombine src/cphasecombine.cpp)
target_link_libraries(cphasecombine ${CCP4_LIBRARIES})

add_executable(cpatterson src/cpatterson.cpp)
target_link_libraries(cpatterson ${CCP4_LIBRARIES})

add_executable(caniso src/caniso.cpp)
target_link_libraries(caniso ${CCP4_LIBRARIES})

add_executable(convert2mtz src/convert2mtz.cpp)
target_link_libraries(convert2mtz ${CCP4_LIBRARIES})

add_executable(csymmatch src/csymmatch.cpp)
target_link_libraries(csymmatch ${CCP4_LIBRARIES})

add_executable(cmaplocal src/cmaplocal.cpp)
target_link_libraries(cmaplocal ${CCP4_LIBRARIES})

add_executable(cncsfrommodel src/cncsfrommodel.cpp)
target_link_libraries(cncsfrommodel ${CCP4_LIBRARIES})

add_executable(cmodeltoseq src/cmodeltoseq.cpp)
target_link_libraries(cmodeltoseq ${CCP4_LIBRARIES})

add_executable(cmapcut src/cmapcut.cpp)
target_link_libraries(cmapcut ${CCP4_LIBRARIES})

add_executable(comit src/comit.cpp)
target_link_libraries(comit ${CCP4_LIBRARIES})

add_executable(cmapcoeff src/cmapcoeff.cpp)
target_link_libraries(cmapcoeff ${CCP4_LIBRARIES})

add_executable(cmtzjoin src/cmtzjoin.cpp)
target_link_libraries(cmtzjoin ${CCP4_LIBRARIES})

add_executable(cmtzsplit src/cmtzsplit.cpp)
target_link_libraries(cmtzsplit ${CCP4_LIBRARIES})

add_executable(ctruncate
 src/ctruncate/ctruncate.cpp
 src/ctruncate/ctruncate_utils.cpp
 src/ctruncate/cpsf_utils.cpp
 src/ctruncate/twinlaws.cpp
 src/ctruncate/ctruncate_truncate.cpp
 src/ctruncate/ctruncate_twin.cpp
 src/ctruncate/ctruncate_parity.cpp
 src/ctruncate/ctruncate_moments.cpp
 src/ctruncate/ctruncate_analyse.cpp
 src/ctruncate/ctruncate_matthews.cpp
 src/ctruncate/ctruncate_wilson.cpp
 src/ctruncate/intensity_scale.cpp
 src/ctruncate/alt_hkl_datatypes.cpp
 src/ctruncate/ctruncate_aniso.cpp
 src/ctruncate/best.cpp
 src/ctruncate/best_rna.cpp)
target_link_libraries(ctruncate ${CCP4_LIBRARIES})

add_executable(cpirate
 src/pirate/cpirate.cpp
 src/pirate/simulate-lib.cpp
 src/pirate/pirate-lib.cpp
 src/pirate/pirate-ncsaver.cpp
 src/pirate/pirate-ncsfind.cpp
 src/pirate/simplex-lib.cpp)
target_link_libraries(cpirate ${CCP4_LIBRARIES})

add_executable(cmakereference
 src/pirate/cmakereference.cpp
 src/pirate/ftplib.c
 src/pirate/compress42.c)
target_link_libraries(cmakereference
                      ${CLIPPER-CIF_LIBRARY} ${CCP4_LIBRARIES} ${SOCKET_LIB})

add_executable(cbuccaneer
 src/buccaneer/buccaneer-build.cpp
 src/buccaneer/buccaneer-prune.cpp
 src/buccaneer/buccaneer-lib.cpp
 src/buccaneer/buccaneer-correct.cpp
 src/buccaneer/buccaneer-sequence.cpp
 src/buccaneer/buccaneer-link.cpp
 src/buccaneer/buccaneer-filter.cpp
 src/buccaneer/buccaneer-util.cpp
 src/buccaneer/buccaneer-ncsbuild.cpp
 src/buccaneer/buccaneer-find.cpp
 src/buccaneer/cbuccaneer.cpp
 src/buccaneer/buccaneer-prep.cpp
 src/buccaneer/simplex-lib.cpp
 src/buccaneer/buccaneer-grow.cpp
 src/buccaneer/buccaneer-prot.cpp
 src/buccaneer/simulate-lib.cpp
 src/buccaneer/buccaneer-join.cpp
 src/buccaneer/buccaneer-known.cpp
 src/buccaneer/buccaneer-merge.cpp
 src/buccaneer/buccaneer-tidy.cpp)
target_link_libraries(cbuccaneer ${CCP4_LIBRARIES})

add_executable(csequins
 src/buccaneer/csequins.cpp
 src/buccaneer/buccaneer-sequence.cpp
 src/buccaneer/simulate-lib.cpp
 src/buccaneer/buccaneer-prot.cpp
 src/buccaneer/buccaneer-lib.cpp
 src/buccaneer/buccaneer-tidy.cpp)
target_link_libraries(csequins ${CCP4_LIBRARIES})

add_executable(cparrot
 src/parrot/cparrot.cpp
 src/parrot/parrot-lib.cpp
 src/parrot/parrot-ncs.cpp
 src/parrot/simplex-lib.cpp
 src/parrot/parrot-ncsfind.cpp
 src/parrot/parrot-ncsaver.cpp
 src/parrot/simulate-lib.cpp)
target_link_libraries(cparrot ${CCP4_LIBRARIES})

add_executable(csloop
 src/sloop/csloop.cpp
 src/sloop/protein_db.cpp
 src/sloop/protein_db_utils.cpp)
target_link_libraries(csloop ${CCP4_LIBRARIES})

add_executable(cnautilus
 src/nautilus/cnautilus.cpp
 src/nautilus/nautilus-sequence.cpp
 src/nautilus/nautilus-tools.cpp
 src/nautilus/nautilus-join.cpp
 src/nautilus/nautilus-ss-find.cpp
 src/nautilus/nautilus-util.cpp
 src/nautilus/nautilus-tidy.cpp
 src/nautilus/nautilus-rebuild-bases.cpp
 src/nautilus/nautilus-target.cpp
 src/nautilus/nucleicacid_db.cpp)
target_link_libraries(cnautilus ${CCP4_LIBRARIES})


install(FILES src/nautilus/data/nautilus_lib.pdb src/sloop/prot500.db
        DESTINATION lib/data)

install(TARGETS ccontacts cextend chltofom cphasematch csfcalc cfft cinvfft csigmaa
                csheetbend cecalc cphasecombine cpatterson caniso convert2mtz
                csymmatch cmaplocal cncsfrommodel cmodeltoseq cmapcut
                comit cmapcoeff cmtzsplit cmtzjoin
                ctruncate cpirate cmakereference
                cbuccaneer csequins
                cparrot csloop cnautilus
        DESTINATION bin)

