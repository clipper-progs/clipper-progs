
# Process this file with autoconf to produce a configure script.

dnl package-name version-number bugs-to:
AC_INIT(clipper_progs, 0.1.0, ccb@dl.ac.uk) 

AC_CANONICAL_TARGET

dnl package-name version-number, the gzip'd tar file is named
dnl package-version.tar.gz.  
AM_INIT_AUTOMAKE(clipper_progs, 0.1.0)

AC_CONFIG_SRCDIR([src/buccaneer])

AM_MAINTAINER_MODE

AC_LANG(C++)

# Checks for programs.
AC_PROG_CXX(xlc++ CC cxx c++ icc g++)
AC_PROG_CC(xlc cc icc gcc)
AC_PROG_F77
AC_HEADER_STDC
AC_PROG_INSTALL
#AC_PROG_RANLIB

AX_PTHREAD(,
    [ AC_MSG_ERROR([Failed to find threading library])
    ] )

AC_ARG_ENABLE(debug,
  AC_HELP_STRING( [--enable-debug], [ enable debugging options ] ),
  [ac_cv_enable_debug=yes ],
  [ac_cv_enable_debug=no ] )

AC_PROGS_OPTIONS

dnl Check for FFTW.
AM_PATH_FFTW(, AC_MSG_ERROR([locating necessary fftw]) )

dnl Optional: Check for MMDB.
AM_PATH_MMDB(,
  [ if test "x$mmdb_prefix" != x ; then
      if test "x$ac_MMDB_CXXFLAGS" != x ; then
        case "$ac_MMDB_CXXFLAGS" in
        *src/mmdb )
        ac_MMDB_LDOPTS=`echo "$ac_MMDB_CXXFLAGS" | sed s#src/mmdb#src#g`
        ;;
        *include/mmdb )
        ac_MMDB_LDOPTS=`echo "$ac_MMDB_CXXFLAGS" | sed s#include/mmdb#lib#g`
        ;;
        *)
        esac
        MMDB_CXXFLAGS="$ac_MMDB_CXXFLAGS"
        MMDB_LIBS=`echo "$ac_MMDB_LDOPTS -lmmdb" | sed s#-I#-L#`
        AC_MSG_WARN([assuming library $MMDB_LIBS headers $MMDB_CXXFLAGS])
        have_mmdb=yes
      else
        MMDB_CXXFLAGS="-I$mmdb_prefix/src -I$mmdb_prefix/include"
        MMDB_LIBS="-L$mmdb_prefix/src -L$mmdb_prefix/lib -lmmdb"
        AC_MSG_WARN([assuming library $MMDB_LIBS headers $MMDB_CXXFLAGS])
        have_mmdb=yes
      fi
    else
      AC_MSG_ERROR([Failed to find mmdb])
    fi
  ]
)

AM_PATH_CCP4(,
  [ if test "x$ccp4_prefix" != x ; then
     if test "x$ac_CCP4_CXXFLAGS" != x ; then 
       case "$ac_CCP4_CXXFLAGS" in
        *include/ccp4 )
        ac_CCP4_LDOPTS=`echo "$ac_CCP4_CXXFLAGS" | sed s#include/ccp4#lib#g`
        ;;
        *)
        ac_CCP4_LDOPTS="$ac_CCP4_CXXFLAGS"
        esac
        ac_CCP4_LDOPTS=`echo "$ac_CCP4_LDOPTS" | sed s#-I#-L#`
        CCP4_CXXFLAGS="$ac_CCP4_CXXFLAGS" 
        CCP4_LIBS="$ac_CCP4_LDOPTS -lccp4c"
        AC_MSG_WARN([assuming library $CCP4_LIBS headers $CCP4_CXXFLAGS])
        have_ccp4=yes
      else
        CCP4_CXXFLAGS="-I$ccp4_prefix/lib/src -I$ccp4_prefix/include/ccp4"
        CCP4_LIBS="-L$ccp4_prefix/lib/src -L$ccp4_prefix/lib -lccp4c"
        AC_MSG_WARN([assuming library $CCP4_LIBS headers $CCP4_CXXFLAGS])
        have_ccp4=yes
      fi
    else
      AC_MSG_ERROR([Failed to find ccp4 libs])
    fi
  ]
)

AM_PATH_CLIPPER(,
     [ if test "x$clipper_prefix" != x ; then
     if test "x$ac_CLIPPER_CXXFLAGS" != x ; then
       case "$ac_CLIPPER_CXXFLAGS" in
        *include/clipper )
        ac_CLIPPER_LDOPTS=`echo "$ac_CLIPPER_CXXFLAGS" | sed s#include/clipper#lib#g`
        ;;
        *)
        ac_CLIPPER_LDOPTS="$ac_CLIPPER_CXXFLAGS"
        esac
        ac_CLIPPER_LDOPTS=`echo "$ac_CLIPPER_LDOPTS" | sed s#-I#-L#`
        CLIPPER_CXXFLAGS="$ac_CLIPPER_CXXFLAGS"
        CLIPPER_LIBS="$ac_CLIPPER_LDOPTS -lclipper-cif -lclipper-ccp4 -lclipper-contrib -lclipper-minimol -lclipper-mmdb -lclipper-core"
        AC_MSG_WARN([assuming library $CLIPPER_LIBS headers $CLIPPER_CXXFLAGS])
        have_clipper=yes
      else
        CLIPPER_CXXFLAGS="-I$clipper_prefix/lib/src"
        CLIPPER_LIBS="-L$clipper_prefix/lib/src -L$clipper_prefix/lib -lclipper-cif -lclipper-ccp4 -lclipper-contrib -lclipper-minimol -lclipper-mmdb -lclipper-core"
        AC_MSG_WARN([assuming library $CLIPPER_LIBS headers $CLIPPER_CXXFLAGS])
        have_clipper=yes
      fi
    else
     AC_MSG_ERROR([Failed to find clipper ])
    fi
   ]    
) 

# Checks for typedefs, structures, and compiler characteristics.
#
# Does clipper use these?  I'm guessing not, commented.
# AC_C_CONST
# AC_C_INLINE
# AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([floor pow rint sqrt])

# test for function in m
AC_CHECK_LIB([m], [sqrt]) 

AC_SUBST(AR_FLAGS)
AM_CONDITIONAL(DEBUG, test "$ac_cv_enable_debug" = yes)
AM_CONDITIONAL(AMDEP, test "x$enable_dependency_tracking" != xno)

AC_SUBST(chltofom_LDFLAGS)
AC_SUBST(cphasematch_LDFLAGS)
AC_SUBST(csfcalc_LDFLAGS)
AC_SUBST(cfft_LDFLAGS)
AC_SUBST(cinvfft_LDFLAGS)
AC_SUBST(csigmaa_LDFLAGS)
AC_SUBST(cecalc_LDFLAGS)
AC_SUBST(cphasecombine_LDFLAGS)
AC_SUBST(cpatterson_LDFLAGS)
AC_SUBST(caniso_LDFLAGS)
AC_SUBST(convert2mtz_LDFLAGS)
AC_SUBST(ctruncate_LDFLAGS)
AC_SUBST(csymmatch_LDFLAGS)
AC_SUBST(cpirate_LDFLAGS)
AC_SUBST(cparrot_LDFLAGS)
AC_SUBST(cbuccaneer_LDFLAGS)
AC_SUBST(csequins_LDFLAGS)
AC_SUBST(cmaplocal_LDFLAGS)
AC_SUBST(cmakereference_LDFLAGS)
AC_SUBST(cmodeltoseq_LDFLAGS)
AC_SUBST(cncsfrommodel_LDFLAGS)
AC_SUBST(csloop_LDFLAGS)

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/pirate/Makefile
                 src/parrot/Makefile
                 src/buccaneer/Makefile
                 src/sloop/Makefile
                 html/Makefile
                 ])

AC_OUTPUT
