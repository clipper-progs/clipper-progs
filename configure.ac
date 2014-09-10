
AC_INIT(clipper-progs, 0.1.0)

AC_CONFIG_AUX_DIR(build-aux)
AC_CONFIG_MACRO_DIR(m4)
AC_CONFIG_SRCDIR(src/buccaneer/buccaneer-build.cpp)
AC_CONFIG_HEADERS(config.h)
AM_INIT_AUTOMAKE([1.11 foreign subdir-objects silent-rules -Wall])

AM_MAINTAINER_MODE dnl this disables (by default) maintainer mode

AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL
AM_PROG_CC_C_O

AC_PROGS_OPTIONS

# pthreads are used by clipper only on Unix. We need to know pthread flags
# for linking with static clipper library.
case $host_os in
  cygwin* | mingw* | pw32* | cegcc*)
    ;;
  *)
    AX_PTHREAD
    ;;
esac

AC_SEARCH_LIBS(cos, m, , AC_MSG_ERROR([need math library]))

# ftplib (included in src/pirate) uses WinSock 2 on Windows
SOCKET_LIB=
AC_SEARCH_LIBS(select, , has_select=yes)
if test x$has_select != xyes ; then
    # 'main' must be used, AC_CHECK_LIB doesn't work well with Windows libs
    AC_CHECK_LIB(ws2_32, main, SOCKET_LIB="-lws2_32",
                               AC_MSG_ERROR([socked library not found.]))
fi

dnl Check for FFTW.
AM_PATH_FFTW(, AC_MSG_ERROR([locating necessary fftw]) )

AC_LANG_PUSH([C++])
AC_CHECK_LIB(ccp4c, ccp4_banner, :,
             AC_MSG_WARN([ccp4c library not found.]))
AC_CHECK_HEADER(ccp4/ccp4_general.h, :,
                AC_MSG_ERROR([libccp4 headers not found.]))

AC_CHECK_LIB(mmdb2, MMDB_CCIF_Init, :,
             AC_MSG_ERROR([mmdb2 library not found.]))
AC_CHECK_HEADER(mmdb2/mmdb_manager.h, :,
                AC_MSG_ERROR([mmdb2 headers not found.]))

AC_CHECK_HEADERS([clipper/clipper.h \
                  clipper/clipper-mmdb.h \
                  clipper/clipper-minimol.h \
                  clipper/clipper-contrib.h \
                  clipper/clipper-ccp4.h \
                  clipper/clipper-cif.h], :,
                AC_MSG_ERROR([clipper headers not found.]))
saved_LIBS="$LIBS"
AC_MSG_CHECKING(for class clipper::String in -lclipper-core)
LIBS="-lclipper-core $LIBS"
AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <clipper/clipper.h>],
                                [clipper::String foo])],
               AC_MSG_RESULT(yes),
               [AC_MSG_RESULT(no)
                AC_MSG_WARN([clipper library not found.])])
LIBS="$saved_LIBS"
AC_LANG_POP([C++])

AC_SUBST(SOCKET_LIB)

AC_OUTPUT(Makefile)
