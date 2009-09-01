# aclocal.m4 generated automatically by aclocal 1.6.3 -*- Autoconf -*-

# Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002
# Free Software Foundation, Inc.
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

# Do all the work for Automake.                            -*- Autoconf -*-

# This macro actually does too much some checks are only needed if
# your package does certain things.  But this isn't really a big deal.

# Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002
# Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# serial 8

# There are a few dirty hacks below to avoid letting `AC_PROG_CC' be
# written in clear, in which case automake, when reading aclocal.m4,
# will think it sees a *use*, and therefore will trigger all it's
# C support machinery.  Also note that it means that autoscan, seeing
# CC etc. in the Makefile, will ask for an AC_PROG_CC use...


AC_PREREQ([2.52])

# Autoconf 2.50 wants to disallow AM_ names.  We explicitly allow
# the ones we care about.
m4_pattern_allow([^AM_[A-Z]+FLAGS$])dnl

# AM_INIT_AUTOMAKE(PACKAGE, VERSION, [NO-DEFINE])
# AM_INIT_AUTOMAKE([OPTIONS])
# -----------------------------------------------
# The call with PACKAGE and VERSION arguments is the old style
# call (pre autoconf-2.50), which is being phased out.  PACKAGE
# and VERSION should now be passed to AC_INIT and removed from
# the call to AM_INIT_AUTOMAKE.
# We support both call styles for the transition.  After
# the next Automake release, Autoconf can make the AC_INIT
# arguments mandatory, and then we can depend on a new Autoconf
# release and drop the old call support.
AC_DEFUN([AM_INIT_AUTOMAKE],
[AC_REQUIRE([AM_SET_CURRENT_AUTOMAKE_VERSION])dnl
 AC_REQUIRE([AC_PROG_INSTALL])dnl
# test to see if srcdir already configured
if test "`cd $srcdir && pwd`" != "`pwd`" &&
   test -f $srcdir/config.status; then
  AC_MSG_ERROR([source directory already configured; run "make distclean" there first])
fi

# Define the identity of the package.
dnl Distinguish between old-style and new-style calls.
m4_ifval([$2],
[m4_ifval([$3], [_AM_SET_OPTION([no-define])])dnl
 AC_SUBST([PACKAGE], [$1])dnl
 AC_SUBST([VERSION], [$2])],
[_AM_SET_OPTIONS([$1])dnl
 AC_SUBST([PACKAGE], [AC_PACKAGE_TARNAME])dnl
 AC_SUBST([VERSION], [AC_PACKAGE_VERSION])])dnl

_AM_IF_OPTION([no-define],,
[AC_DEFINE_UNQUOTED(PACKAGE, "$PACKAGE", [Name of package])
 AC_DEFINE_UNQUOTED(VERSION, "$VERSION", [Version number of package])])dnl

# Some tools Automake needs.
AC_REQUIRE([AM_SANITY_CHECK])dnl
AC_REQUIRE([AC_ARG_PROGRAM])dnl
AM_MISSING_PROG(ACLOCAL, aclocal-${am__api_version})
AM_MISSING_PROG(AUTOCONF, autoconf)
AM_MISSING_PROG(AUTOMAKE, automake-${am__api_version})
AM_MISSING_PROG(AUTOHEADER, autoheader)
AM_MISSING_PROG(MAKEINFO, makeinfo)
AM_MISSING_PROG(AMTAR, tar)
AM_PROG_INSTALL_SH
AM_PROG_INSTALL_STRIP
# We need awk for the "check" target.  The system "awk" is bad on
# some platforms.
AC_REQUIRE([AC_PROG_AWK])dnl
AC_REQUIRE([AC_PROG_MAKE_SET])dnl

_AM_IF_OPTION([no-dependencies],,
[AC_PROVIDE_IFELSE([AC_PROG_][CC],
                  [_AM_DEPENDENCIES(CC)],
                  [define([AC_PROG_][CC],
                          defn([AC_PROG_][CC])[_AM_DEPENDENCIES(CC)])])dnl
AC_PROVIDE_IFELSE([AC_PROG_][CXX],
                  [_AM_DEPENDENCIES(CXX)],
                  [define([AC_PROG_][CXX],
                          defn([AC_PROG_][CXX])[_AM_DEPENDENCIES(CXX)])])dnl
])
])

# Copyright 2002  Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA

# AM_AUTOMAKE_VERSION(VERSION)
# ----------------------------
# Automake X.Y traces this macro to ensure aclocal.m4 has been
# generated from the m4 files accompanying Automake X.Y.
AC_DEFUN([AM_AUTOMAKE_VERSION],[am__api_version="1.6"])

# AM_SET_CURRENT_AUTOMAKE_VERSION
# -------------------------------
# Call AM_AUTOMAKE_VERSION so it can be traced.
# This function is AC_REQUIREd by AC_INIT_AUTOMAKE.
AC_DEFUN([AM_SET_CURRENT_AUTOMAKE_VERSION],
	 [AM_AUTOMAKE_VERSION([1.6.3])])

# Helper functions for option handling.                    -*- Autoconf -*-

# Copyright 2001, 2002  Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# serial 2

# _AM_MANGLE_OPTION(NAME)
# -----------------------
AC_DEFUN([_AM_MANGLE_OPTION],
[[_AM_OPTION_]m4_bpatsubst($1, [[^a-zA-Z0-9_]], [_])])

# _AM_SET_OPTION(NAME)
# ------------------------------
# Set option NAME.  Presently that only means defining a flag for this option.
AC_DEFUN([_AM_SET_OPTION],
[m4_define(_AM_MANGLE_OPTION([$1]), 1)])

# _AM_SET_OPTIONS(OPTIONS)
# ----------------------------------
# OPTIONS is a space-separated list of Automake options.
AC_DEFUN([_AM_SET_OPTIONS],
[AC_FOREACH([_AM_Option], [$1], [_AM_SET_OPTION(_AM_Option)])])

# _AM_IF_OPTION(OPTION, IF-SET, [IF-NOT-SET])
# -------------------------------------------
# Execute IF-SET if OPTION is set, IF-NOT-SET otherwise.
AC_DEFUN([_AM_IF_OPTION],
[m4_ifset(_AM_MANGLE_OPTION([$1]), [$2], [$3])])

#
# Check to make sure that the build environment is sane.
#

# Copyright 1996, 1997, 2000, 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# serial 3

# AM_SANITY_CHECK
# ---------------
AC_DEFUN([AM_SANITY_CHECK],
[AC_MSG_CHECKING([whether build environment is sane])
# Just in case
sleep 1
echo timestamp > conftest.file
# Do `set' in a subshell so we don't clobber the current shell's
# arguments.  Must try -L first in case configure is actually a
# symlink; some systems play weird games with the mod time of symlinks
# (eg FreeBSD returns the mod time of the symlink's containing
# directory).
if (
   set X `ls -Lt $srcdir/configure conftest.file 2> /dev/null`
   if test "$[*]" = "X"; then
      # -L didn't work.
      set X `ls -t $srcdir/configure conftest.file`
   fi
   rm -f conftest.file
   if test "$[*]" != "X $srcdir/configure conftest.file" \
      && test "$[*]" != "X conftest.file $srcdir/configure"; then

      # If neither matched, then we have a broken ls.  This can happen
      # if, for instance, CONFIG_SHELL is bash and it inherits a
      # broken ls alias from the environment.  This has actually
      # happened.  Such a system could not be considered "sane".
      AC_MSG_ERROR([ls -t appears to fail.  Make sure there is not a broken
alias in your environment])
   fi

   test "$[2]" = conftest.file
   )
then
   # Ok.
   :
else
   AC_MSG_ERROR([newly created file is older than distributed files!
Check your system clock])
fi
AC_MSG_RESULT(yes)])

#  -*- Autoconf -*-


# Copyright 1997, 1999, 2000, 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# serial 3

# AM_MISSING_PROG(NAME, PROGRAM)
# ------------------------------
AC_DEFUN([AM_MISSING_PROG],
[AC_REQUIRE([AM_MISSING_HAS_RUN])
$1=${$1-"${am_missing_run}$2"}
AC_SUBST($1)])


# AM_MISSING_HAS_RUN
# ------------------
# Define MISSING if not defined so far and test if it supports --run.
# If it does, set am_missing_run to use it, otherwise, to nothing.
AC_DEFUN([AM_MISSING_HAS_RUN],
[AC_REQUIRE([AM_AUX_DIR_EXPAND])dnl
test x"${MISSING+set}" = xset || MISSING="\${SHELL} $am_aux_dir/missing"
# Use eval to expand $SHELL
if eval "$MISSING --run true"; then
  am_missing_run="$MISSING --run "
else
  am_missing_run=
  AC_MSG_WARN([`missing' script is too old or missing])
fi
])

# AM_AUX_DIR_EXPAND

# Copyright 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# For projects using AC_CONFIG_AUX_DIR([foo]), Autoconf sets
# $ac_aux_dir to `$srcdir/foo'.  In other projects, it is set to
# `$srcdir', `$srcdir/..', or `$srcdir/../..'.
#
# Of course, Automake must honor this variable whenever it calls a
# tool from the auxiliary directory.  The problem is that $srcdir (and
# therefore $ac_aux_dir as well) can be either absolute or relative,
# depending on how configure is run.  This is pretty annoying, since
# it makes $ac_aux_dir quite unusable in subdirectories: in the top
# source directory, any form will work fine, but in subdirectories a
# relative path needs to be adjusted first.
#
# $ac_aux_dir/missing
#    fails when called from a subdirectory if $ac_aux_dir is relative
# $top_srcdir/$ac_aux_dir/missing
#    fails if $ac_aux_dir is absolute,
#    fails when called from a subdirectory in a VPATH build with
#          a relative $ac_aux_dir
#
# The reason of the latter failure is that $top_srcdir and $ac_aux_dir
# are both prefixed by $srcdir.  In an in-source build this is usually
# harmless because $srcdir is `.', but things will broke when you
# start a VPATH build or use an absolute $srcdir.
#
# So we could use something similar to $top_srcdir/$ac_aux_dir/missing,
# iff we strip the leading $srcdir from $ac_aux_dir.  That would be:
#   am_aux_dir='\$(top_srcdir)/'`expr "$ac_aux_dir" : "$srcdir//*\(.*\)"`
# and then we would define $MISSING as
#   MISSING="\${SHELL} $am_aux_dir/missing"
# This will work as long as MISSING is not called from configure, because
# unfortunately $(top_srcdir) has no meaning in configure.
# However there are other variables, like CC, which are often used in
# configure, and could therefore not use this "fixed" $ac_aux_dir.
#
# Another solution, used here, is to always expand $ac_aux_dir to an
# absolute PATH.  The drawback is that using absolute paths prevent a
# configured tree to be moved without reconfiguration.

# Rely on autoconf to set up CDPATH properly.
AC_PREREQ([2.50])

AC_DEFUN([AM_AUX_DIR_EXPAND], [
# expand $ac_aux_dir to an absolute path
am_aux_dir=`cd $ac_aux_dir && pwd`
])

# AM_PROG_INSTALL_SH
# ------------------
# Define $install_sh.

# Copyright 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

AC_DEFUN([AM_PROG_INSTALL_SH],
[AC_REQUIRE([AM_AUX_DIR_EXPAND])dnl
install_sh=${install_sh-"$am_aux_dir/install-sh"}
AC_SUBST(install_sh)])

# AM_PROG_INSTALL_STRIP

# Copyright 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# One issue with vendor `install' (even GNU) is that you can't
# specify the program used to strip binaries.  This is especially
# annoying in cross-compiling environments, where the build's strip
# is unlikely to handle the host's binaries.
# Fortunately install-sh will honor a STRIPPROG variable, so we
# always use install-sh in `make install-strip', and initialize
# STRIPPROG with the value of the STRIP variable (set by the user).
AC_DEFUN([AM_PROG_INSTALL_STRIP],
[AC_REQUIRE([AM_PROG_INSTALL_SH])dnl
# Installed binaries are usually stripped using `strip' when the user
# run `make install-strip'.  However `strip' might not be the right
# tool to use in cross-compilation environments, therefore Automake
# will honor the `STRIP' environment variable to overrule this program.
dnl Don't test for $cross_compiling = yes, because it might be `maybe'.
if test "$cross_compiling" != no; then
  AC_CHECK_TOOL([STRIP], [strip], :)
fi
INSTALL_STRIP_PROGRAM="\${SHELL} \$(install_sh) -c -s"
AC_SUBST([INSTALL_STRIP_PROGRAM])])

# serial 4						-*- Autoconf -*-

# Copyright 1999, 2000, 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.


# There are a few dirty hacks below to avoid letting `AC_PROG_CC' be
# written in clear, in which case automake, when reading aclocal.m4,
# will think it sees a *use*, and therefore will trigger all it's
# C support machinery.  Also note that it means that autoscan, seeing
# CC etc. in the Makefile, will ask for an AC_PROG_CC use...



# _AM_DEPENDENCIES(NAME)
# ----------------------
# See how the compiler implements dependency checking.
# NAME is "CC", "CXX", "GCJ", or "OBJC".
# We try a few techniques and use that to set a single cache variable.
#
# We don't AC_REQUIRE the corresponding AC_PROG_CC since the latter was
# modified to invoke _AM_DEPENDENCIES(CC); we would have a circular
# dependency, and given that the user is not expected to run this macro,
# just rely on AC_PROG_CC.
AC_DEFUN([_AM_DEPENDENCIES],
[AC_REQUIRE([AM_SET_DEPDIR])dnl
AC_REQUIRE([AM_OUTPUT_DEPENDENCY_COMMANDS])dnl
AC_REQUIRE([AM_MAKE_INCLUDE])dnl
AC_REQUIRE([AM_DEP_TRACK])dnl

ifelse([$1], CC,   [depcc="$CC"   am_compiler_list=],
       [$1], CXX,  [depcc="$CXX"  am_compiler_list=],
       [$1], OBJC, [depcc="$OBJC" am_compiler_list='gcc3 gcc'],
       [$1], GCJ,  [depcc="$GCJ"  am_compiler_list='gcc3 gcc'],
                   [depcc="$$1"   am_compiler_list=])

AC_CACHE_CHECK([dependency style of $depcc],
               [am_cv_$1_dependencies_compiler_type],
[if test -z "$AMDEP_TRUE" && test -f "$am_depcomp"; then
  # We make a subdir and do the tests there.  Otherwise we can end up
  # making bogus files that we don't know about and never remove.  For
  # instance it was reported that on HP-UX the gcc test will end up
  # making a dummy file named `D' -- because `-MD' means `put the output
  # in D'.
  mkdir conftest.dir
  # Copy depcomp to subdir because otherwise we won't find it if we're
  # using a relative directory.
  cp "$am_depcomp" conftest.dir
  cd conftest.dir

  am_cv_$1_dependencies_compiler_type=none
  if test "$am_compiler_list" = ""; then
     am_compiler_list=`sed -n ['s/^#*\([a-zA-Z0-9]*\))$/\1/p'] < ./depcomp`
  fi
  for depmode in $am_compiler_list; do
    # We need to recreate these files for each test, as the compiler may
    # overwrite some of them when testing with obscure command lines.
    # This happens at least with the AIX C compiler.
    echo '#include "conftest.h"' > conftest.c
    echo 'int i;' > conftest.h
    echo "${am__include} ${am__quote}conftest.Po${am__quote}" > confmf

    case $depmode in
    nosideeffect)
      # after this tag, mechanisms are not by side-effect, so they'll
      # only be used when explicitly requested
      if test "x$enable_dependency_tracking" = xyes; then
	continue
      else
	break
      fi
      ;;
    none) break ;;
    esac
    # We check with `-c' and `-o' for the sake of the "dashmstdout"
    # mode.  It turns out that the SunPro C++ compiler does not properly
    # handle `-M -o', and we need to detect this.
    if depmode=$depmode \
       source=conftest.c object=conftest.o \
       depfile=conftest.Po tmpdepfile=conftest.TPo \
       $SHELL ./depcomp $depcc -c conftest.c -o conftest.o >/dev/null 2>&1 &&
       grep conftest.h conftest.Po > /dev/null 2>&1 &&
       ${MAKE-make} -s -f confmf > /dev/null 2>&1; then
      am_cv_$1_dependencies_compiler_type=$depmode
      break
    fi
  done

  cd ..
  rm -rf conftest.dir
else
  am_cv_$1_dependencies_compiler_type=none
fi
])
AC_SUBST([$1DEPMODE], [depmode=$am_cv_$1_dependencies_compiler_type])
])


# AM_SET_DEPDIR
# -------------
# Choose a directory name for dependency files.
# This macro is AC_REQUIREd in _AM_DEPENDENCIES
AC_DEFUN([AM_SET_DEPDIR],
[rm -f .deps 2>/dev/null
mkdir .deps 2>/dev/null
if test -d .deps; then
  DEPDIR=.deps
else
  # MS-DOS does not allow filenames that begin with a dot.
  DEPDIR=_deps
fi
rmdir .deps 2>/dev/null
AC_SUBST([DEPDIR])
])


# AM_DEP_TRACK
# ------------
AC_DEFUN([AM_DEP_TRACK],
[AC_ARG_ENABLE(dependency-tracking,
[  --disable-dependency-tracking Speeds up one-time builds
  --enable-dependency-tracking  Do not reject slow dependency extractors])
if test "x$enable_dependency_tracking" != xno; then
  am_depcomp="$ac_aux_dir/depcomp"
  AMDEPBACKSLASH='\'
fi
AM_CONDITIONAL([AMDEP], [test "x$enable_dependency_tracking" != xno])
AC_SUBST([AMDEPBACKSLASH])
])

# Generate code to set up dependency tracking.   -*- Autoconf -*-

# Copyright 1999, 2000, 2001, 2002 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

#serial 2

# _AM_OUTPUT_DEPENDENCY_COMMANDS
# ------------------------------
AC_DEFUN([_AM_OUTPUT_DEPENDENCY_COMMANDS],
[for mf in $CONFIG_FILES; do
  # Strip MF so we end up with the name of the file.
  mf=`echo "$mf" | sed -e 's/:.*$//'`
  # Check whether this is an Automake generated Makefile or not.
  # We used to match only the files named `Makefile.in', but
  # some people rename them; so instead we look at the file content.
  # Grep'ing the first line is not enough: some people post-process
  # each Makefile.in and add a new line on top of each file to say so.
  # So let's grep whole file.
  if grep '^#.*generated by automake' $mf > /dev/null 2>&1; then
    dirpart=`AS_DIRNAME("$mf")`
  else
    continue
  fi
  grep '^DEP_FILES *= *[[^ @%:@]]' < "$mf" > /dev/null || continue
  # Extract the definition of DEP_FILES from the Makefile without
  # running `make'.
  DEPDIR=`sed -n -e '/^DEPDIR = / s///p' < "$mf"`
  test -z "$DEPDIR" && continue
  # When using ansi2knr, U may be empty or an underscore; expand it
  U=`sed -n -e '/^U = / s///p' < "$mf"`
  test -d "$dirpart/$DEPDIR" || mkdir "$dirpart/$DEPDIR"
  # We invoke sed twice because it is the simplest approach to
  # changing $(DEPDIR) to its actual value in the expansion.
  for file in `sed -n -e '
    /^DEP_FILES = .*\\\\$/ {
      s/^DEP_FILES = //
      :loop
	s/\\\\$//
	p
	n
	/\\\\$/ b loop
      p
    }
    /^DEP_FILES = / s/^DEP_FILES = //p' < "$mf" | \
       sed -e 's/\$(DEPDIR)/'"$DEPDIR"'/g' -e 's/\$U/'"$U"'/g'`; do
    # Make sure the directory exists.
    test -f "$dirpart/$file" && continue
    fdir=`AS_DIRNAME(["$file"])`
    AS_MKDIR_P([$dirpart/$fdir])
    # echo "creating $dirpart/$file"
    echo '# dummy' > "$dirpart/$file"
  done
done
])# _AM_OUTPUT_DEPENDENCY_COMMANDS


# AM_OUTPUT_DEPENDENCY_COMMANDS
# -----------------------------
# This macro should only be invoked once -- use via AC_REQUIRE.
#
# This code is only required when automatic dependency tracking
# is enabled.  FIXME.  This creates each `.P' file that we will
# need in order to bootstrap the dependency handling code.
AC_DEFUN([AM_OUTPUT_DEPENDENCY_COMMANDS],
[AC_CONFIG_COMMANDS([depfiles],
     [test x"$AMDEP_TRUE" != x"" || _AM_OUTPUT_DEPENDENCY_COMMANDS],
     [AMDEP_TRUE="$AMDEP_TRUE" ac_aux_dir="$ac_aux_dir"])
])

# Copyright 2001 Free Software Foundation, Inc.             -*- Autoconf -*-

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# serial 2

# AM_MAKE_INCLUDE()
# -----------------
# Check to see how make treats includes.
AC_DEFUN([AM_MAKE_INCLUDE],
[am_make=${MAKE-make}
cat > confinc << 'END'
doit:
	@echo done
END
# If we don't find an include directive, just comment out the code.
AC_MSG_CHECKING([for style of include used by $am_make])
am__include="#"
am__quote=
_am_result=none
# First try GNU make style include.
echo "include confinc" > confmf
# We grep out `Entering directory' and `Leaving directory'
# messages which can occur if `w' ends up in MAKEFLAGS.
# In particular we don't look at `^make:' because GNU make might
# be invoked under some other name (usually "gmake"), in which
# case it prints its new name instead of `make'.
if test "`$am_make -s -f confmf 2> /dev/null | fgrep -v 'ing directory'`" = "done"; then
   am__include=include
   am__quote=
   _am_result=GNU
fi
# Now try BSD make style include.
if test "$am__include" = "#"; then
   echo '.include "confinc"' > confmf
   if test "`$am_make -s -f confmf 2> /dev/null`" = "done"; then
      am__include=.include
      am__quote="\""
      _am_result=BSD
   fi
fi
AC_SUBST(am__include)
AC_SUBST(am__quote)
AC_MSG_RESULT($_am_result)
rm -f confinc confmf
])

# AM_CONDITIONAL                                              -*- Autoconf -*-

# Copyright 1997, 2000, 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# serial 5

AC_PREREQ(2.52)

# AM_CONDITIONAL(NAME, SHELL-CONDITION)
# -------------------------------------
# Define a conditional.
AC_DEFUN([AM_CONDITIONAL],
[ifelse([$1], [TRUE],  [AC_FATAL([$0: invalid condition: $1])],
        [$1], [FALSE], [AC_FATAL([$0: invalid condition: $1])])dnl
AC_SUBST([$1_TRUE])
AC_SUBST([$1_FALSE])
if $2; then
  $1_TRUE=
  $1_FALSE='#'
else
  $1_TRUE='#'
  $1_FALSE=
fi
AC_CONFIG_COMMANDS_PRE(
[if test -z "${$1_TRUE}" && test -z "${$1_FALSE}"; then
  AC_MSG_ERROR([conditional \"$1\" was never defined.
Usually this means the macro was only invoked conditionally.])
fi])])

# Add --enable-maintainer-mode option to configure.
# From Jim Meyering

# Copyright 1996, 1998, 2000, 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# serial 1

AC_DEFUN([AM_MAINTAINER_MODE],
[AC_MSG_CHECKING([whether to enable maintainer-specific portions of Makefiles])
  dnl maintainer-mode is disabled by default
  AC_ARG_ENABLE(maintainer-mode,
[  --enable-maintainer-mode enable make rules and dependencies not useful
                          (and sometimes confusing) to the casual installer],
      USE_MAINTAINER_MODE=$enableval,
      USE_MAINTAINER_MODE=no)
  AC_MSG_RESULT([$USE_MAINTAINER_MODE])
  AM_CONDITIONAL(MAINTAINER_MODE, [test $USE_MAINTAINER_MODE = yes])
  MAINT=$MAINTAINER_MODE_TRUE
  AC_SUBST(MAINT)dnl
]
)

# ===========================================================================
#           http://www.nongnu.org/autoconf-archive/ax_pthread.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_PTHREAD([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro figures out how to build C programs using POSIX threads. It
#   sets the PTHREAD_LIBS output variable to the threads library and linker
#   flags, and the PTHREAD_CFLAGS output variable to any special C compiler
#   flags that are needed. (The user can also force certain compiler
#   flags/libs to be tested by setting these environment variables.)
#
#   Also sets PTHREAD_CC to any special C compiler that is needed for
#   multi-threaded programs (defaults to the value of CC otherwise). (This
#   is necessary on AIX to use the special cc_r compiler alias.)
#
#   NOTE: You are assumed to not only compile your program with these flags,
#   but also link it with them as well. e.g. you should link with
#   $PTHREAD_CC $CFLAGS $PTHREAD_CFLAGS $LDFLAGS ... $PTHREAD_LIBS $LIBS
#
#   If you are only building threads programs, you may wish to use these
#   variables in your default LIBS, CFLAGS, and CC:
#
#          LIBS="$PTHREAD_LIBS $LIBS"
#          CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
#          CC="$PTHREAD_CC"
#
#   In addition, if the PTHREAD_CREATE_JOINABLE thread-attribute constant
#   has a nonstandard name, defines PTHREAD_CREATE_JOINABLE to that name
#   (e.g. PTHREAD_CREATE_UNDETACHED on AIX).
#
#   ACTION-IF-FOUND is a list of shell commands to run if a threads library
#   is found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it
#   is not found. If ACTION-IF-FOUND is not specified, the default action
#   will define HAVE_PTHREAD.
#
#   Please let the authors know if this macro fails on any platform, or if
#   you have any other suggestions or comments. This macro was based on work
#   by SGJ on autoconf scripts for FFTW (http://www.fftw.org/) (with help
#   from M. Frigo), as well as ac_pthread and hb_pthread macros posted by
#   Alejandro Forero Cuervo to the autoconf macro repository. We are also
#   grateful for the helpful feedback of numerous users.
#
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

AC_DEFUN([AX_PTHREAD], [
AC_REQUIRE([AC_CANONICAL_HOST])
AC_LANG_SAVE
AC_LANG_C
ax_pthread_ok=no

# We used to check for pthread.h first, but this fails if pthread.h
# requires special compiler flags (e.g. on True64 or Sequent).
# It gets checked for in the link test anyway.

# First of all, check if the user has set any of the PTHREAD_LIBS,
# etcetera environment variables, and if threads linking works using
# them:
if test x"$PTHREAD_LIBS$PTHREAD_CFLAGS" != x; then
        save_CFLAGS="$CFLAGS"
        CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
        save_LIBS="$LIBS"
        LIBS="$PTHREAD_LIBS $LIBS"
        AC_MSG_CHECKING([for pthread_join in LIBS=$PTHREAD_LIBS with CFLAGS=$PTHREAD_CFLAGS])
        AC_TRY_LINK_FUNC(pthread_join, ax_pthread_ok=yes)
        AC_MSG_RESULT($ax_pthread_ok)
        if test x"$ax_pthread_ok" = xno; then
                PTHREAD_LIBS=""
                PTHREAD_CFLAGS=""
        fi
        LIBS="$save_LIBS"
        CFLAGS="$save_CFLAGS"
fi

# We must check for the threads library under a number of different
# names; the ordering is very important because some systems
# (e.g. DEC) have both -lpthread and -lpthreads, where one of the
# libraries is broken (non-POSIX).

# Create a list of thread flags to try.  Items starting with a "-" are
# C compiler flags, and other items are library names, except for "none"
# which indicates that we try without any flags at all, and "pthread-config"
# which is a program returning the flags for the Pth emulation library.

ax_pthread_flags="pthreads none -Kthread -kthread lthread -pthread -pthreads -mthreads pthread --thread-safe -mt pthread-config"

# The ordering *is* (sometimes) important.  Some notes on the
# individual items follow:

# pthreads: AIX (must check this before -lpthread)
# none: in case threads are in libc; should be tried before -Kthread and
#       other compiler flags to prevent continual compiler warnings
# -Kthread: Sequent (threads in libc, but -Kthread needed for pthread.h)
# -kthread: FreeBSD kernel threads (preferred to -pthread since SMP-able)
# lthread: LinuxThreads port on FreeBSD (also preferred to -pthread)
# -pthread: Linux/gcc (kernel threads), BSD/gcc (userland threads)
# -pthreads: Solaris/gcc
# -mthreads: Mingw32/gcc, Lynx/gcc
# -mt: Sun Workshop C (may only link SunOS threads [-lthread], but it
#      doesn't hurt to check since this sometimes defines pthreads too;
#      also defines -D_REENTRANT)
#      ... -mt is also the pthreads flag for HP/aCC
# pthread: Linux, etcetera
# --thread-safe: KAI C++
# pthread-config: use pthread-config program (for GNU Pth library)

case "${host_cpu}-${host_os}" in
        *solaris*)

        # On Solaris (at least, for some versions), libc contains stubbed
        # (non-functional) versions of the pthreads routines, so link-based
        # tests will erroneously succeed.  (We need to link with -pthreads/-mt/
        # -lpthread.)  (The stubs are missing pthread_cleanup_push, or rather
        # a function called by this macro, so we could check for that, but
        # who knows whether they'll stub that too in a future libc.)  So,
        # we'll just look for -pthreads and -lpthread first:

        ax_pthread_flags="-pthreads pthread -mt -pthread $ax_pthread_flags"
        ;;
esac

if test x"$ax_pthread_ok" = xno; then
for flag in $ax_pthread_flags; do

        case $flag in
                none)
                AC_MSG_CHECKING([whether pthreads work without any flags])
                ;;

                -*)
                AC_MSG_CHECKING([whether pthreads work with $flag])
                PTHREAD_CFLAGS="$flag"
                ;;

		pthread-config)
		AC_CHECK_PROG(ax_pthread_config, pthread-config, yes, no)
		if test x"$ax_pthread_config" = xno; then continue; fi
		PTHREAD_CFLAGS="`pthread-config --cflags`"
		PTHREAD_LIBS="`pthread-config --ldflags` `pthread-config --libs`"
		;;

                *)
                AC_MSG_CHECKING([for the pthreads library -l$flag])
                PTHREAD_LIBS="-l$flag"
                ;;
        esac

        save_LIBS="$LIBS"
        save_CFLAGS="$CFLAGS"
        LIBS="$PTHREAD_LIBS $LIBS"
        CFLAGS="$CFLAGS $PTHREAD_CFLAGS"

        # Check for various functions.  We must include pthread.h,
        # since some functions may be macros.  (On the Sequent, we
        # need a special flag -Kthread to make this header compile.)
        # We check for pthread_join because it is in -lpthread on IRIX
        # while pthread_create is in libc.  We check for pthread_attr_init
        # due to DEC craziness with -lpthreads.  We check for
        # pthread_cleanup_push because it is one of the few pthread
        # functions on Solaris that doesn't have a non-functional libc stub.
        # We try pthread_create on general principles.
        AC_TRY_LINK([#include <pthread.h>],
                    [pthread_t th; pthread_join(th, 0);
                     pthread_attr_init(0); pthread_cleanup_push(0, 0);
                     pthread_create(0,0,0,0); pthread_cleanup_pop(0); ],
                    [ax_pthread_ok=yes])

        LIBS="$save_LIBS"
        CFLAGS="$save_CFLAGS"

        AC_MSG_RESULT($ax_pthread_ok)
        if test "x$ax_pthread_ok" = xyes; then
                break;
        fi

        PTHREAD_LIBS=""
        PTHREAD_CFLAGS=""
done
fi

# Various other checks:
if test "x$ax_pthread_ok" = xyes; then
        save_LIBS="$LIBS"
        LIBS="$PTHREAD_LIBS $LIBS"
        save_CFLAGS="$CFLAGS"
        CFLAGS="$CFLAGS $PTHREAD_CFLAGS"

        # Detect AIX lossage: JOINABLE attribute is called UNDETACHED.
	AC_MSG_CHECKING([for joinable pthread attribute])
	attr_name=unknown
	for attr in PTHREAD_CREATE_JOINABLE PTHREAD_CREATE_UNDETACHED; do
	    AC_TRY_LINK([#include <pthread.h>], [int attr=$attr; return attr;],
                        [attr_name=$attr; break])
	done
        AC_MSG_RESULT($attr_name)
        if test "$attr_name" != PTHREAD_CREATE_JOINABLE; then
            AC_DEFINE_UNQUOTED(PTHREAD_CREATE_JOINABLE, $attr_name,
                               [Define to necessary symbol if this constant
                                uses a non-standard name on your system.])
        fi

        AC_MSG_CHECKING([if more special flags are required for pthreads])
        flag=no
        case "${host_cpu}-${host_os}" in
            *-aix* | *-freebsd* | *-darwin*) flag="-D_THREAD_SAFE";;
            *solaris* | *-osf* | *-hpux*) flag="-D_REENTRANT";;
        esac
        AC_MSG_RESULT(${flag})
        if test "x$flag" != xno; then
            PTHREAD_CFLAGS="$flag $PTHREAD_CFLAGS"
        fi

        LIBS="$save_LIBS"
        CFLAGS="$save_CFLAGS"

        # More AIX lossage: must compile with xlc_r or cc_r
	if test x"$GCC" != xyes; then
          AC_CHECK_PROGS(PTHREAD_CC, xlc_r cc_r, ${CC})
        else
          PTHREAD_CC=$CC
	fi
else
        PTHREAD_CC="$CC"
fi

AC_SUBST(PTHREAD_LIBS)
AC_SUBST(PTHREAD_CFLAGS)
AC_SUBST(PTHREAD_CC)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$ax_pthread_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_PTHREAD,1,[Define if you have POSIX threads libraries and header files.]),[$1])
        :
else
        ax_pthread_ok=no
        $2
fi
AC_LANG_RESTORE
])dnl AX_PTHREAD
# AC_PROGS_OPTIONS      
# ------------------
# target compilation optons
#
AC_DEFUN([AC_PROGS_OPTIONS],
[
#specifics for various machines
test "${target_os:+set}" = set || target_os="$host_os"
case "$target_os" in
  *osf* | *64* )
    if test "`basename $CXX`" = cxx; then
      case "$CXXFLAGS" in
        *strict_ansi* ) ;;
        * )
          CXXFLAGS="$CXXFLAGS -ieee -std strict_ansi -alternative_tokens -timplicit_local -no_implicit_include"
        esac
    fi
  ;;
  *linux* )
    if test "`basename $CXX`" = icpc ; then
      case `$CXX -dumpversion 2>&1` in
        11.0 )
          case "$CXXFLAGS" in
            *-fno-inline*) ;;
            *) CXXFLAGS="$CXXFLAGS -fno-inline"
          esac
        ;;
        *)
      esac
    fi
   ;;
  *irix* )
    if test "`basename $CXX`" = CC; then
      case "$CXXFLAGS" in
        *LANG:std* ) ;;
        * )
          CXXFLAGS="$CXXFLAGS -LANG:std"
        esac
    fi
  ;;
  *darwin* )
    if test "`basename $CXX`" = gcc || test "`basename $CXX`" = g++; then
      case `$CXX -v 2>&1` in 
       *3.1*)
# problem with PIC relocation tables for 3.1
        case "$CXXFLAGS" in
         *-O* | *-O1* | *-O2* | *-O3* )
          CXXFLAGS=`echo $CXXFLAGS | sed s%-O[[\ 123]]%-O0%g`
          ;;
         *-O0* ) ;;
         * )
          CXXFLAGS="$CXXFLAGS -O0"
        esac
        ;;
      *) 
      esac
    fi
  ;;
  *solaris* )
    if test "`basename $CXX`" = CC; then
      AR=CC
      AR_FLAGS="-xar -o"
    fi
  ;;
  * )
esac

] )  # 


# AM_PATH_FFTW([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])
# ----------------------------------------------------------
# set up for FFTW
#
AC_DEFUN([AM_PATH_FFTW],
[
AC_PROVIDE([AM_PATH_FFTW])

AC_ARG_WITH(fftw,
 AC_HELP_STRING([--with-fftw=PFX], [Prefix where FFTW has been installed] ),
 [
    test "$withval" = no && AC_MSG_ERROR([fftw is a required package])
    test "$withval" = yes || fftw_prefix="$withval" 
    with_fftw=yes ],
 [ with_fftw=yes ] )

if test $with_fftw = yes ; then
#user override
AS_IF([test "x$FFTW_LIBS" != x && test "x$FFTW_CXXFLAGS" != x ],
[
  have_fftw=yes
],
[
saved_LIBS="$LIBS"
saved_CXXFLAGS="$CXXFLAGS"
FFTW_LIBS=""
FFTW_CXXFLAGS=""

if test x$fftw_prefix != x; then
	# very likely the majority of cases, we will have been configured with:
	# --with-fftw=/some/thing
	#

 # should be ac_FFTW_CXXFLAGS="-I$FFTW_prefix/include"
 #  
 ac_FFTW_CXXFLAGS="-I$fftw_prefix/include"
 #
 # Similarly for fftw, the uninstalled library position is simply in
 # $fftw_prefix, but the installed is in the standard prefixed subdirectory.
 #
 # SGI compiler CC (CXX=CC) needs -lm to link maths library, but 
 # GCC c++ does not.
 #
 ac_FFTW_LDOPTS="-L$fftw_prefix/lib"
else
 # the compiler looks in the "standard" places for FFTW.  In real life,
 # it would be quite possible that FFTW would not be installed in
 # /usr/include, /usr/lib etc. so the defaults will not usually find
 # the right dependencies.
 ac_FFTW_CXXFLAGS=""
 ac_FFTW_LDOPTS=""
fi #dnl test fftw_prefix

fftwname="fftw"
rfftwname="rfftw"

AC_MSG_CHECKING([for fftw_print_max_memory_usage in $fftwname])

	LIBS="$ac_FFTW_LDOPTS $saved_LIBS -l$rfftwname -l$fftwname"
	CXXFLAGS="$ac_FFTW_CXXFLAGS $saved_CXXFLAGS"
	#
	# AC_TRY_LINK uses the c compiler (set by AC_LANG), so we will
	# temporarily reassign $CC to the c++ compiler.
 	#
	AC_LANG_PUSH(C++)
	AC_TRY_LINK([#include <$fftwname.h>] ,[  fftw_print_max_memory_usage();  ], have_fftw=yes, have_fftw=no)
	if test x$have_fftw=xyes; then
	   AC_TRY_LINK(
[#include <$fftwname.h>] ,[  
   fftw_real *fftwp = 0;
   float *fftp = 0;
   fftp = fftwp;
          ], 
           have_fftw=yes, have_fftw=no)
	fi
	AC_MSG_RESULT($have_fftw)

if test $have_fftw = no; then

  fftwname="sfftw"
  rfftwname="srfftw"
  AC_MSG_CHECKING([for fftw_print_max_memory_usage in $fftwname])
  LIBS="$ac_FFTW_LDOPTS $saved_LIBS -l$rfftwname -l$fftwname"
  CXXFLAGS="$ac_FFTW_CXXFLAGS $saved_CXXFLAGS"
  AC_TRY_LINK([#include <$fftwname.h>] ,[  fftw_print_max_memory_usage();  ], have_fftw=yes, have_fftw=no)
  if test x$have_fftw=xyes; then
    AC_TRY_LINK(
      [#include <$fftwname.h>] ,[
   fftw_real *fftwp = 0;
   float *fftp = 0;
   fftp = fftwp;
      ],
      have_fftw=yes, have_fftw=no)
    fi
  AC_MSG_RESULT($have_fftw)
fi

 AC_LANG_POP(C++)
 LIBS="$saved_LIBS"
 CXXFLAGS="$saved_CXXFLAGS"

]) #dnl user override

AS_IF([test x$have_fftw = xyes],
  [
    test "x$FFTW_CXXFLAGS" = x && FFTW_CXXFLAGS="$ac_FFTW_CXXFLAGS"
    test "x$FFTW_LIBS" = x && FFTW_LIBS="$ac_FFTW_LDOPTS -l$rfftwname -l$fftwname"
    ifelse([$1], , :, [$1])
  ],
  [
    AC_MSG_ERROR([If fftw exist on you system, are you sure you are using the fftw libraries that was configured with --enable-float?])
    ifelse([$2], , :, [$2])
  ])

fi # --with-fftw

AC_SUBST(FFTW_CXXFLAGS)
AC_SUBST(FFTW_LIBS)

])


# Note that mmdb as distributed by Eugene does not install and does
# conform to normal unix standard notions of the behaviour of include 
# and libraries.
#
# So we will depend on you doing something by hand, and that is to link
# (or copy) mmdb.a to libmmdb.a (and that it is in the top directory of 
# mmdb).
#

 
# AM_PATH_MMDB([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])
# ----------------------------------------------------------
# set up for MMDB
#
AC_DEFUN([AM_PATH_MMDB],
[
AC_PROVIDE([AM_PATH_MMDB])

AC_ARG_WITH(mmdb,
  AC_HELP_STRING( [--with-mmdb=PFX], [use mmdb library (default NO) and set prefix] ),
  [
    test "$withval" = no || with_mmdb=yes 
    test "$withval" = yes || mmdb_prefix="$withval" ],
  [ with_mmdb="$enable_mmdb" 
    test $enable_mmdbold = yes && with_mmdb=yes
    test $enable_cif = yes && with_mmdb=yes
    test $enable_minimol = yes && with_mmdb=yes] ) 

if test x$with_mmdb = xyes ; then
#user override
AS_IF([test "x$MMDB_LIBS" != x && test "x$MMDB_CXXFLAGS" != x ],
[
  have_mmdb=yes
],
[
saved_LIBS="$LIBS"
saved_CXXFLAGS="$CXXFLAGS"
MMDB_CXXFLAGS=""
MMDB_LIBS=""

if test x$mmdb_prefix != x; then
	# very likely the majority of cases, we will try to configure with:
	# --with-mmdb=/some/thing
	#


 # should ideally be MMDB_CXXFLAGS="-I$MMDB_prefix/include", and the like
 # when MMDB and dependencies get installed
 #  
ac_mmdb_dirs='
.
include
lib
src
lib/src
lib/src/mmdb'
for ac_dir in $ac_mmdb_dirs; do
  if test -r "$mmdb_prefix/$ac_dir/mmdb/mmdb_manager.h"; then
    ac_MMDB_CXXFLAGS="-I$mmdb_prefix/$ac_dir"
    break
    fi
  done
 #
 # SGI compiler CC (CXX=CC) needs -lm to link maths library, but 
 # GCC c++ does not.
 #
for ac_dir in $ac_mmdb_dirs; do
  for ac_extension in a so sl dylib; do
  if test -r "$mmdb_prefix/$ac_dir/libmmdb.$ac_extension"; then
    ac_MMDB_LDOPTS="-L$mmdb_prefix/$ac_dir -lmmdb"
    break 2
    fi
  done
  done

else
 # the compiler looks in the "standard" places for MMDB.  In real life,
 # it would be quite unlikely that MMDB would be installed in /usr/include, 
 # /usr/lib etc. so this code will not usually find the right dependencies.
 ac_MMDB_CXXFLAGS=""
 ac_MMDB_LDOPTS=""
fi

AC_MSG_CHECKING([for CMMDBManager in MMDB])
	LIBS="$ac_MMDB_LDOPTS $saved_LIBS"
	CXXFLAGS="$ac_MMDB_CXXFLAGS $saved_CXXFLAGS"
	#
	# AC_TRY_LINK uses the c compiler (set by AC_LANG), so we will
	# temporarily reassign $CC to the c++ compiler.
 	#
	AC_LANG_PUSH(C++)
	AC_TRY_LINK([#include "mmdb/mmdb_manager.h"] ,[ CMMDBManager a;  ], have_mmdb=yes, have_mmdb=no)
	AC_LANG_POP(C++)  # the language we have just quit

 LIBS="$saved_LIBS"
 CXXFLAGS="$saved_CXXFLAGS"

]) # user override

AS_IF([test "x$have_mmdb" = xyes],
  [  test "x$MMDB_CXXFLAGS" = x && MMDB_CXXFLAGS=$ac_MMDB_CXXFLAGS
     test "x$MMDB_LIBS" = x && MMDB_LIBS=$ac_MMDB_LDOPTS
     AC_MSG_RESULT($have_mmdb)
     ifelse([$1], , :, [$1])],
  [  AC_MSG_RESULT($have_mmdb)
     ifelse([$2], , :, [$2])] )

fi #dnl --with-mmdb 

AC_SUBST(MMDB_CXXFLAGS)
AC_SUBST(MMDB_LIBS)
 
])


# AM_PATH_CCP4([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])
AC_DEFUN([AM_PATH_CCP4],
[
AC_PROVIDE([AM_PATH_CCP4])

AC_ARG_WITH(ccp4,
  AC_HELP_STRING( [--with-ccp4=PFX], [location of ccp4c] ),
  [
    test "$withval" = no || with_ccp4=yes 
    test "$withval" = yes || ccp4_prefix="$withval" ],
  [ with_ccp4=yes ] ) 

if test x$with_ccp4 = xyes ; then
#user override
AS_IF([test "x$CCP4_LIBS" != x && test "x$CCP4_CXXFLAGS" != x ],
[
  have_ccp4=yes
],
[
AC_MSG_CHECKING([for ccp4_errno in CCP4])

saved_LIBS="$LIBS"
saved_CXXFLAGS="$CXXFLAGS"
CCP4_LIBS=""
CCP4_CXXFLAGS=""

if test "x$ccp4_prefix" != x; then
ac_ccp4_dirs='
.
include
include/ccp4
lib
lib/src'
for ac_dir in $ac_ccp4_dirs; do
  if test -r "$ccp4_prefix/$ac_dir/ccp4_errno.h"; then
    ac_CCP4_CXXFLAGS="-I$ccp4_prefix/$ac_dir"
    break
    fi
  done
for ac_dir in $ac_ccp4_dirs; do
  for ac_extension in a so sl dylib; do
  if test -r "$ccp4_prefix/$ac_dir/libccp4c.$ac_extension"; then
    ac_CCP4_LDOPTS="-L$ccp4_prefix/$ac_dir -lccp4c"
    break 2
    fi
  done
  done
else
 ac_CCP4_CXXFLAGS=""
 ac_CCP4_LDOPTS="-lccp4c"
fi


LIBS="$ac_CCP4_LDOPTS $saved_LIBS"
CXXFLAGS="$ac_CCP4_CXXFLAGS $saved_CXXFLAGS"
#
# AC_TRY_LINK uses the c compiler (set by AC_LANG), so we will
# temporarily reassign $CC to the c++ compiler.
#
AC_LANG_PUSH(C++)
AC_TRY_LINK([#include "ccp4_errno.h"],
  [int a = ccp4_errno;  CCP4::ccp4_error("conftest"); ], have_ccp4=yes, have_ccp4=no)
AC_LANG_POP(C++)  # the language we have just quit
AC_MSG_RESULT($have_ccp4)

 LIBS="$saved_LIBS"
 CXXFLAGS="$saved_CXXFLAGS"

]) # user override

AS_IF([test x$have_ccp4 = xyes],
 [
   test "x$CCP4_CXXFLAGS" = x && CCP4_CXXFLAGS="$ac_CCP4_CXXFLAGS"
   test "x$CCP4_LIBS" = x && CCP4_LIBS="$ac_CCP4_LDOPTS"
   ifelse([$1], , :, [$1]) ],
 [
   ifelse([$2], , :, [$2]) ]
)

fi #dnl --with-ccp4

AC_SUBST(CCP4_CXXFLAGS)
AC_SUBST(CCP4_LIBS)
])



# AM_PATH_CLIPPER([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])
AC_DEFUN([AM_PATH_CLIPPER],
[
AC_PROVIDE([AM_PATH_CLIPPER])

AC_ARG_WITH(clipper,
  AC_HELP_STRING( [--with-clipper=PFX], [location of clipper] ),
  [
    test "$withval" = no || with_clipper=yes 
    test "$withval" = yes || clipper_prefix="$withval" ],
  [ with_clipper=yes ] ) 

if test x$with_clipper = xyes ; then
#user override
AS_IF([test "x$CLIPPER_LIBS" != x && test "x$CLIPPER_CXXFLAGS" != x ],
[
  have_clipper=yes
],
[
AC_MSG_CHECKING([for clipper_errno in CLIPPER])

saved_LIBS="$LIBS"
saved_CXXFLAGS="$CXXFLAGS"
CLIPPER_LIBS=""
CLIPPER_CXXFLAGS=""

if test "x$clipper_prefix" != x; then
ac_clipper_dirs='
.
include
include/clipper
lib/clipper/clipper/
lib'
for ac_dir in $ac_clipper_dirs; do
  if test -r "$clipper_prefix/$ac_dir/clipper/clipper.h"; then
    ac_CLIPPER_CXXFLAGS="-I$clipper_prefix/$ac_dir"
    break
    fi
  done
for ac_dir in $ac_clipper_dirs; do
  for ac_extension in a so sl dylib; do
  if test -r "$clipper_prefix/$ac_dir/libclipper-core.$ac_extension"; then
    ac_CLIPPER_LDOPTS="-L$clipper_prefix/$ac_dir -lclipper-cif -lclipper-ccp4 -lclipper-contrib -lclipper-minimol -lclipper-mmdb -lclipper-core"
    break 2
    fi
  done
  done
else
 ac_CLIPPER_CXXFLAGS=""
 ac_CLIPPER_LDOPTS="-lclipper-cif -lclipper-ccp4 -lclipper-contrib -lclipper-minimol -lclipper-mmdb -lclipper-core"
fi


LIBS="$ac_CLIPPER_LDOPTS $saved_LIBS"
CXXFLAGS="$ac_CLIPPER_CXXFLAGS $saved_CXXFLAGS"
#
# AC_TRY_LINK uses the c compiler (set by AC_LANG), so we will
# temporarily reassign $CC to the c++ compiler.
#
AC_LANG_PUSH(C++)
AC_TRY_LINK([#include "clipper/clipper.h"],
  [clipper::String a("test");  ], have_clipper=yes, have_clipper=no)
AC_LANG_POP(C++)  # the language we have just quit
AC_MSG_RESULT($have_clipper)

 LIBS="$saved_LIBS"
 CXXFLAGS="$saved_CXXFLAGS"

]) # user override

AS_IF([test x$have_clipper = xyes],
 [
   test "x$CLIPPER_CXXFLAGS" = x && CLIPPER_CXXFLAGS="$ac_CLIPPER_CXXFLAGS"
   test "x$CLIPPER_LIBS" = x && CLIPPER_LIBS="$ac_CLIPPER_LDOPTS"
   ifelse([$1], , :, [$1]) ],
 [
   ifelse([$2], , :, [$2]) ]
)

fi #dnl --with-clipper

AC_SUBST(CLIPPER_CXXFLAGS)
AC_SUBST(CLIPPER_LIBS)
])


