#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.56])
AC_INIT([ffp], [3.19], [gsims1997@yahoo.com],[ffp],[http://sourceforge/projects/ffp-phylogeny])
#  Define authors of the package
AC_DEFINE([PACKAGE_AUTHORS], ["Gregory E. Sims"], [Authors of this package.])
AC_DEFINE([COPY_YEAR], ["2009-2012"], [Copyright year.])
AC_DEFINE([PACKAGE_URL], ["http://sourceforge/projects/ffp-phylogeny"], [Package URL.])
#  Insert copyright notice in configure script
AC_COPYRIGHT
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE(ffp, 3.19)
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])
AC_CANONICAL_HOST

case $host_os in
	mac*)
	AC_SUBST(ISOSX,[yes])
	;;
	*)
	;;
esac	

# Enable additional argument to configure
AC_ARG_ENABLE(gui,[ --disable-gui    Disable compilation of ffpgui, requires the perl Tk module ],
	      [case "${enableval}" in
	       	yes | no ) WITH_GUI="${enableval}" ;;
		*) AC_MSG_ERROR(Bad value ${enableval} for --disable-gui) ;;
	esac], [WITH_GUI="yes"]
)

AM_CONDITIONAL([WITH_GUI], [test "x$WITH_GUI" = "xyes"])
# Check if 32-bit -- then enable large file support
AC_SYS_LARGEFILE

# Checks for programs.
AC_PROG_CC
AC_PATH_PROG([MKTEMP], [mktemp],[no])
if test "$MKTEMP" = "no" ; then
   AC_MSG_ERROR('mktemp' missing from path. Please install, or add to path.)
fi
AC_PATH_PROG([GETOPT], [getopt],[no])
if test "$GETOPT" = "no" ; then
   AC_MSG_ERROR('getopt' missing from path. Please install, or add to path.)
fi

if test "x$WITH_GUI" = "xyes" ; then 
	# Checks for perl module DBI
	AC_MSG_CHECKING(for perl module Tk)
	if  `/usr/bin/env perl  -MTk -e 1 2> /dev/null`
	then
		AC_MSG_RESULT(yes)
	else
		AC_MSG_RESULT(no)
	# May need alternate instructions for cygwin.
		AC_MSG_ERROR(Perl module Tk missing. Install using: perl -MCPAN -e 'install Tk'. Or use configure with '--disable-gui' to disable ffpgui build. )
	fi
fi


# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])

# Checks for header files.
AC_CHECK_HEADERS([limits.h stdio.h time.h stdlib.h string.h ctype.h float.h getopt.h math.h stdbool.h])
# Check for struct members
AC_CHECK_MEMBERS([struct stat.st_blksize])
# check for size types
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE

# Zero out CFLAGS -- use the params in src/Makefile.am
CFLAGS=


# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memmove pow sqrt memset strupr])

AC_CONFIG_FILES([Makefile
                 doc/Makefile
                 examples/Makefile
                 man/Makefile
                 scripts/Makefile
                 src/Makefile
		 tests/Makefile])
AC_OUTPUT
