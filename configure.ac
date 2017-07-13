AC_INIT([HomomorphicEncryption], 0.1)

dnl Louis Aslett (aslett@stats.ox.ac.uk)
dnl September 2014

dnl Allow user to specify location of libraries etc
gmp_home=""
AC_ARG_WITH([gmp],
            AC_HELP_STRING([--with-gmp=MGR],
                           [specify the folder containing GMP include and lib directories, e.g. /usr/local]),
            [gmp=$withval])
if test [ -n "$gmp" ] ; then
GMP_HOME="${gmp}"
else
  if test [ -n "${GMP_HOME}" ] ; then
GMP_HOME="${GMP_HOME}"
else
  if test [ -e "/usr/local/lib/libgmp.so" -o -e "/usr/local/lib/libgmp.dylib" ] ; then
GMP_HOME="/usr/local"
else
  if test [ -e "/usr/lib/libgmp.so" -o -e "/usr/lib/libgmp.so" ] ; then
GMP_HOME="/usr"
fi
fi
fi
fi

mpfr_home=""
AC_ARG_WITH([mpfr],
            AC_HELP_STRING([--with-mpfr=MGR],
                           [specify the folder containing MPFR include and lib directories, e.g. /usr/local]),
            [mpfr=$withval])
if test [ -n "$mpfr" ] ; then
MPFR_HOME="${mpfr}"
else
  if test [ -n "${MPFR_HOME}" ] ; then
MPFR_HOME="${MPFR_HOME}"
else
  if test [ -e "/usr/local/lib/libmpfr.so" -o -e "/usr/local/lib/libmpfr.dylib" ] ; then
MPFR_HOME="/usr/local"
else
  if test [ -e "/usr/lib/libmpfr.so" -o -e "/usr/lib/libmpfr.dylib" ] ; then
MPFR_HOME="/usr"
fi
fi
fi
fi

flint_home=""
AC_ARG_WITH([flint],
            AC_HELP_STRING([--with-flint=MGR],
                           [specify the folder containing FLINT include and lib directories, e.g. /usr/local]),
            [flint=$withval])
if test [ -n "$flint" ] ; then
FLINT_HOME="${flint}"
else
  if test [ -n "${FLINT_HOME}" ] ; then
FLINT_HOME="${FLINT_HOME}"
else
  if test [ -e "/usr/local/include/flint/flint.h" ] ; then
FLINT_HOME="/usr/local"
else
  if test [ -e "/usr/include/flint/flint.h" ] ; then
FLINT_HOME="/usr"
fi
fi
fi
fi

dnl Catch non lib directories
dnl Ubuntu etc
if test [ -e "${FLINT_HOME}/lib/libflint.so" -o -e "${FLINT_HOME}/lib/libflint.dylib" ] ; then
  LIBDIR="lib"
else
  dnl Fedora
  if test [ -e "${FLINT_HOME}/lib64/libflint.so" -o -e "${FLINT_HOME}/lib64/libflint.dylib" ] ; then
    LIBDIR="lib64"
  fi
fi

dnl Check for not universally available functions
AC_CHECK_FUNCS([malloc_stats])
AC_CHECK_HEADERS([malloc.h])
AC_CONFIG_HEADER(config.h)

dnl Name the substitution variable and instruct configure to construct Makevars
dnl from Makevars.in
AC_SUBST(GMP_HOME)
AC_SUBST(MPFR_HOME)
AC_SUBST(FLINT_HOME)
AC_SUBST(LIBDIR)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
