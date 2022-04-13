#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# Original configurations provided by Dr. Stéphane Guindon. 
# Modifications done by Marcelo S Zanetti.

AC_PREREQ(2.61)
AC_INIT([IgPHYML],[1.1.4 04022],[kenneth.hoehn@gmail.com])
AC_CONFIG_SRCDIR([src/simu.c]) #configure source directory
AC_CONFIG_HEADER([config.h]) #create header file

AM_INIT_AUTOMAKE
AC_CANONICAL_HOST
test -e "/etc/debian_version" && echo 1 || echo 0

# Checks for programs.
AC_PROG_CC

# Checks for m and log libraries.
AC_CHECK_LIB([m], [log])

# Checks for stdlib.h, stdarg.h, string.h, and float.h header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_SETVBUF_REVERSED
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([floor pow rint sqrt strchr strstr])

# Added for cross-compilation (thanks to Alban Mancheron)
case $host_os in
  darwin*)  CFLAGS+=" -I$MACOS_SDK/usr/include/gcc/  \
                      -I$MACOS_BASE/include/  \
                      -arch i386 -arch ppc  \
                      -Wno-long-long -Wno-variadic-macros"
            AC_DEFINE([MACOS], [1], [Running on a mac?])
            AM_CONDITIONAL([MACOS], [test 1 = 1]);;
dnl             LDFLAGS+=" -arch i386 -arch ppc";;
  mingw32*) CFLAGS+=" -std=gnu99"
            AM_CONDITIONAL([MACOS], [test 1 = 0]);;
  linux*|*) AM_CONDITIONAL([MACOS], [test 1 = 0]);;
esac

# End of cross-compilation macros
# CFLAGS 
# -g                   = debugging
# -pg                  = profiling
# -fopenmp             = openMP support 
# -O3                  = optimization level (fast option) 
# -Wfloat-equal        = Warn about = and <> comparisons of real numbers. 
# -fomit-frame-pointer = Don't keep the frame pointer in a register for functions that don't need one. This avoids the instructions to save, set up and restore frame pointers; it also makes an extra register available in many functions. It also makes debugging impossible on some machines.
# -funroll-loops       = Unroll loops whose number of iterations can be determined at compile time or upon entry to the loop. -funroll-loops implies -frerun-cse-after-loop. It also turns on complete loop peeling (i.e. complete removal of loops with small constant number of iterations). This option makes code larger, and may or may not make it run faster. 


# Recommended settings
CFLAGS="-O3 -fomit-frame-pointer -funroll-loops -fopenmp"
#CFLAGS="-O3 -fomit-frame-pointer -funroll-loops"

AC_ARG_ENABLE([blas_omp],[  --enable-blas_omp    Compile blas_omp],[blas_omp=yes])
AM_CONDITIONAL([WANT_BLAS_OMP], [test "$blas_omp" = yes]) 
if test "$blas_omp" == "yes"
then
    AC_DEFINE([BLAS_OMP], [1], [Define if BLAS and OMP is used.])
    #AC_DEFINE([TEST2], [1], [Define if TEST2 is used.])
    BLASFLAGS=''
  my_save_cflags="$CFLAGS"
  CFLAGS='-L/usr/lib64/atlas/ -L/usr/lib/atlas-base/ -llapack'
  AM_CFLAGS=''
  AC_MSG_CHECKING([whether CC supports -llapack])
  AC_LINK_IFELSE([AC_LANG_SOURCE([[int main() { return 0; }]])],
    [AC_MSG_RESULT([yes])]
    [AM_CFLAGS=' -llapack'],
    [AC_MSG_RESULT([no])]) 
  CFLAGS="$my_save_cflags"
  BLASFLAGS+="$AM_CFLAGS"
  
  my_save_cflags="$CFLAGS"
  CFLAGS=' -L/usr/lib64/atlas/ -L/usr/lib/atlas-base/ -lblas'
  AM_CFLAGS=''
  AC_MSG_CHECKING([whether CC supports -lblas])
  AC_LINK_IFELSE([AC_LANG_SOURCE([[int main() { return 0; }]])],
    [AC_MSG_RESULT([yes])]
    [AM_CFLAGS=' -lblas'],
    [AC_MSG_RESULT([no])])
  CFLAGS="$my_save_cflags"
  BLASFLAGS+="$AM_CFLAGS"
  
  my_save_cflags="$CFLAGS"
  CFLAGS='-L/usr/lib64/atlas/ -L/usr/lib/atlas-base/ -lcblas'
  AM_CFLAGS=''
  AC_MSG_CHECKING([whether CC supports -lcblas])
  AC_LINK_IFELSE([AC_LANG_SOURCE([[int main() { return 0; }]])],
    [AC_MSG_RESULT([yes])]
    [AM_CFLAGS=' -lcblas'],
    [AC_MSG_RESULT([no])])
  CFLAGS="$my_save_cflags"
  BLASFLAGS+="$AM_CFLAGS"
  #CFLAGS+=$BLASFLAGS

  my_save_cflags="$CFLAGS"
  CFLAGS=' -L/usr/lib64/atlas/ -L/usr/lib/atlas-base/ -latlas'
  AM_CFLAGS=''
  AC_MSG_CHECKING([whether CC supports -latlas])
  AC_LINK_IFELSE([AC_LANG_SOURCE([[int main() { return 0; }]])],
    [AC_MSG_RESULT([yes])]
    [AM_CFLAGS=' -latlas'],
    [AC_MSG_RESULT([no])]) 
  CFLAGS="$my_save_cflags"
  BLASFLAGS+="$AM_CFLAGS"
  
  my_save_cflags="$CFLAGS"
  CFLAGS=' -L/usr/lib64/atlas/ -L/usr/lib/atlas-base/ -lsatlas'
  AM_CFLAGS=''
  AC_MSG_CHECKING([whether CC supports -lsatlas])
  AC_LINK_IFELSE([AC_LANG_SOURCE([[int main() { return 0; }]])],
    [AC_MSG_RESULT([yes])]
    [AM_CFLAGS=' -lsatlas'],
    [AC_MSG_RESULT([no])]) 
  CFLAGS="$my_save_cflags"
  BLASFLAGS+="$AM_CFLAGS"

  AC_SUBST(igphyml_LDADD,$BLASFLAGS)
fi

AC_ARG_ENABLE([omp],[  --enable-omp    Compile omp],[omp=yes])
AM_CONDITIONAL([WANT_OMP], [test "$omp" = yes]) 
if test "$omp" == "yes"
then
    AC_DEFINE([OMP], [1], [Define if OMP is used.])
fi

#AC_ARG_ENABLE([blas],[  --enable-blas    Compile blas],[blas=yes])
#AM_CONDITIONAL([WANT_BLAS], [test "$blas" = yes]) 
#if test "$blas" == "yes"
#then
#    AC_DEFINE([BLAS], [1], [Define if BLAS is used.])
#fi

AC_ARG_ENABLE([fedora],[  --enable-fedora    Compile Fedora],[fedora=yes])
AM_CONDITIONAL([FEDORA], [test "$fedora" = yes]) 

AC_CONFIG_FILES([Makefile src/Makefile])
echo $CFLAGS
AC_OUTPUT
