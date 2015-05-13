AC_INIT(Piranha, version-1.2.1)

# Allow the user to specify where the BAMTools headers and libraries are
# located.
AC_ARG_WITH(stuff, [  --with-bam_tools_headers            specify the location of BAMTools headers if they're not in a default search path])
AC_ARG_WITH(stuff, [  --with-bam_tools_library            specify the location of the BAMTools libraries if they're not in a default search path])

# do tilde expansion on the paths given by the user and echo them back
if test "$with_bam_tools_headers" != ""
then 
  with_bam_tools_headers=`echo "$with_bam_tools_headers" | sed "s|~|${HOME}|g"`
  echo "provided BAMTools headers directory is $with_bam_tools_headers"
fi 
if test "$with_bam_tools_library" != ""
then
  with_bam_tools_library=`echo "$with_bam_tools_library" | sed "s|~|${HOME}|g"`
  echo "provided BAMTools library directory is $with_bam_tools_library"
fi

# set language to cpp
AC_PROG_CXX
AC_LANG(C++)

# check for MacOS and append the architecture flag
if test `uname` == "Darwin"
then
  CXXFLAGS="$CXXFLAGS -arch x86_64"
  echo "checking whether we're compiling on Darwin/MacOS... yes"
else
  echo "checking whether we're compiling on Darwin/MacOS... no"
fi

# check to see if the user gave us a location for the BAMTools headers
BAM_HEADERS_DIR=""
if test "$with_bam_tools_headers" == "yes"
then
  AC_MSG_ERROR([You need to provide a directory if you're going to use the --with-bam_tools_headers option. For example: with-bam_tools_headers=/some/path])
elif test "$with_bam_tools_headers" != ""
then
  BAM_HEADERS_DIR="$with_bam_tools_headers"
else
  AC_MSG_NOTICE([no BAMTools headers directory specified])
fi

# check to see if the user gave us a location for the BAMTools libraries
BAM_LIBRARIES_DIR=""
if test "$with_bam_tools_library" == "yes"
then
  AC_MSG_ERROR([You need to provide a directory if you're going to use the --with-bam_tools_library option. For example: with-bam_tools_library=/some/path])
elif test "$with_bam_tools_library" != ""
then
  BAM_LIBRARIES_DIR="$with_bam_tools_library"
else
  AC_MSG_NOTICE([no BAMTools library directory specified])
fi

# check for BAMTools but force 'as-needed' linking to be off
# this works on Ubuntu, but not on MacOS
SAVED_LDFLAGS=$LDFLAGS
SAVED_CXXFLAGS=$CXXFLAGS
CXXFLAGS="$CXXFLAGS -Wl,--no-as-needed"
if test "$BAM_LIBRARIES_DIR" != ""
then
  LDFLAGS="$LDFLAGS -L$BAM_LIBRARIES_DIR -lbamtools"
else 
  LDFLAGS="$LDFLAGS -lbamtools"
fi
if test "$BAM_HEADERS_DIR" != ""
then
  CXXFLAGS="$CXXFLAGS -I$BAM_HEADERS_DIR"
fi
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([#include "api/BamAlignment.h"],
    [BamTools::BamAlignment dummy])],
  [ac_cv_lib_bamtools_BamAlignment=yes],
  [])
LDFLAGS=$SAVED_LDFLAGS
CXXFLAGS=$SAVED_CXXFLAGS

# check for BAM tools using the default linking order (this works on 
# MacOS, but not on Ubuntu)  
SAVED_LDFLAGS=$LDFLAGS
SAVED_CXXFLAGS=$CXXFLAGS
if test "$BAM_LIBRARIES_DIR" != ""
then
  LDFLAGS="$LDFLAGS -L$BAM_LIBRARIES_DIR -lbamtools"
else 
  LDFLAGS="$LDFLAGS -lbamtools"
fi
if test "$BAM_HEADERS_DIR" != ""
then
  CXXFLAGS="$CXXFLAGS -I$BAM_HEADERS_DIR"
fi
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([#include "api/BamAlignment.h"],
    [BamTools::BamAlignment dummy])],
  [ac_cv_lib_bamtools_BamAlignment=yes],
  [])
LDFLAGS=$SAVED_LDFLAGS
CXXFLAGS=$SAVED_CXXFLAGS

# check for gsl
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])
ac_cv_lib_gsl_matrix=no
SAVED_LDFLAGS=$LDFLAGS
LDFLAGS="$LDFLAGS -lgsl"
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([#include "gsl/gsl_matrix_double.h"],
    [gsl_matrix_calloc(5,5)])],
  [ac_cv_lib_gsl_matrix=yes],
  [])
LDFLAGS=$SAVED_LDFLAGS
if test "$ac_cv_lib_gsl_matrix" != yes
then 
  echo "checking for linking against gsl matrix... no"
else
  echo "checking for linking against gsl matrix... yes"
fi

GSLERROR="could not find or link to GSL -- Piranha requires GSL and \
cannot be built without it. GSL can be downloaded from \
http://www.gnu.org/software/gsl/ -- please ensure GSL is \
correctly installed on your system."
if test "$ac_cv_lib_m_cos" != yes 
then
  AC_MSG_ERROR([$GSLERROR])
fi
if test "$ac_cv_lib_gsl_gsl_blas_dgemm" != yes
then
  AC_MSG_ERROR([$GSLERROR])
fi 
if test "$ac_cv_lib_gslcblas_cblas_dgemm" != yes
then
  AC_MSG_ERROR([$GSLERROR])
fi
if test "$ac_cv_lib_gsl_matrix" != yes
then
  AC_MSG_ERROR([$GSLERROR])
fi

if test "$ac_cv_lib_bamtools_BamAlignment" != yes
then
  echo "checking for BAMTools... no"
  AC_MSG_WARN([could not find BAMTools, disabling BAM support])
  AC_SUBST(bamsupport, NO_BAM_SUPPORT)
  AC_SUBST(bamtest_exclude_start_token, ["<!--"])
  AC_SUBST(bamtest_exclude_end_token, ["-->"])
else
  echo "checking for BAMTools... yes"
  AC_SUBST(bamsupport, BAM_SUPPORT)
  AC_SUBST(bamlib, -lbamtools)
  # only if the user told us where the BAM library was do we add it to 
  # the places where the executable should look for libraries.
  if test "$BAM_LIBRARIES_DIR" != ""
  then
    AC_SUBST(rpathflags, [-Wl,-rpath,`cd $BAM_LIBRARIES_DIR; pwd`])
    AC_SUBST(bamlibdir, -L$BAM_LIBRARIES_DIR)
  fi
  # this is ok, because we don't add the -I until the makefile, where 
  # it's added to a bunch of things, so this can be empty. Should make -L
  # and -l work like this too... 
  AC_SUBST(bamhead, $BAM_HEADERS_DIR)
fi

AC_OUTPUT(Makefile)
AC_OUTPUT(src/Makefile)
AC_OUTPUT(src/test/Makefile)
AC_OUTPUT(src/progs/Makefile)
AC_OUTPUT(src/common/Makefile)
AC_OUTPUT(src/common/config.hpp)
AC_OUTPUT(src/test/regressionTests.xml)
