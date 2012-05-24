AC_INIT(myconfig, version-0.1)

AC_PROG_CXX
AC_LANG(C++)

# Allow the user to specify where the BAMTools headers and libraries are
# located.
AC_ARG_WITH(stuff, [  --with-bam_tools_headers            specify the location of BAMTools headers if they're not in a default search path])
AC_ARG_WITH(stuff, [  --with-bam_tools_library            specify the location of the BAMTools libraries if they're not in a default search path])

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
  AC_MSG_NOTICE([AC_MSG_ERROR([You need to provide a directory if you're going to use the --with-bam_tools_library option. For example: with-bam_tools_library=/some/path])])
elif test "$with_bam_tools_library" != ""
then
  BAM_LIBRARIES_DIR="$with_bam_tools_library"
else
  AC_MSG_NOTICE([no BAMTools library directory specified])
fi

# check for gsl
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

# check for BAMTools
SAVED_LDFLAGS=$LDFLAGS
SAVED_CXXFLAGS=$CXXFLAGS
if test "$BAM_LIBRARIES_DIR" != ""
then
  LDFLAGS="$LDFLAGS -L$BAM_LIBRARIES_DIR -lbamtools"
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

if test "$ac_cv_lib_gsl_gsl_blas_dgemm" == no
then
  AC_MSG_ERROR([could not find GSL])
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
  AC_SUBST(bamlibdir, -L$BAM_LIBRARIES_DIR)
  AC_SUBST(bamhead, $BAM_HEADERS_DIR)
fi

AC_OUTPUT(src/progs/Makefile)
AC_OUTPUT(src/common/config.hpp)
AC_OUTPUT(src/test/regressionTests.xml)