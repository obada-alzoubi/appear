#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([nan-toolbox], [3.1.4], [alois.schloegl@gmail.com])
AC_CONFIG_SRCDIR([train.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lblas':
AC_CHECK_LIB([blas], [caxpy])
# FIXME: Replace `main' with a function in `-lgomp':
AC_CHECK_LIB([gomp], [omp_get_num_procs])
# FIXME: Replace `main' with a function in `-llibmat':
AC_CHECK_LIB([libmat], [matClose])
# FIXME: Replace `main' with a function in `-llibmex':
AC_CHECK_LIB([libmex], [mexPrintf])
# FIXME: Replace `main' with a function in `-llibmx':
AC_CHECK_LIB([libmx], [mxGetNumberOfFields])
# FIXME: Replace `main' with a function in `-llinear':
AC_CHECK_LIB([linear], [train])
# FIXME: Replace `main' with a function in `-lpthread':
AC_CHECK_LIB([pthread], [pthread_exit])
# FIXME: Replace `main' with a function in `-lsvm':
AC_CHECK_LIB([svm], [svm_train])


# Checks for header files.
AC_CHECK_HEADERS([float.h limits.h locale.h stdint.h stdlib.h string.h sys/param.h])
AC_CHECK_HEADERS([libsvm/svm.h], 
	[AC_DEFINE([HAVE_EXTERNAL_LIBSVM],[1],[libsvm-dev is installed.])
	 AC_SUBST([HAVE_EXTERNAL_LIBSVM],"1")
	 ])
AC_CHECK_HEADERS([linear.h],
	[AC_DEFINE([HAVE_EXTERNAL_LIBLINEAR],[1],[liblinear-dev is installed.])
	 AC_SUBST([HAVE_EXTERNAL_LIBLINEAR],"1")
	 ])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_INT8_T
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor memset mkdir setlocale sqrt strchr strcspn strdup strncasecmp strrchr strtol strtoul])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
