#ifndef __OCTROPE_INTERNAL_H
#define __OCTROPE_INTERNAL_H

#if (__cplusplus || c_plusplus)
extern "C" {
#endif 
  
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_STDIO_H
#include <stdio.h>
#endif

#ifdef HAVE_FLOAT_H
#include <float.h>
#endif

#ifndef DBL_MAX
#  ifdef FLT_MAX
#    define DBL_MAX FLT_MAX
#  else
#    define DBL_MAX 3.40282347e+38F
#  endif
#endif

#if (__cplusplus || c_plusplus)
};
#endif
#endif

