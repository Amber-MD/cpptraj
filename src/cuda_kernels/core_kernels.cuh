#ifndef INC_CORE_KERNELS_CUH
#define INC_CORE_KERNELS_CUH
#if defined(__HIP_PLATFORM_HCC__)
#include <hip/hip_runtime.h>
#endif
// ----- Device kernel definitions ---------------------------------------------
// No imaging
__global__ void kClosestDistsToPt_NoImage(double*,const double *,const double*,double,int,int,int);
__global__ void kClosestDistsToAtoms_NoImage(double*,const double*,const double *,double,int,int,int,int);
// Orthorhombic imaging
__global__ void kClosestDistsToPt_Ortho(double*,const double*,const double*,double,const double*,int,int,int);
__global__ void kClosestDistsToAtoms_Ortho(double*,const double*,const double*,double,const double*,int,int,int,int);
// Non-orthorhombic imaging
__global__ void kClosestDistsToPt_Nonortho(double*,const double*,const double*,double,const double*,const double*,int,int,int);
__global__ void kClosestDistsToAtoms_Nonortho(double*,const double*,const double*,double,const double*,const double*,int,int,int,int);
#endif
