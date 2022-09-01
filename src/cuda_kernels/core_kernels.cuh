#ifndef INC_CORE_KERNELS_CUH
#define INC_CORE_KERNELS_CUH
#include "../Gpu.h"
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
// RDF no imaging
__global__ void kBinDistances_nonOverlap_NoImage(int*, const CpptrajGpu::FpType*, int, const CpptrajGpu::FpType*, int,
                                               CpptrajGpu::FpType, CpptrajGpu::FpType);
// RDF ortho imaging
__global__ void kBinDistances_nonOverlap_Ortho(int*, const CpptrajGpu::FpType*, int, const CpptrajGpu::FpType*, int,
                                               const CpptrajGpu::FpType*, CpptrajGpu::FpType, CpptrajGpu::FpType);
// RDF nonortho imaging
__global__ void kBinDistances_nonOverlap_nonOrtho(int*, const CpptrajGpu::FpType*, int, const CpptrajGpu::FpType*, int,
                                                  const CpptrajGpu::FpType*, const CpptrajGpu::FpType*, CpptrajGpu::FpType, CpptrajGpu::FpType);
#endif
