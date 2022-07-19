#ifndef INC_KERNEL_RDF_CUH
#define INC_KERNEL_RDF_CUH
#include "../Gpu.h"
#include "../ImageOption.h"
int Cpptraj_GPU_RDF(unsigned long*, int, CpptrajGpu::FpType, CpptrajGpu::FpType,
                    const CpptrajGpu::FpType*, int, const CpptrajGpu::FpType*, int,
                    ImageOption::Type, const CpptrajGpu::FpType*,
                    const CpptrajGpu::FpType*, const CpptrajGpu::FpType*);
#endif
