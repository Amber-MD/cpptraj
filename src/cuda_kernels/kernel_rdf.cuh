#ifndef INC_KERNEL_RDF_CUH
#define INC_KERNEL_RDF_CUH
#include "../ImageOption.h"
void Cpptraj_GPU_RDF(int*, const double*, int, const double*, int,
                     ImageOption::Type, const double*, const double*, const double*);
#endif
