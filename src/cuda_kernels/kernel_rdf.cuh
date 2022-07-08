#ifndef INC_KERNEL_RDF_CUH
#define INC_KERNEL_RDF_CUH
#include "../ImageOption.h"
int Cpptraj_GPU_RDF(unsigned long*, int, double, double, const double*, int, const double*, int,
                     ImageOption::Type, const double*, const double*, const double*);
#endif
