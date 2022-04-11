#ifndef INC_CUDA_KERNELS_CUH
#define INC_CUDA_KERNELS_CUH
#include "../ImageOption.h"
// CUDA kernel wrappers
void Action_Closest_Center(const double*,double*,const double*,double,int,int,ImageOption::Type,const double*,const double*,const double*);
void Action_Closest_NoCenter(const double*,double*,const double*,double,int,int,int,ImageOption::Type,const double*,const double*,const double*);
#endif
