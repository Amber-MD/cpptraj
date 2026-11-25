#ifndef INC_CUDA_KERNELS_CUH
#define INC_CUDA_KERNELS_CUH
#include "../ImageOption.h"
#include <cstddef>
// CUDA kernel wrappers
void Action_Closest_Center(const double*,double*,const double*,double,int,int,ImageOption::Type,const double*,const double*,const double*);
// Reusable device buffers for Action_Closest CUDA calculations.
struct ClosestCudaWorkspace {
  double* distances;
  double* solventCoords;
  double* soluteCoords;
  double* box;
  double* ucell;
  double* recip;
  size_t distancesCapacity;
  size_t solventCapacity;
  size_t soluteCapacity;

  ClosestCudaWorkspace();
  ~ClosestCudaWorkspace();
  // Ensure sufficient device capacity; return false on allocation failure.
  bool EnsureCapacity(size_t nMols, size_t nSolventCoords,
                      size_t nSoluteCoords, ImageOption::Type type);

  ClosestCudaWorkspace(const ClosestCudaWorkspace&) = delete;
  ClosestCudaWorkspace& operator=(const ClosestCudaWorkspace&) = delete;
};
bool Action_Closest_NoCenter(const double*, double*, const double*,
                             double, int, int, int, ImageOption::Type,
                             const double*, const double*, const double*,
                             ClosestCudaWorkspace&);
#endif
