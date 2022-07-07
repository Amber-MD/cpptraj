#include "kernel_rdf.cuh"
#include "../CpptrajStdio.h"
#if defined(__HIP_PLATFORM_HCC__)
#include <hip/hip_runtime.h>
#include "../HipDefinitions.h"
#endif

#define BLOCKDIM 32

static inline int calc_nblocks(int ntotal, int nthreadsPerBlock)
{
  int nblocks = ntotal / nthreadsPerBlock;
  if ( (ntotal % nthreadsPerBlock) != 0 )
    nblocks++;
  return nblocks;
}

/** Calculate distances between pairs of atoms and bin them into a 1D histogram. */
void Cpptraj_GPU_RDF(unsigned long* bins,
                     const double* xyz1, int N1,
                     const double* xyz2, int N2,
                     ImageOption::Type imageType,
                     const double* box, const double* ucell, const double* recip)
{
  double* device_xyz1;
  cudaMalloc(((void**)(&device_xyz1)), N1 * 3 * sizeof(double));

  double* device_xyz2;
  cudaMalloc(((void**)(&device_xyz2)), N2 * 3 * sizeof(double));

  double *boxDev;
  double *ucellDev, *recipDev;
  if (imageType == ImageOption::ORTHO) {
    cudaMalloc(((void**)(&boxDev)), 3 * sizeof(double));
    cudaMemcpy(boxDev,box, 3 * sizeof(double), cudaMemcpyHostToDevice);
  } else if (imageType == ImageOption::NONORTHO) {
    cudaMalloc(((void**)(&ucellDev)), 9 * sizeof(double));
    cudaMalloc(((void**)(&recipDev)), 9 * sizeof(double));
    cudaMemcpy(ucellDev,ucell, 9 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(recipDev,recip, 9 * sizeof(double), cudaMemcpyHostToDevice);
  }

  // Determine number of blocks
  dim3 threadsPerBlock(BLOCKDIM, BLOCKDIM);
  dim3 numBlocks(calc_nblocks(N1, threadsPerBlock.x), calc_nblocks(N2, threadsPerBlock.y));
  mprintf("#Atoms = %i, %i; Threads per block = %i, %i;  #Blocks = %i, %i\n",
          N1, N2, threadsPerBlock.x, threadsPerBlock.y, numBlocks.x, numBlocks.y);

}
