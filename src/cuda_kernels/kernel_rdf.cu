#include "kernel_rdf.cuh"
#include "core_kernels.cuh"
#include "../CpptrajStdio.h"
#include "../Gpu.h"
#if defined(__HIP_PLATFORM_HCC__)
#include <hip/hip_runtime.h>
#include "../HipDefinitions.h"
#endif

static inline int calc_nblocks(int ntotal, int nthreadsPerBlock)
{
  int nblocks = ntotal / nthreadsPerBlock;
  if ( (ntotal % nthreadsPerBlock) != 0 )
    nblocks++;
  return nblocks;
}

/** Report any cuda errors. */
static inline int Cuda_check(cudaError_t err, const char* desc) {
  //cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    mprintf("Warning: CUDA Error %s: %s\n", desc, cudaGetErrorString(err));
    mprinterr("Error: CUDA Error %s: %s\n", desc, cudaGetErrorString(err));
    //return 1;
  }
  return 0;
}

/** Calculate distances between pairs of atoms and bin them into a 1D histogram. */
int Cpptraj_GPU_RDF(unsigned long* bins, int nbins, double maximum2, double one_over_spacing,
                     const double* xyz1, int N1,
                     const double* xyz2, int N2,
                     ImageOption::Type imageType,
                     const double* box, const double* ucell, const double* recip)
{
  int* device_rdf;
  Cuda_check(cudaMalloc(((void**)(&device_rdf)), nbins * sizeof(int)), "Allocating rdf bins");
  Cuda_check(cudaMemset( device_rdf, 0, nbins*sizeof(int) ), "Setting rdf bins to 0");

  double* device_xyz1;
  Cuda_check(cudaMalloc(((void**)(&device_xyz1)), N1 * 3 * sizeof(double)), "Allocating xyz1");
  Cuda_check(cudaMemcpy(device_xyz1, xyz1, N1 * 3 * sizeof(double), cudaMemcpyHostToDevice), "Copying xyz1");

  double* device_xyz2;
  Cuda_check(cudaMalloc(((void**)(&device_xyz2)), N2 * 3 * sizeof(double)), "Allocating xyz2");
  Cuda_check(cudaMemcpy(device_xyz2, xyz2, N2 * 3 * sizeof(double), cudaMemcpyHostToDevice), "Copying xyz2");

  double *boxDev;
  double *ucellDev, *recipDev;
  if (imageType == ImageOption::ORTHO) {
    Cuda_check(cudaMalloc(((void**)(&boxDev)), 3 * sizeof(double)), "Allocating box");
    Cuda_check(cudaMemcpy(boxDev,box, 3 * sizeof(double), cudaMemcpyHostToDevice), "Copying box");
  } else if (imageType == ImageOption::NONORTHO) {
    Cuda_check(cudaMalloc(((void**)(&ucellDev)), 9 * sizeof(double)), "Allocating ucell");
    Cuda_check(cudaMalloc(((void**)(&recipDev)), 9 * sizeof(double)), "Allocating frac");
    Cuda_check(cudaMemcpy(ucellDev,ucell, 9 * sizeof(double), cudaMemcpyHostToDevice), "Copying ucell");
    Cuda_check(cudaMemcpy(recipDev,recip, 9 * sizeof(double), cudaMemcpyHostToDevice), "Copying frac");
  }

  // Determine number of blocks
  unsigned int BLOCKDIM = CpptrajGpu::MaxBlockDim_2D();

  dim3 threadsPerBlock(BLOCKDIM, BLOCKDIM);
  dim3 numBlocks(calc_nblocks(N1, threadsPerBlock.x), calc_nblocks(N2, threadsPerBlock.y));
  mprintf("#Atoms = %i, %i; Threads per block = %i, %i;  #Blocks = %i, %i\n",
          N1, N2, threadsPerBlock.x, threadsPerBlock.y, numBlocks.x, numBlocks.y);

  // Launch kernel
  switch (imageType) {
    case ImageOption::NONORTHO:
      kBinDistances_nonOverlap_nonOrtho<<<numBlocks, threadsPerBlock>>>(
        device_rdf, device_xyz1, N1, device_xyz2, N2, recipDev, ucellDev, maximum2, one_over_spacing);
      break;
    default:
      mprinterr("Internal Error: kernel_rdf: Unhandled image type.\n");
      return 1;
  }
  // Error check
  Cuda_check(cudaGetLastError(), "kernel launch");

  // Copy the result back
  int* local_bins = new int[ nbins ];
  cudaMemcpy(local_bins, device_rdf, nbins*sizeof(int), cudaMemcpyDeviceToHost);
  for (int ibin = 0; ibin != nbins; ibin++) {
    //mprintf("DEBUG:\tBin %i = %i (%i)\n", ibin, local_bins[ibin], device_rdf[ibin]);
    //mprintf("DEBUG:\tBin %i = %i\n", ibin, local_bins[ibin]);
    bins[ibin] += local_bins[ibin];
  }
  delete[] local_bins;
  // Free device memory
  cudaFree(device_rdf);
  cudaFree(device_xyz1);
  cudaFree(device_xyz2);
  if (imageType == ImageOption::ORTHO)
    cudaFree(boxDev);
  else if (imageType == ImageOption::NONORTHO) {
    cudaFree(ucellDev);
    cudaFree(recipDev);
  } 
  return 0;
}
