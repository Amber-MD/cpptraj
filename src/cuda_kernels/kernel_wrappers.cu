#include "kernel_wrappers.cuh"
#ifdef DEBUG_CUDA
#include <cstdio>
#endif
#if defined(__HIP_PLATFORM_HCC__)
#include <hip/hip_runtime.h>
#include "../HipDefinitions.h"
#endif
#include "core_kernels.cuh"

#define BLOCKDIM 512
// Initialize reusable CUDA device buffers.
ClosestCudaWorkspace::ClosestCudaWorkspace() :
  distances(0),
  solventCoords(0),
  soluteCoords(0),
  box(0),
  ucell(0),
  recip(0),
  distancesCapacity(0),
  solventCapacity(0),
  soluteCapacity(0)
{}

// Release all device buffers owned by this workspace.
ClosestCudaWorkspace::~ClosestCudaWorkspace() {
  if (distances != 0) cudaFree(distances);
  if (solventCoords != 0) cudaFree(solventCoords);
  if (soluteCoords != 0) cudaFree(soluteCoords);
  if (box != 0) cudaFree(box);
  if (ucell != 0) cudaFree(ucell);
  if (recip != 0) cudaFree(recip);
}
// Grow a device buffer only when additional capacity is required.
// Keep the existing allocation if cudaMalloc fails.
static bool EnsureDoubleBuffer(double*& ptr, size_t& capacity,
                               size_t required)
{
  if (required <= capacity)
    return true;

  double* newPtr = 0;
  if (cudaMalloc((void**)&newPtr, required * sizeof(double)) != cudaSuccess)
    return false;

  if (ptr != 0)
    cudaFree(ptr);

  ptr = newPtr;
  capacity = required;
  return true;
}

// Allocate reusable buffers for closest-distance calculations.
bool ClosestCudaWorkspace::EnsureCapacity(
    size_t nMols, size_t nSolventCoords,
    size_t nSoluteCoords, ImageOption::Type type)
{
  if (!EnsureDoubleBuffer(distances, distancesCapacity, nMols))
    return false;

  if (!EnsureDoubleBuffer(solventCoords, solventCapacity, nSolventCoords))
    return false;

  if (!EnsureDoubleBuffer(soluteCoords, soluteCapacity, nSoluteCoords))
    return false;

  if (type == ImageOption::ORTHO && box == 0) {
    if (cudaMalloc((void**)&box, 3 * sizeof(double)) != cudaSuccess)
      return false;
  }

  if (type == ImageOption::NONORTHO) {
    if (ucell == 0) {
      if (cudaMalloc((void**)&ucell, 9 * sizeof(double)) != cudaSuccess)
        return false;
    }
    if (recip == 0) {
      if (cudaMalloc((void**)&recip, 9 * sizeof(double)) != cudaSuccess)
        return false;
    }
  }

  return true;
}

/** Calculate the closest distances between atoms in solvent molecules and
  * the given point.
  * \param SolventMols_ Coordinates for each atom of each solvent molecule.
  * \param D_ Output distances for each molecule.
  * \param maskCenter Point to calculate distances to.
  * \param maxD Maximum possible distance.
  * \param NMols Number of solvent molecules.
  * \param NAtoms Number of atoms in each solvent molecule.
  * \param type Imaging type (none, ortho, non-ortho).
  * \param box Box coordinates.
  * \param ucell Unit cell matrix.
  * \param recip Fractional cell matrix.
  */
void Action_Closest_Center(const double *SolventMols_, double *D_, const double* maskCenter,
                           double maxD, int NMols, int NAtoms, ImageOption::Type type,
                           const double* box, const double* ucell, const double* recip)
{
  #ifdef DEBUG_CUDA
  cudaEvent_t start_event, stop_event;
  float time_gpu;
  #endif
  double *devI2Ptr;
  double *devI1Ptr;
  double *devO1Ptr;
  double *boxDev;
  double *ucellDev, *recipDev;

  cudaMalloc(((void **)(&devO1Ptr)),NMols * sizeof(double ));

  cudaMalloc(((void **)(&devI1Ptr)),3 * sizeof(double ));
  cudaMemcpy(devI1Ptr,maskCenter,3 * sizeof(double ),cudaMemcpyHostToDevice);

  cudaMalloc(((void **)(&devI2Ptr)),NMols * NAtoms * 3 * sizeof(double ));
  cudaMemcpy(devI2Ptr,SolventMols_,NMols * NAtoms * 3 * sizeof(double ),cudaMemcpyHostToDevice);

  if (type == ImageOption::ORTHO)
  {
    cudaMalloc(((void**)(&boxDev)), 3 * sizeof(double));
    cudaMemcpy(boxDev,box, 3 * sizeof(double), cudaMemcpyHostToDevice);
  }
  if (type == ImageOption::NONORTHO)
  {
    cudaMalloc(((void**)(&ucellDev)), 9 * sizeof(double));
    cudaMalloc(((void**)(&recipDev)), 9 * sizeof(double));
    cudaMemcpy(ucellDev,ucell, 9 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(recipDev,recip, 9 * sizeof(double), cudaMemcpyHostToDevice);
  }

  int active_size  =  BLOCKDIM/NAtoms * NAtoms;
  int NBlocks = ceil(float(NMols)/ (BLOCKDIM));

  dim3 dimGrid0 = dim3(NBlocks,1);
  dim3 dimBlock0 = dim3(BLOCKDIM,1);

  #ifdef DEBUG_CUDA
  printf("NMols =  %d, NAtoms = %d\n", NMols, NAtoms);
  printf("active_size =  %d\n", active_size);
  printf("NBlocks =  %d\n", NBlocks);
  printf("sizeof(double) = %d\n", sizeof(double));
  printf("About to launch kernel.\n");

  cudaEventCreate(&start_event);
  cudaEventCreate(&stop_event);
  cudaEventRecord(start_event, 0);
  #endif

  switch (type) {
    case ImageOption::NO_IMAGE:
      kClosestDistsToPt_NoImage<<<dimGrid0,dimBlock0>>>(devO1Ptr,devI1Ptr, devI2Ptr, maxD, NMols, NAtoms,active_size);
      break;
    case ImageOption::ORTHO:
      kClosestDistsToPt_Ortho<<<dimGrid0,dimBlock0>>>(devO1Ptr,devI1Ptr, devI2Ptr, maxD,boxDev, NMols, NAtoms,active_size);
      break;
    case ImageOption::NONORTHO:
      kClosestDistsToPt_Nonortho<<<dimGrid0,dimBlock0>>>(devO1Ptr,devI1Ptr, devI2Ptr, maxD,ucellDev, recipDev, NMols, NAtoms,active_size);
  }

  //cudaThreadSynchronize();

  #ifdef DEBUG_CUDA
  cudaEventRecord(stop_event, 0);
  cudaEventSynchronize(stop_event);
  cudaEventElapsedTime(&time_gpu,start_event, stop_event );

  printf("Done with kernel CUDA Kernel Time: %.2f\n", time_gpu);
  #endif

  cudaMemcpy(D_,devO1Ptr,NMols * sizeof(double ),cudaMemcpyDeviceToHost);
  cudaFree(devO1Ptr);
  cudaFree(devI1Ptr);
  cudaFree(devI2Ptr);
  if (type == ImageOption::ORTHO)
    cudaFree(boxDev);
  if (type == ImageOption::NONORTHO)
  {
    cudaFree(ucellDev);
    cudaFree(recipDev);
  }
}

/** Calculate the closest distances between atoms in solvent molecules and
  * each solute atom.
  * \param SolventMols_ Coordinates for each atom of each solvent molecule.
  * \param D_ Output distances for each molecule.
  * \param Solute_atoms Coordinates for each solute atom.
  * \param maxD Maximum possible distance.
  * \param NMols Number of solvent molecules.
  * \param NAtoms Number of atoms in each solvent molecule.
  * \param NSAtoms Number of solute atoms.
  * \param type Imaging type (none, ortho, non-ortho).
  * \param box Box coordinates.
  * \param ucell Unit cell matrix.
  * \param recip Fractional cell matrix.
  */
bool Action_Closest_NoCenter(const double *SolventMols_, double *D_, const double *Solute_atoms,
                             double maxD, int NMols, int NAtoms, int NSAtoms, ImageOption::Type type,
                             const double* box, const double* ucell, const double* recip,
                             ClosestCudaWorkspace& workspace)
{
  #ifdef DEBUG_CUDA
  cudaEvent_t start_event, stop_event;
  float time_gpu;
  #endif
  // Allocate device buffers only when more capacity is required.
  if (!workspace.EnsureCapacity(NMols, NMols * NAtoms * 3,
                                NSAtoms * 3, type))
    return false;

  double* devO1Ptr = workspace.distances;
  double* devI2Ptr = workspace.solventCoords;
  double* devI3Ptr = workspace.soluteCoords;
  double* boxDev = workspace.box;
  double* ucellDev = workspace.ucell;
  double* recipDev = workspace.recip;

  if (cudaMemcpy(devI2Ptr, SolventMols_,
                 NMols * NAtoms * 3 * sizeof(double),
                 cudaMemcpyHostToDevice) != cudaSuccess)
    return false;

  if (cudaMemcpy(devI3Ptr, Solute_atoms,
                 NSAtoms * 3 * sizeof(double),
                 cudaMemcpyHostToDevice) != cudaSuccess)
    return false;

  if (type == ImageOption::ORTHO) {
    if (cudaMemcpy(boxDev, box, 3 * sizeof(double),
                   cudaMemcpyHostToDevice) != cudaSuccess)
      return false;
  }

  if (type == ImageOption::NONORTHO) {
    if (cudaMemcpy(ucellDev, ucell, 9 * sizeof(double),
                   cudaMemcpyHostToDevice) != cudaSuccess)
      return false;

    if (cudaMemcpy(recipDev, recip, 9 * sizeof(double),
                   cudaMemcpyHostToDevice) != cudaSuccess)
      return false;
  }
  int active_size  =  BLOCKDIM/NAtoms * NAtoms;
  int NBlocks =  ceil(NMols * NAtoms / float(active_size));

  dim3 dimGrid0 = dim3(NBlocks,1);
  dim3 dimBlock0 = dim3(BLOCKDIM,1);

  #ifdef DEBUG_CUDA
  printf("NMols =  %d, NAtoms = %d\n", NMols, NAtoms);
  printf("active_size =  %d\n", active_size);
  printf("NBlocks =  %d\n", NBlocks);
  printf("sizeof(double) = %d\n", sizeof(double));
  printf("About to launch kernel.\n");

  cudaEventCreate(&start_event);
  cudaEventCreate(&stop_event);
  cudaEventRecord(start_event, 0);
  #endif

  switch (type) {
    case ImageOption::NO_IMAGE:
      kClosestDistsToAtoms_NoImage<<<dimGrid0,dimBlock0>>>(devO1Ptr, devI2Ptr,devI3Ptr, maxD, NMols, NAtoms,NSAtoms,active_size);
      break;
    case ImageOption::ORTHO:
      kClosestDistsToAtoms_Ortho<<<dimGrid0,dimBlock0>>>(devO1Ptr, devI2Ptr,devI3Ptr, maxD, boxDev,  NMols, NAtoms,NSAtoms,active_size);
      break;
    case ImageOption::NONORTHO:
      kClosestDistsToAtoms_Nonortho<<<dimGrid0,dimBlock0>>>(devO1Ptr, devI2Ptr,devI3Ptr, maxD, ucellDev, recipDev,  NMols, NAtoms,NSAtoms,active_size);
    break;
  }
  // Detect errors launching the CUDA kernel.
  if (cudaGetLastError() != cudaSuccess)
    return false;
  //cudaThreadSynchronize();

  #ifdef DEBUG_CUDA
  cudaEventRecord(stop_event, 0);
  cudaEventSynchronize(stop_event);
  cudaEventElapsedTime(&time_gpu,start_event, stop_event );

  printf("Done with kernel CUDA Kernel Time: %.2f\n", time_gpu);
  #endif

  // Copy distances back to the host. Device memory belongs to workspace.
  if (cudaMemcpy(D_, devO1Ptr, NMols * sizeof(double),
                 cudaMemcpyDeviceToHost) != cudaSuccess)
    return false;

  return true;
}
