#ifdef DEBUG_CUDA
#include <cstdio>
#endif
#include "../DistRoutines.h"

#define BLOCKDIM 512

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
// -----------------------------------------------------------------------------

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
                           double maxD, int NMols, int NAtoms, ImagingType type,
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

  if (type == ORTHO)
  {
    cudaMalloc(((void**)(&boxDev)), 3 * sizeof(double));
    cudaMemcpy(boxDev,box, 3 * sizeof(double), cudaMemcpyHostToDevice);
  }
  if (type == NONORTHO)
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
    case NOIMAGE:
      kClosestDistsToPt_NoImage<<<dimGrid0,dimBlock0>>>(devO1Ptr,devI1Ptr, devI2Ptr, maxD, NMols, NAtoms,active_size);
      break;
    case ORTHO:
      kClosestDistsToPt_Ortho<<<dimGrid0,dimBlock0>>>(devO1Ptr,devI1Ptr, devI2Ptr, maxD,boxDev, NMols, NAtoms,active_size);
      break;
    case NONORTHO:
      kClosestDistsToPt_Nonortho<<<dimGrid0,dimBlock0>>>(devO1Ptr,devI1Ptr, devI2Ptr, maxD,ucellDev, recipDev, NMols, NAtoms,active_size);
  }

  cudaThreadSynchronize();

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
  if (type == ORTHO)
    cudaFree(boxDev);
  if (type == NONORTHO)
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
void Action_Closest_NoCenter(const double *SolventMols_, double *D_, const double *Solute_atoms,
                             double maxD, int NMols, int NAtoms, int NSAtoms, ImagingType type,
                             const double* box, const double* ucell, const double* recip)
{
  #ifdef DEBUG_CUDA
  cudaEvent_t start_event, stop_event;
  float time_gpu;
  #endif
  double *devI3Ptr;
  double *devI2Ptr;
  double *devO1Ptr;
  double *boxDev;
  double *ucellDev, *recipDev;

  cudaMalloc(((void **)(&devO1Ptr)),NMols * sizeof(double ));

  cudaMalloc(((void **)(&devI2Ptr)),NMols * NAtoms * 3 * sizeof(double ));
  cudaMemcpy(devI2Ptr,SolventMols_,NMols * NAtoms * 3 * sizeof(double ),cudaMemcpyHostToDevice);
  
  cudaMalloc(((void **)(&devI3Ptr)), NSAtoms * 3 * sizeof(double ));
  cudaMemcpy(devI3Ptr,Solute_atoms,NSAtoms * 3 * sizeof(double ),cudaMemcpyHostToDevice);

  if (type == ORTHO)
  {
    cudaMalloc(((void**)(&boxDev)), 3 * sizeof(double));
    cudaMemcpy(boxDev,box, 3 * sizeof(double), cudaMemcpyHostToDevice);
  }
  if (type == NONORTHO)
  {
    cudaMalloc(((void**)(&ucellDev)), 9 * sizeof(double));
    cudaMalloc(((void**)(&recipDev)), 9 * sizeof(double));
    cudaMemcpy(ucellDev,ucell, 9 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(recipDev,recip, 9 * sizeof(double), cudaMemcpyHostToDevice);
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
    case NOIMAGE:
      kClosestDistsToAtoms_NoImage<<<dimGrid0,dimBlock0>>>(devO1Ptr, devI2Ptr,devI3Ptr, maxD, NMols, NAtoms,NSAtoms,active_size);
      break;
    case ORTHO:
      kClosestDistsToAtoms_Ortho<<<dimGrid0,dimBlock0>>>(devO1Ptr, devI2Ptr,devI3Ptr, maxD, boxDev,  NMols, NAtoms,NSAtoms,active_size);
      break;
    case NONORTHO:
      kClosestDistsToAtoms_Nonortho<<<dimGrid0,dimBlock0>>>(devO1Ptr, devI2Ptr,devI3Ptr, maxD, ucellDev, recipDev,  NMols, NAtoms,NSAtoms,active_size);
    break;
  }
  
  cudaThreadSynchronize();

  #ifdef DEBUG_CUDA
  cudaEventRecord(stop_event, 0);
  cudaEventSynchronize(stop_event);
  cudaEventElapsedTime(&time_gpu,start_event, stop_event );

  printf("Done with kernel CUDA Kernel Time: %.2f\n", time_gpu);
  #endif
  
  cudaMemcpy(D_,devO1Ptr,NMols * sizeof(double ),cudaMemcpyDeviceToHost);
  cudaFree(devO1Ptr);
  cudaFree(devI2Ptr);
  cudaFree(devI3Ptr);
  if (type == ORTHO)
    cudaFree(boxDev);
  if (type == NONORTHO)
  {
    cudaFree(ucellDev);
    cudaFree(recipDev);
  }
}
