#include "GistCudaSetup.cuh"
#include "GistCudaCalc.cuh"

#if defined(__HIP_PLATFORM_HCC__)
#include <hip/hip_runtime.h>
#include "../HipDefinitions.h"
#endif

/**
 * Allocate memory on the GPU.
 * @param array: The pointer to the array, which will be allocated on the GPU.
 * @param size: An integer giving the size of the array, which will be allocated.
 * @throws: CudaException if a problem occurs.
 */
__host__
void allocateCuda(void **array, int size) {
  // Check if the array is actually free, if not, it will be freed 
  // (fun fact: checking is not necessary, one could also simply free the memory).
  if ((*array) != NULL) {
    cudaFree(*array);
  }
  // If something goes wrong, throw exception
  if (cudaMalloc(array, size) != cudaSuccess) {
    throw CudaException();
  }
}

/**
 * Copy memory from the CPU to the GPU.
 * @param array: The array from which the values shall be copied.
 * @param array_c: The array on the device, to which the values shall be copied.
 * @param size: The size of the stuff which will be copied.
 * @throws: CudaException if something goes wrong.
 */
__host__
void copyMemoryToDevice(void *array, void *array_c, int size) {
  // If something goes wrong, throw exception
  // In this case only copying can go wrong.
  if (cudaMemcpy(array_c, array, size, cudaMemcpyHostToDevice) != cudaSuccess) {
    throw CudaException();
  }
}

/**
 * A simple helper function that copies a lot of stuff to the GPU (as structs).
 * @param charge: An array holding the charges for the different atoms.
 * @param atomtype: An array holding the integers for the atom types of the different atoms.
 * @param solvent: An array of boolean values, holding the information whether a certain atom is solvent or solute.
 * @param atomNumber: The total number of atoms.
 * @param atomProps_c: A pointer to an array on the GPU, which will hold the atom properties.
 * @param ljA: An array holding the lennard-jones parameter A for each atom type pair.
 * @param ljB: An array holding the lennard-jones parameter B for each atom type pair.
 * @param length: The length of the two aforementioned arrays (ljA & ljB).
 * @param lJparams_c: A pointer to an array on the GPU, which will hold the lj parameters.
 * @throws: CudaException if something bad happens.
 */
__host__
void copyMemoryToDeviceStruct(float *charge, int *atomtype, bool *solvent, int *molecule, int atomNumber, void **atomProps_c,
                              float *ljA, float *ljB, int length, void **lJparams_c) {
  // Check if the two arrays are free. Again, this could be removed (but will stay!)
  if ((*atomProps_c) != NULL) {
    cudaFree(*atomProps_c);
  }
  if ((*lJparams_c) != NULL) {
    cudaFree(*lJparams_c);
  }
  // Allocate the necessary memory on the GPU.
  if (cudaMalloc(atomProps_c, atomNumber * sizeof(AtomProperties)) != cudaSuccess) {
    throw CudaException();
  }
  if (cudaMalloc(lJparams_c, length * sizeof(ParamsLJ)) != cudaSuccess) {
    throw CudaException();
  }

  // Create an array for the lennard-jones parameters.
  ParamsLJ *ljp = (ParamsLJ *) malloc (length * sizeof(ParamsLJ));
  // Add the lennard-jones parameters to the array.
  for (int i = 0; i < length; ++i) {
    ljp[i] = ParamsLJ(ljA[i], ljB[i]);
  }

  // Create an array for the atom properties.
  AtomProperties *array = (AtomProperties *)malloc(atomNumber * sizeof(AtomProperties));
  // Add the properties into the array.
  for (int i = 0; i < atomNumber; ++i) {
    array[i] = AtomProperties(charge[i], atomtype[i], solvent[i], molecule[i]);
  }
  // Copy the memory from the host to the device.
  if (cudaMemcpy((*atomProps_c), array, atomNumber * sizeof(AtomProperties), cudaMemcpyHostToDevice) != cudaSuccess) {
    throw CudaException();
  }
  if (cudaMemcpy((*lJparams_c), ljp, length * sizeof(ParamsLJ), cudaMemcpyHostToDevice) != cudaSuccess) {
    throw CudaException();
  }

  // Free the two arrays (so that no memory leak occurs).
  free(ljp);
  free(array);
}

/**
 * Free an array on the CUDA capable device.
 * @param array: The array you want to free.
 */
__host__
void freeCuda(void *array) {
  cudaFree(array);
}


// This is coded C-like, but uses exceptions.
/**
 * This starts the cuda kernel, thus it is actually a quite long function.
 */
__host__
std::vector<std::vector<float> > doActionCudaEnergy(const double *coords, int *NBindex_c, int ntypes, void *parameter, void *molecule_c,
                            int boxinfo, float *recip_o_box, float *ucell, int maxAtoms, float *min_c, float *max_c, int headAtomType, 
                            float neighbourCut2, int *result_o, int *result_n, float *result_w_c, float *result_s_c,
                            int *result_O_c, int *result_N_c, bool doorder) {
  Coordinates *coords_c   = NULL;
  float *recip_b_c  = NULL;
  float *ucell_c    = NULL;
  
  

  float *result_A = (float *) calloc(maxAtoms, sizeof(float));
  float *result_s = (float *) calloc(maxAtoms, sizeof(float));
  // TODO: Fix this, test is actually a quite bad name here!
  Coordinates *coord_array = (Coordinates *) calloc(maxAtoms, sizeof(Coordinates));
  
  // Casting
  AtomProperties *sender = (AtomProperties *) molecule_c;
  ParamsLJ *lennardJonesParams = (ParamsLJ *) parameter;
  
  // Create Boxinfo and Unit cell. This is actually very important for the speed (otherwise
  // there would be LOTS of access to non-local variables).
  BoxInfo boxinf;
  if (boxinfo != 0) {
    boxinf = BoxInfo(recip_o_box, boxinfo);
  }
  UnitCell ucellN;
  if (boxinfo == 2) {
    ucellN = UnitCell(ucell);
  }
  
  // Add the coordinates to the array.
  // TODO: Fix Test here also!
  for (int i = 0; i < maxAtoms; ++i) {
    coord_array[i] = Coordinates(&coords[i * 3]);
  }

  // vectors that will return the necessary information.
  std::vector<std::vector<float> > result;
  std::vector<float> result_esw;
  std::vector<float> result_eww;

  // Allocate space on the GPU
  if (cudaMalloc(&coords_c, maxAtoms * sizeof(Coordinates)) != cudaSuccess) {
    free(result_A); free(result_s); free(coord_array);
    throw CudaException();
  }


  // Copy the data to the GPU
  if (cudaMemcpy(coords_c, coord_array, maxAtoms * sizeof(Coordinates), cudaMemcpyHostToDevice) != cudaSuccess) {
    cudaFree(coords_c); cudaFree(recip_b_c); cudaFree(ucell_c);
    free(result_A); free(result_s); free(coord_array);
    throw CudaException();
  }
  if (cudaMemcpy(result_w_c, result_A, maxAtoms * sizeof(float), cudaMemcpyHostToDevice) != cudaSuccess) {
    cudaFree(coords_c); cudaFree(recip_b_c); cudaFree(ucell_c);
    free(result_A); free(result_s); free(coord_array);
    throw CudaException();
  }
  if (cudaMemcpy(result_s_c, result_s, maxAtoms * sizeof(float), cudaMemcpyHostToDevice) != cudaSuccess) {
    cudaFree(coords_c); cudaFree(recip_b_c); cudaFree(ucell_c);
    free(result_A); free(result_s); free(coord_array);
    throw CudaException();
  }

  // If the doorder calculation is used, it needs to calculate everything differently, so the slow version is used
  // (this is about 10% slower).
  if (doorder) {
    cudaCalcEnergySlow<<< (maxAtoms + SLOW_BLOCKSIZE) / SLOW_BLOCKSIZE, SLOW_BLOCKSIZE >>> (coords_c, NBindex_c, ntypes, lennardJonesParams, sender,
                                                                                            boxinf, ucellN, maxAtoms, result_w_c, result_s_c, min_c, max_c,
                                                                                            headAtomType, neighbourCut2, result_O_c, result_N_c);
  } else {
    // Uses a 2D array, which is nice for memory access.
    dim3 threadsPerBlock(BLOCKSIZE, BLOCKSIZE);
    dim3 numBlocks((maxAtoms + threadsPerBlock.x) / threadsPerBlock.x, (maxAtoms + threadsPerBlock.y) / threadsPerBlock.y);
    // The actual call of the device function
    cudaCalcEnergy<<<numBlocks, threadsPerBlock>>> (coords_c, NBindex_c, ntypes, lennardJonesParams, sender,
                                                                      boxinf, ucellN, maxAtoms, result_w_c, result_s_c, min_c, max_c,
                                                                      headAtomType, neighbourCut2, result_O_c, result_N_c);
    // Check if there was an error.
    cudaError_t cudaError = cudaGetLastError();
    if (cudaError != cudaSuccess) {
      printf("returned %s\n", cudaGetErrorString(cudaError));
    }
  }
  // Return the results of the calculation to the main memory
  if (cudaMemcpy(result_A, result_w_c, maxAtoms * sizeof(float), cudaMemcpyDeviceToHost) != cudaSuccess) {
    cudaFree(coords_c); cudaFree(recip_b_c); cudaFree(ucell_c);
    free(result_A); free(result_s); free(coord_array);
    throw CudaException();
  }  
  

  if (cudaMemcpy(result_s, result_s_c, maxAtoms * sizeof(float), cudaMemcpyDeviceToHost) != cudaSuccess) {
    cudaFree(coords_c); cudaFree(recip_b_c); cudaFree(ucell_c);
    free(result_A); free(result_s); free(coord_array);
    throw CudaException();
  }


  
  if (cudaMemcpy(result_o, result_O_c, maxAtoms * 4 * sizeof(int), cudaMemcpyDeviceToHost) != cudaSuccess) {
    cudaFree(coords_c); cudaFree(recip_b_c); cudaFree(ucell_c);
    free(result_A); free(result_s); free(coord_array);
    throw CudaException();
  }
  
  if (cudaMemcpy(result_n, result_N_c, maxAtoms * sizeof(int), cudaMemcpyDeviceToHost) != cudaSuccess) {
    cudaFree(coords_c); cudaFree(recip_b_c); cudaFree(ucell_c);
    free(result_A); free(result_s); free(coord_array);
    throw CudaException();
  }

  for (int i = 0; i < maxAtoms; ++i) {
    result_eww.push_back(result_A[i]);
    result_esw.push_back(result_s[i]);
  }

  result.push_back(result_eww);
  result.push_back(result_esw);

  // Free everything used in here.
  cudaFree(coords_c); cudaFree(recip_b_c); cudaFree(ucell_c);
  free(result_A); free(result_s); free(coord_array);
  
  return result;
}
