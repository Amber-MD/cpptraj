
#include <cstdio>
#include "../DistRoutines.h"


#define BLOCKDIM 512



// device kernel def
__global__ void Action_noImage_center_GPU(double *D_,double *maskCenter,double *SolventMols_,double maxD, int Nmols , int NAtoms, int active_size);
__global__ void Action_noImage_no_center_GPU(double *D_,double *SolventMols_,double *Solute_atoms ,double maxD, int Nmols , int NAtoms,int NSAtoms , int active_size);


//for imaging with ortho
__global__ void Action_ImageOrtho_center_GPU(double *D_,double *maskCenter,double *SolventMols_,double maxD, double *box, int Nmols , int NAtoms, int active_size);
__global__ void Action_ImageOrtho_no_center_GPU(double *D_,double *SolventMols_,double *Solute_atoms ,double maxD, double *box, int Nmols , int NAtoms,int NSAtoms , int active_size);

//for imaging with NonOrtho
__global__ void Action_ImageNonOrtho_center_GPU(double *D_,double *maskCenter,double *SolventMols_,double maxD, double *ucell, double *recip ,int Nmols , int NAtoms, int active_size);
__global__ void Action_ImageNonOrtho_no_center_GPU(double *D_,double *SolventMols_,double *Solute_atoms ,double maxD, double *ucell, double *recip, int Nmols , int NAtoms,int NSAtoms , int active_size);

////////////////////////





void Action_NoImage_Center(double *SolventMols_,double *D_, double maskCenter[3],double maxD,int  NMols, int NAtoms, float &time_gpu, ImagingType type, double box[3], double ucell[9], double recip[9])
{


  #ifdef DEBUG_CUDA
  cudaEvent_t start_event, stop_event;
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

  if(type == NOIMAGE)
    Action_noImage_center_GPU<<<dimGrid0,dimBlock0>>>(devO1Ptr,devI1Ptr, devI2Ptr, maxD, NMols, NAtoms,active_size);
  else if (type == ORTHO)
    Action_ImageOrtho_center_GPU<<<dimGrid0,dimBlock0>>>(devO1Ptr,devI1Ptr, devI2Ptr, maxD,boxDev, NMols, NAtoms,active_size);
  else if (type == NONORTHO )
    Action_ImageNonOrtho_center_GPU<<<dimGrid0,dimBlock0>>>(devO1Ptr,devI1Ptr, devI2Ptr, maxD,ucellDev, recipDev, NMols, NAtoms,active_size);
  else
    printf("kernel_wrapper: error in Imagingtype\n");

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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Action_NoImage_no_Center(double *SolventMols_,double *D_, double *Solute_atoms,double maxD,int  NMols, int NAtoms,int NSAtoms, float &time_gpu, ImagingType type,double box[3], double ucell[9], double recip[9])
{


  #ifdef DEBUG_CUDA
  cudaEvent_t start_event, stop_event;
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

  if(type == NOIMAGE)
    Action_noImage_no_center_GPU<<<dimGrid0,dimBlock0>>>(devO1Ptr, devI2Ptr,devI3Ptr, maxD, NMols, NAtoms,NSAtoms,active_size);
  else if(type == ORTHO)
    Action_ImageOrtho_no_center_GPU<<<dimGrid0,dimBlock0>>>(devO1Ptr, devI2Ptr,devI3Ptr, maxD, boxDev,  NMols, NAtoms,NSAtoms,active_size);
  else if (type == NONORTHO)
    Action_ImageNonOrtho_no_center_GPU<<<dimGrid0,dimBlock0>>>(devO1Ptr, devI2Ptr,devI3Ptr, maxD, ucellDev, recipDev,  NMols, NAtoms,NSAtoms,active_size);
  else
    printf("kernel_wrapper: error in type no center version\n");
  
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
