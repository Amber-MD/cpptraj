#include "core_kernels.cuh"
#include "NonOrtho_dist2.cuh"
#include "ortho_dist2.cuh"
//#include <cstdio> // DEBUG
#define BLOCKDIM 512
#define RSIZE 512

// -----------------------------------------------------------------------------
//try thread coarsening 
/** Calculate the closest distances of atoms of solvent molecules to a point. */
__global__ void kClosestDistsToPt_NoImage(double* D_, const double* maskCenter,
                                          const double* SolventMols_,
                                          double maxD, int Nmols, int NAtoms,
                                          int active_size)
{
  //__shared__ double dist_array[BLOCKDIM];

  //int mol  =  (blockIdx.x * active_size + threadIdx.x)/NAtoms; 
  //int atom  = (blockIdx.x * active_size + threadIdx.x) - (mol * NAtoms);
  //int mol_in_block = threadIdx.x/NAtoms;

  int mol = blockIdx.x*BLOCKDIM + threadIdx.x;

  //advantage of register
  double a0 = maskCenter[0];
  double a1 = maskCenter[1];
  double a2 = maskCenter[2];

  if ( mol < Nmols )
  {
    int sIndex =  mol*NAtoms*3;
    double min_val  = maxD;
    for(int offset  = 0 ; offset < NAtoms*3 ; offset+=3 )
    {
      //double x = a0 - SolventMols_[sIndex++];
      //double y = a1 - SolventMols_[sIndex++];
      //double z = a2 - SolventMols_[sIndex++];

      double x = a0 - SolventMols_[sIndex+ offset + 0 ];
      double y = a1 - SolventMols_[sIndex+offset + 1];
      double z = a2 - SolventMols_[sIndex+offset + 2];

      min_val  =  min(min_val, x*x + y*y + z*z);
    }

    D_[mol] = min_val;
  }
}

// -----------------------------------------------------------------------------
/** Calculate closest distances of atoms of solvent molecules to solute atoms.
  */
__global__ void kClosestDistsToAtoms_NoImage(double* D_,
                                             const double* SolventMols_,
                                             const double* Solute_atoms,
                                             double maxD, int Nmols, int NAtoms,
                                             int NSAtoms, int active_size)
{
  __shared__ double dist_array[BLOCKDIM];
  //__shared__ double sAtom_shared[RSIZE];

  int mol  =  (blockIdx.x * active_size + threadIdx.x)/NAtoms; 
  int atom  = (blockIdx.x * active_size + threadIdx.x) - (mol * NAtoms);
  //int mol_in_block = threadIdx.x/NAtoms;

  //handling the chunks for  solute_atoms
  int chunksize,start,end, NChunks,i,j;

  if(NSAtoms*3 > RSIZE)
  {
    chunksize = (RSIZE/3)*3;
    NChunks = ceil(double(NSAtoms*3)/chunksize);
    start = 0;
    end = chunksize;
  }
  else
  {
    chunksize = NSAtoms*3;
    NChunks = 1;
    start = 0;
    end = NSAtoms*3;
  }

  // if(threadIdx.x == 0 && blockIdx.x == 0 )
  //   printf("chunkszize = %d ; Nchunk =  %d; start = %d; end = %d\n ",
  //     chunksize,NChunks,start,end);

  if ( threadIdx.x < active_size && mol*NAtoms + atom < Nmols*NAtoms )
  {
    // if(atom == 0 )
    //   D_[mol] = maxD;
    //__syncthreads(); 
    double min_val  = maxD;
    double dist;
    int sIndex =  mol*NAtoms*3 + atom*3;
    double a0 = SolventMols_[sIndex + 0];
    double a1 = SolventMols_[sIndex + 1];
    double a2 = SolventMols_[sIndex + 2];

    //this is to imporve cache hits! (in the old days this would be thrown in shared mem)
    for(i  = 0 ; i  < NChunks ; i++)
    {
      //copying to shared
      //if (threadIdx.x < (end - start))
      //  sAtom_shared[threadIdx.x] = Solute_atoms[start + threadIdx.x];

      //__syncthreads();

      //TODO - add skew per thread 
      for (j = start ; j < end; j+=3 )
      {
        //int offset = start + (j + threadIdx.x)%(end - start);
        double x = Solute_atoms[j + 0]  - a0;
        double y = Solute_atoms[j + 1]  - a1;
        double z = Solute_atoms[j + 2]  - a2;
        dist =  x*x + y*y + z*z;
        //if (mol ==  11)
        //  printf("min  = %f\n",min_val);
        min_val = min(min_val,dist);


      }

      start = end;
      end = min(end + chunksize, NSAtoms*3);
    }

    dist_array[threadIdx.x] = min_val;
    //if (threadIdx.x == 0)
    //  printf("min_val  = %f\n",min_val);
    //printf(" dist  =  %f\n", Dist);

    __syncthreads();

    //first thread
    //naive approach to a reduction algorithm
    //this works if NAtoms is small other wise you need split
    //and do some of log(n) parallel reduction 
    //min_val  = maxD;
    if( atom ==0 )
    {
      for(i  = 0 ; i < NAtoms ; i++ ){
        //sIndex = mol*NAtoms*3 + i*3;
        //if (dist_array[threadIdx.x + i]  < min_val) 
        //  min_val = dist_array[threadIdx.x + i] ;
        min_val =  min(min_val, dist_array[threadIdx.x + i]);
      }
      D_[mol] = min_val;
    }
  //if(tx == 0 && bx == 0 )
  //  printf("end of kernel");
  }
}

// -----------------------------------------------------------------------------
/** Calculate the closest distances of atoms of solvent molecules to a point.
  * Perform orthorhombic imaging.
  */
__global__ void kClosestDistsToPt_Ortho(double *D_, const double* maskCenter,
                                        const double* SolventMols_, double maxD,
                                        const double *box, int Nmols,
                                        int NAtoms, int active_size)
{
  //__shared__ double dist_array[BLOCKDIM];

  //int mol  =  (blockIdx.x * active_size + threadIdx.x)/NAtoms; 
  //int atom  = (blockIdx.x * active_size + threadIdx.x) - (mol * NAtoms);
  //int mol_in_block = threadIdx.x/NAtoms;

  int mol = blockIdx.x*BLOCKDIM + threadIdx.x;

  //advantage of register
  double a0 = maskCenter[0];
  double a1 = maskCenter[1];
  double a2 = maskCenter[2];

  if ( mol < Nmols )
  {
    int sIndex =  mol*NAtoms*3;
    double min_val  = maxD;
    double dist;
    for(int offset  = 0 ; offset < NAtoms ; offset++ )
    {
      dist = ortho_dist2<double>( a0, a1, a2,
                                  SolventMols_[sIndex], SolventMols_[sIndex+1], SolventMols_[sIndex+2],
                                  box );
      sIndex += 3;

      if (box[0]==0.0 || box[1]==0.0 || box[2]==0.0)
        dist= -1.0;

      min_val  =  min(min_val, dist);
    }

    D_[mol] = min_val;
  }
}

//------------------------------------------------------------------------------
/** Calculate closest distances of atoms of solvent molecules to solute atoms.
  * Perform orthorhombic imaging.
  */
__global__ void kClosestDistsToAtoms_Ortho(double* D_, const double* SolventMols_,
                                           const double* Solute_atoms, double maxD,
                                           const double* box, int Nmols, int NAtoms,
                                           int NSAtoms, int active_size)
{
  __shared__ double dist_array[BLOCKDIM];
  //__shared__ double sAtom_shared[RSIZE];

  int mol  =  (blockIdx.x * active_size + threadIdx.x)/NAtoms; 
  int atom  = (blockIdx.x * active_size + threadIdx.x) - (mol * NAtoms);
  //int mol_in_block = threadIdx.x/NAtoms;

  //handling the chunks for  solute_atoms
  int chunksize,start,end, NChunks,i,j;

  if(NSAtoms*3 > RSIZE)
  {
    chunksize = (RSIZE/3)*3;
    NChunks = ceil(double(NSAtoms*3)/chunksize);
    start = 0;
    end = chunksize;
  }
  else
  {
    chunksize = NSAtoms*3;
    NChunks = 1;
    start = 0;
    end = NSAtoms*3;
  }

  // if(threadIdx.x == 0 && blockIdx.x == 0 )
  //   printf("chunkszize = %d ; Nchunk =  %d; start = %d; end = %d\n ",
  //     chunksize,NChunks,start,end);

  if ( threadIdx.x < active_size && mol*NAtoms + atom < Nmols*NAtoms )
  {
    // if(atom == 0 )
    //   D_[mol] = maxD;
    //__syncthreads(); 
    double min_val  = maxD;
    double dist;
    int sIndex =  mol*NAtoms*3 + atom*3;
    double a0 = SolventMols_[sIndex + 0];
    double a1 = SolventMols_[sIndex + 1];
    double a2 = SolventMols_[sIndex + 2];

    for(i  = 0 ; i  < NChunks ; i++)
    {
      //copying to shared
      //if (threadIdx.x < (end - start))
      //  sAtom_shared[threadIdx.x] = Solute_atoms[start + threadIdx.x];

      //__syncthreads();

      //TODO - add skew per thread 
      for (j = start ; j < end; j+=3 )
      {
        //int offset = start + (j + threadIdx.x)%(end - start);
        dist = ortho_dist2<double>( Solute_atoms[j], Solute_atoms[j+1], Solute_atoms[j+2],
                                    a0, a1, a2,
                                    box );

        if (box[0]==0.0 || box[1]==0.0 || box[2]==0.0)
          dist = -1.0;

        //if (mol ==  11)
        //  printf("min  = %f\n",min_val);
        min_val = min(min_val,dist);
      }

      start = end;
      end = min(end + chunksize, NSAtoms*3);
    }

    dist_array[threadIdx.x] = min_val;
    //if (threadIdx.x == 0)
    //  printf("min_val  = %f\n",min_val);
    //printf(" dist  =  %f\n", Dist);

    __syncthreads();

    //first thread
    //naive approach to a reduction algorithm
    //this works if NAtoms is small other wise you need split
    //and do some of log(n) parallel reduction 
    //min_val  = maxD;
    if( atom ==0 )
    {
      for(i  = 0 ; i < NAtoms ; i++ ){
        //sIndex = mol*NAtoms*3 + i*3;
        //if (dist_array[threadIdx.x + i]  < min_val) 
        //  min_val = dist_array[threadIdx.x + i] ;
        min_val =  min(min_val, dist_array[threadIdx.x + i]);
      }
      D_[mol] = min_val;
    }

  //if(tx == 0 && bx == 0 )
  //  printf("end of kernel");
  }
}

// -----------------------------------------------------------------------------
/** Calculate the closest distances of atoms of solvent molecules to a point.
  * Perform non-orthorhombic imaging.
  */
__global__ void kClosestDistsToPt_Nonortho(double* D_, const double* maskCenter,
                                           const double* SolventMols_,
                                           double maxD, const double* ucell,
                                           const double *recip, int Nmols,
                                           int NAtoms, int active_size)
{
  //__shared__ double dist_array[BLOCKDIM];

  //int mol  =  (blockIdx.x * active_size + threadIdx.x)/NAtoms; 
  //int atom  = (blockIdx.x * active_size + threadIdx.x) - (mol * NAtoms);
  //int mol_in_block = threadIdx.x/NAtoms;

  int mol = blockIdx.x*BLOCKDIM + threadIdx.x;

  //advantage of register
  double a0 = recip[0]*maskCenter[0] + recip[1]*maskCenter[1] + recip[2]*maskCenter[2];
  double a1 = recip[3]*maskCenter[0] + recip[4]*maskCenter[1] + recip[5]*maskCenter[2];
  double a2 = recip[6]*maskCenter[0] + recip[7]*maskCenter[1] + recip[8]*maskCenter[2];

  if ( mol < Nmols )
  {
    int sIndex =  mol*NAtoms*3;
    double min_val  = maxD;
    for(int offset  = 0 ; offset < NAtoms*3 ; offset+=3 )
    {
      double x =  recip[0]*SolventMols_[sIndex + offset + 0] + recip[1]*SolventMols_[sIndex + offset + 1] + recip[2]*SolventMols_[sIndex + offset + 2];
      double y =  recip[3]*SolventMols_[sIndex + offset + 0] + recip[4]*SolventMols_[sIndex + offset + 1] + recip[5]*SolventMols_[sIndex + offset + 2];
      double z =  recip[6]*SolventMols_[sIndex + offset + 0] + recip[7]*SolventMols_[sIndex + offset + 1] + recip[8]*SolventMols_[sIndex + offset + 2];
      double dist  = NonOrtho_dist2<double>(x,y,z,a0,a1,a2,ucell);
      // if (mol ==  0)
      //   printf("dist  = %f\n",dist);

      min_val  =  min(min_val, dist);
    }

    D_[mol] = min_val;
  }
}

// -----------------------------------------------------------------------------
/** Calculate closest distances of atoms of solvent molecules to solute atoms.
  * Perform non-orthorhombic imaging.
  */
__global__ void kClosestDistsToAtoms_Nonortho(double*D_,
                                              const double* SolventMols_,
                                              const double* Solute_atoms,
                                              double maxD, const double *ucell,
                                              const double* recip, int Nmols,
                                              int NAtoms, int NSAtoms,
                                              int active_size)
{
  __shared__ double dist_array[BLOCKDIM];
  //__shared__ double sAtom_shared[RSIZE];

  int mol  =  (blockIdx.x * active_size + threadIdx.x)/NAtoms; 
  int atom  = (blockIdx.x * active_size + threadIdx.x) - (mol * NAtoms);
  //int mol_in_block = threadIdx.x/NAtoms;

  //handling the chunks for  solute_atoms
  int chunksize,start,end, NChunks,i,j;

  if(NSAtoms*3 > RSIZE)
  {
    chunksize = (RSIZE/3)*3;
    NChunks = ceil(double(NSAtoms*3)/chunksize);
    start = 0;
    end = chunksize;
  }
  else
  {
    chunksize = NSAtoms*3;
    NChunks = 1;
    start = 0;
    end = NSAtoms*3;
  }

  // if(threadIdx.x == 0 && blockIdx.x == 0 )
  //   printf("chunkszize = %d ; Nchunk =  %d; start = %d; end = %d\n ",
  //     chunksize,NChunks,start,end);

  if ( threadIdx.x < active_size && mol*NAtoms + atom < Nmols*NAtoms )
  {
    // if(atom == 0 )
    //   D_[mol] = maxD;
    //__syncthreads(); 
    double min_val  = maxD;
    double dist;
    int sIndex =  mol*NAtoms*3 + atom*3;

    double a0 = recip[0]*SolventMols_[sIndex + 0] + recip[1]*SolventMols_[sIndex + 1] + recip[2]*SolventMols_[sIndex + 2];
    double a1 = recip[3]*SolventMols_[sIndex + 0] + recip[4]*SolventMols_[sIndex + 1] + recip[5]*SolventMols_[sIndex + 2];
    double a2 = recip[6]*SolventMols_[sIndex + 0] + recip[7]*SolventMols_[sIndex + 1] + recip[8]*SolventMols_[sIndex + 2];

    for(i  = 0 ; i  < NChunks ; i++)
    {
      //copying to shared
      //if (threadIdx.x < (end - start))
      //  sAtom_shared[threadIdx.x] = Solute_atoms[start + threadIdx.x];

      //__syncthreads();

      //TODO - add skew per thread 
      for (j = start ; j < end; j+=3 )
      {
        //int offset = start + (j + threadIdx.x)%(end - start);

        double x = recip[0]*Solute_atoms[j + 0]  + recip[1]*Solute_atoms[j + 1]  + recip[2]*Solute_atoms[j + 2] ;
        double y = recip[3]*Solute_atoms[j + 0]  + recip[4]*Solute_atoms[j + 1]  + recip[5]*Solute_atoms[j + 2] ;
        double z = recip[6]*Solute_atoms[j + 0]  + recip[7]*Solute_atoms[j + 1]  + recip[8]*Solute_atoms[j + 2] ;

        dist =  NonOrtho_dist2<double>(x,y,z,a0,a1,a2,ucell);
        //if (mol ==  11)
        //  printf("min  = %f\n",min_val);
        min_val = min(min_val,dist);
      }

      start = end;
      end = min(end + chunksize, NSAtoms*3);
    }

    dist_array[threadIdx.x] = min_val;
    //if (threadIdx.x == 0)
    //  printf("min_val  = %f\n",min_val);
    //printf(" dist  =  %f\n", Dist);

    __syncthreads();

    //first thread
    //naive approach to a reduction algorithm
    //this works if NAtoms is small other wise you need split
    //and do some of log(n) parallel reduction 
    //min_val  = maxD;
    if( atom ==0 )
    {
      for(i  = 0 ; i < NAtoms ; i++ ){
        //sIndex = mol*NAtoms*3 + i*3;
        //if (dist_array[threadIdx.x + i]  < min_val) 
        //  min_val = dist_array[threadIdx.x + i] ;
        min_val =  min(min_val, dist_array[threadIdx.x + i]);
      }
      D_[mol] = min_val;
    }

  //if(tx == 0 && bx == 0 )
  //  printf("end of kernel");
  }
}

// -----------------------------------------------------------------------------
/** Bin distances from two non-overlapping sets of coords. */
__global__ void kBinDistances_nonOverlap_NoImage(int* RDF,
                                               const double* xyz1, int N1, const double* xyz2, int N2,
                                               double maximum2, double one_over_spacing)
{
  int a1 = blockIdx.x * blockDim.x + threadIdx.x;
  int a2 = blockIdx.y * blockDim.y + threadIdx.y;

  if (a1 < N1 && a2 < N2) {
    int idx1 = a1 * 3;
    double a1x = xyz1[idx1  ];
    double a1y = xyz1[idx1+1];
    double a1z = xyz1[idx1+2];

    int idx2 = a2 * 3;
    double x = a1x - xyz2[idx2  ];
    double y = a1y - xyz2[idx2+1];
    double z = a1z - xyz2[idx2+2];

    double dist2 = (x*x) + (y*y) + (z*z); 
    if (dist2 > 0 && dist2 <= maximum2) {
      double dist = sqrt(dist2);
      int histIdx = (int) (dist * one_over_spacing);
      //printf("DEBUG: a1= %i  a2= %i  dist= %f  bin=%i\n", a1+1, a2+1, dist, histIdx);
      //printf("DEBUG: xyz1= %f %f %f\n", a1x, a1y, a1z);
      //printf("DEBUG: a1= %i  a2= %i  dist= %f  bin=%i  xyz1=%f %f %f  xyz2=%f %f %f\n", a1+1, a2+1, dist, histIdx,
      //       a1x, a1y, a1z, a2x, a2y, a2z);
      atomicAdd( RDF + histIdx, 1 );
    }
  }
}

/** Bin distances from two non-overlapping sets of coords. */
__global__ void kBinDistances_nonOverlap_Ortho(int* RDF,
                                               const double* xyz1, int N1, const double* xyz2, int N2,
                                               const double* box,
                                               double maximum2, double one_over_spacing)
{
  int a1 = blockIdx.x * blockDim.x + threadIdx.x;
  int a2 = blockIdx.y * blockDim.y + threadIdx.y;

  if (a1 < N1 && a2 < N2) {
    int idx1 = a1 * 3;
    double a1x = xyz1[idx1  ];
    double a1y = xyz1[idx1+1];
    double a1z = xyz1[idx1+2];

    int idx2 = a2 * 3;
    double a2x = xyz2[idx2  ];
    double a2y = xyz2[idx2+1];
    double a2z = xyz2[idx2+2];

    double dist2 = ortho_dist2<double>(a1x, a1y, a1z, a2x, a2y, a2z, box);
    if (dist2 > 0 && dist2 <= maximum2) {
      double dist = sqrt(dist2);
      int histIdx = (int) (dist * one_over_spacing);
      //printf("DEBUG: a1= %i  a2= %i  dist= %f  bin=%i\n", a1+1, a2+1, dist, histIdx);
      //printf("DEBUG: xyz1= %f %f %f\n", a1x, a1y, a1z);
      //printf("DEBUG: a1= %i  a2= %i  dist= %f  bin=%i  xyz1=%f %f %f  xyz2=%f %f %f\n", a1+1, a2+1, dist, histIdx,
      //       a1x, a1y, a1z, a2x, a2y, a2z);
      atomicAdd( RDF + histIdx, 1 );
    }
  }
}

/** Bin distances from two non-overlapping sets of coords. */
__global__ void kBinDistances_nonOverlap_nonOrtho(int* RDF,
                                                  const double* xyz1, int N1, const double* xyz2, int N2,
                                                  const double* frac, const double* ucell,
                                                  double maximum2, double one_over_spacing)
{
  int a1 = blockIdx.x * blockDim.x + threadIdx.x;
  int a2 = blockIdx.y * blockDim.y + threadIdx.y;

  if (a1 < N1 && a2 < N2) {
    int idx1 = a1 * 3;
    double a1x = xyz1[idx1  ];
    double a1y = xyz1[idx1+1];
    double a1z = xyz1[idx1+2];
    double f1x = frac[0]*a1x + frac[1]*a1y + frac[2]*a1z;
    double f1y = frac[3]*a1x + frac[4]*a1y + frac[5]*a1z;
    double f1z = frac[6]*a1x + frac[7]*a1y + frac[8]*a1z;

    int idx2 = a2 * 3;
    double a2x = xyz2[idx2  ];
    double a2y = xyz2[idx2+1];
    double a2z = xyz2[idx2+2];
    double f2x = frac[0]*a2x + frac[1]*a2y + frac[2]*a2z;
    double f2y = frac[3]*a2x + frac[4]*a2y + frac[5]*a2z;
    double f2z = frac[6]*a2x + frac[7]*a2y + frac[8]*a2z;

    double dist2 =  NonOrtho_dist2<double>(f2x, f2y, f2z, f1x ,f1y, f1z, ucell);
    if (dist2 > 0 && dist2 <= maximum2) {
      double dist = sqrt(dist2);
      int histIdx = (int) (dist * one_over_spacing);
      //printf("DEBUG: a1= %i  a2= %i  dist= %f  bin=%i\n", a1+1, a2+1, dist, histIdx);
      //printf("DEBUG: xyz1= %f %f %f\n", a1x, a1y, a1z);
      //printf("DEBUG: a1= %i  a2= %i  dist= %f  bin=%i  xyz1=%f %f %f  xyz2=%f %f %f\n", a1+1, a2+1, dist, histIdx,
      //       a1x, a1y, a1z, a2x, a2y, a2z);
      atomicAdd( RDF + histIdx, 1 );
    }
  }
}
