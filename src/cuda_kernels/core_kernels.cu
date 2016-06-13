#define BLOCKDIM 512
#define RSIZE 512

// Forward declaration for non-orthorhombic distance calc.
__device__ double NonOrtho_dist(double a0, double a1, double a2,
        double b0, double b1, double b2,
        const double *ucell);

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
      double x = a0 - SolventMols_[sIndex++];
      double y = a1 - SolventMols_[sIndex++];
      double z = a2 - SolventMols_[sIndex++];

      // Get rid of sign info
      if (x<0) x=-x;
      if (y<0) y=-y;
      if (z<0) z=-z;
      // Get rid of multiples of box lengths 
      //TODO  WIERD that should be a way to simplify it
      while (x > box[0]) x = x - box[0];
      while (y > box[1]) y = y - box[1];
      while (z > box[2]) z = z - box[2];
        // Find shortest distance in periodic reference
      double D = box[0] - x;
      if (D < x) x = D;
      D = box[1] - y;
      if (D < y) y = D;  
      D = box[2] - z;
      if (D < z) z = D;

      //Dist = x*x + y*y + z*z;
      dist = x*x + y*y + z*z;
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
        double x = Solute_atoms[j + 0]  - a0;
        double y = Solute_atoms[j + 1]  - a1;
        double z = Solute_atoms[j + 2]  - a2;

        // Get rid of sign info
        if (x<0) x=-x;
        if (y<0) y=-y;
        if (z<0) z=-z;
        // Get rid of multiples of box lengths 
        //TODO  WIERD that should be a way to simplify it
        while (x > box[0]) x = x - box[0];
        while (y > box[1]) y = y - box[1];
        while (z > box[2]) z = z - box[2];

        //below is actually slower! 
        //x = x - box[0]*((int)x/box[0]);
        //y = y - box[0]*((int)y/box[1]);
        //z = z - box[0]*((int)z/box[2]);
        // Find shortest distance in periodic reference
        double D = box[0] - x;
        if (D < x) x = D;
        D = box[1] - y;
        if (D < y) y = D;  
        D = box[2] - z;
        if (D < z) z = D;

        //Dist = x*x + y*y + z*z;
        dist =  x*x + y*y + z*z;
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
      double dist  = NonOrtho_dist(x,y,z,a0,a1,a2,ucell);
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

        dist =  NonOrtho_dist(x,y,z,a0,a1,a2,ucell);
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
/** \return Shortest imaged distance between given coordinates in fractional space.
  * NOTE: This function is complicated hence we will put into a __device__ only function.
  */
__device__ double NonOrtho_dist(double a0, double a1, double a2,
                                double b0, double b1, double b2,
                                const double *ucell)
{
  int ixyz[3];
  double minIn  = -1.0;

   //double closest2
  // The floor() calls serve to bring each point back in the main unit cell.
  double fx = a0 - floor(a0);
  double fy = a1 - floor(a1);
  double fz = a2 - floor(a2); 
  double f2x = b0 - floor(b0);
  double f2y = b1 - floor(b1);
  double f2z = b2 - floor(b2);
  // f2 back in Cartesian space
  double X_factor = (f2x*ucell[0] + f2y*ucell[3] + f2z*ucell[6]);
  double Y_factor = (f2x*ucell[1] + f2y*ucell[4] + f2z*ucell[7]);
  double Z_factor = (f2x*ucell[2] + f2y*ucell[5] + f2z*ucell[8]);
  // Precompute some factors
  double fxm1 = fx - 1.0;
  double fxp1 = fx + 1.0;
  double fym1 = fy - 1.0;
  double fyp1 = fy + 1.0;
  double fzm1 = fz - 1.0;
  double fzp1 = fz + 1.0;

  double fxm1u0 = fxm1 * ucell[0];
  double fxu0   = fx   * ucell[0];
  double fxp1u0 = fxp1 * ucell[0];
  double fxm1u1 = fxm1 * ucell[1];
  double fxu1   = fx   * ucell[1];
  double fxp1u1 = fxp1 * ucell[1];
  double fxm1u2 = fxm1 * ucell[2];
  double fxu2   = fx   * ucell[2];
  double fxp1u2 = fxp1 * ucell[2];

  double fym1u3 = fym1 * ucell[3];
  double fyu3   = fy   * ucell[3];
  double fyp1u3 = fyp1 * ucell[3];
  double fym1u4 = fym1 * ucell[4];
  double fyu4   = fy   * ucell[4];
  double fyp1u4 = fyp1 * ucell[4];
  double fym1u5 = fym1 * ucell[5];
  double fyu5   = fy   * ucell[5];
  double fyp1u5 = fyp1 * ucell[5];

  double fzm1u6 = fzm1 * ucell[6];
  double fzu6   = fz   * ucell[6];
  double fzp1u6 = fzp1 * ucell[6];
  double fzm1u7 = fzm1 * ucell[7];
  double fzu7   = fz   * ucell[7];
  double fzp1u7 = fzp1 * ucell[7];
  double fzm1u8 = fzm1 * ucell[8];
  double fzu8   = fz   * ucell[8];
  double fzp1u8 = fzp1 * ucell[8];

  // Calc ix iy iz = 0 case
  double x = (fxu0 + fyu3 + fzu6) - X_factor;
  double y = (fxu1 + fyu4 + fzu7) - Y_factor;
  double z = (fxu2 + fyu5 + fzu8) - Z_factor;
  // DEBUG
  //mprintf("DEBUG: a2: %g %g %g\n",(fxu0 + fyu3 + fzu6), (fxu1 + fyu4 + fzu7), (fxu2 + fyu5 + fzu8));
  //mprintf("DEBUG: a1: %g %g %g\n", X_factor, Y_factor, Z_factor);
  double min = (x*x) + (y*y) + (z*z);

  if (minIn > 0.0 && minIn < min) min = minIn;

  ixyz[0] = 0;
  ixyz[1] = 0;
  ixyz[2] = 0;

  // -1 -1 -1
  x = (fxm1u0 + fym1u3 + fzm1u6) - X_factor;
  y = (fxm1u1 + fym1u4 + fzm1u7) - Y_factor;
  z = (fxm1u2 + fym1u5 + fzm1u8) - Z_factor;
  double D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] = -1; ixyz[2] = -1; }
  // -1 -1  0
  x = (fxm1u0 + fym1u3 + fzu6  ) - X_factor;
  y = (fxm1u1 + fym1u4 + fzu7  ) - Y_factor;
  z = (fxm1u2 + fym1u5 + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] = -1; ixyz[2] =  0; }
  // -1 -1 +1
  x = (fxm1u0 + fym1u3 + fzp1u6) - X_factor;
  y = (fxm1u1 + fym1u4 + fzp1u7) - Y_factor;
  z = (fxm1u2 + fym1u5 + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] = -1; ixyz[2] =  1; }
  // -1  0 -1
  x = (fxm1u0 + fyu3   + fzm1u6) - X_factor;
  y = (fxm1u1 + fyu4   + fzm1u7) - Y_factor;
  z = (fxm1u2 + fyu5   + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] =  0; ixyz[2] = -1; }
  // -1  0  0
  x = (fxm1u0 + fyu3   + fzu6  ) - X_factor;
  y = (fxm1u1 + fyu4   + fzu7  ) - Y_factor;
  z = (fxm1u2 + fyu5   + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] =  0; ixyz[2] =  0; }
  // -1  0 +1
  x = (fxm1u0 + fyu3   + fzp1u6) - X_factor;
  y = (fxm1u1 + fyu4   + fzp1u7) - Y_factor;
  z = (fxm1u2 + fyu5   + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] =  0; ixyz[2] =  1; }
  // -1 +1 -1
  x = (fxm1u0 + fyp1u3 + fzm1u6) - X_factor;
  y = (fxm1u1 + fyp1u4 + fzm1u7) - Y_factor;
  z = (fxm1u2 + fyp1u5 + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] =  1; ixyz[2] = -1; }
  // -1 +1  0
  x = (fxm1u0 + fyp1u3 + fzu6  ) - X_factor;
  y = (fxm1u1 + fyp1u4 + fzu7  ) - Y_factor;
  z = (fxm1u2 + fyp1u5 + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] =  1; ixyz[2] =  0; }
  // -1 +1 +1
  x = (fxm1u0 + fyp1u3 + fzp1u6) - X_factor;
  y = (fxm1u1 + fyp1u4 + fzp1u7) - Y_factor;
  z = (fxm1u2 + fyp1u5 + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] =  1; ixyz[2] =  1; }

  //  0 -1 -1
  x = (fxu0   + fym1u3 + fzm1u6) - X_factor;
  y = (fxu1   + fym1u4 + fzm1u7) - Y_factor;
  z = (fxu2   + fym1u5 + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] = -1; ixyz[2] = -1; }
  //  0 -1  0
  x = (fxu0   + fym1u3 + fzu6  ) - X_factor;
  y = (fxu1   + fym1u4 + fzu7  ) - Y_factor;
  z = (fxu2   + fym1u5 + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] = -1; ixyz[2] =  0; }
  //  0 -1 +1
  x = (fxu0   + fym1u3 + fzp1u6) - X_factor;
  y = (fxu1   + fym1u4 + fzp1u7) - Y_factor;
  z = (fxu2   + fym1u5 + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] = -1; ixyz[2] =  1; }
  //  0  0 -1
  x = (fxu0   + fyu3   + fzm1u6) - X_factor;
  y = (fxu1   + fyu4   + fzm1u7) - Y_factor;
  z = (fxu2   + fyu5   + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] =  0; ixyz[2] = -1; }
  //  0  0  0
  //  0  0 +1
  x = (fxu0   + fyu3   + fzp1u6) - X_factor;
  y = (fxu1   + fyu4   + fzp1u7) - Y_factor;
  z = (fxu2   + fyu5   + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] =  0; ixyz[2] =  1; }
  //  0 +1 -1
  x = (fxu0   + fyp1u3 + fzm1u6) - X_factor;
  y = (fxu1   + fyp1u4 + fzm1u7) - Y_factor;
  z = (fxu2   + fyp1u5 + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] =  1; ixyz[2] = -1; }
  //  0 +1  0
  x = (fxu0   + fyp1u3 + fzu6  ) - X_factor;
  y = (fxu1   + fyp1u4 + fzu7  ) - Y_factor;
  z = (fxu2   + fyp1u5 + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] =  1; ixyz[2] =  0; }
  //  0 +1 +1
  x = (fxu0   + fyp1u3 + fzp1u6) - X_factor;
  y = (fxu1   + fyp1u4 + fzp1u7) - Y_factor;
  z = (fxu2   + fyp1u5 + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] =  1; ixyz[2] =  1; }

  // +1 -1 -1
  x = (fxp1u0 + fym1u3 + fzm1u6) - X_factor;
  y = (fxp1u1 + fym1u4 + fzm1u7) - Y_factor;
  z = (fxp1u2 + fym1u5 + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] = -1; ixyz[2] = -1; }
  // +1 -1  0
  x = (fxp1u0 + fym1u3 + fzu6  ) - X_factor;
  y = (fxp1u1 + fym1u4 + fzu7  ) - Y_factor;
  z = (fxp1u2 + fym1u5 + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] = -1; ixyz[2] =  0; }
  // +1 -1 +1
  x = (fxp1u0 + fym1u3 + fzp1u6) - X_factor;
  y = (fxp1u1 + fym1u4 + fzp1u7) - Y_factor;
  z = (fxp1u2 + fym1u5 + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] = -1; ixyz[2] =  1; }
  // +1  0 -1
  x = (fxp1u0 + fyu3   + fzm1u6) - X_factor;
  y = (fxp1u1 + fyu4   + fzm1u7) - Y_factor;
  z = (fxp1u2 + fyu5   + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] =  0; ixyz[2] = -1; }
  // +1  0  0
  x = (fxp1u0 + fyu3   + fzu6  ) - X_factor;
  y = (fxp1u1 + fyu4   + fzu7  ) - Y_factor;
  z = (fxp1u2 + fyu5   + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] =  0; ixyz[2] =  0; }
  // +1  0 +1
  x = (fxp1u0 + fyu3   + fzp1u6) - X_factor;
  y = (fxp1u1 + fyu4   + fzp1u7) - Y_factor;
  z = (fxp1u2 + fyu5   + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] =  0; ixyz[2] =  1; }
  // +1 +1 -1
  x = (fxp1u0 + fyp1u3 + fzm1u6) - X_factor;
  y = (fxp1u1 + fyp1u4 + fzm1u7) - Y_factor;
  z = (fxp1u2 + fyp1u5 + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] =  1; ixyz[2] = -1; }
  // +1 +1  0
  x = (fxp1u0 + fyp1u3 + fzu6  ) - X_factor;
  y = (fxp1u1 + fyp1u4 + fzu7  ) - Y_factor;
  z = (fxp1u2 + fyp1u5 + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] =  1; ixyz[2] =  0; }
  // +1 +1 +1
  x = (fxp1u0 + fyp1u3 + fzp1u6) - X_factor;
  y = (fxp1u1 + fyp1u4 + fzp1u7) - Y_factor;
  z = (fxp1u2 + fyp1u5 + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] =  1; ixyz[2] =  1; }

  //if (closest2 != 0.0 && min < closest2) return (min);
//  this->ClosestImage(a1, a2, ixyz);
//  fprintf(stdout,"DEBUG: Predict  = %2i %2i %2i\n",ixyz[0],ixyz[1],ixyz[2]);

//  ix = ixyz[0];
//  iy = ixyz[1];
//  iz = ixyz[2];

//D = sqrt(min);
//  fprintf(stdout,"DEBUG: MinDist  = %2i %2i %2i = %8.3f\n", ixmin, iymin, izmin, D);
//  printf("---------------------------------------------------------------\n");
  return(min);

}
