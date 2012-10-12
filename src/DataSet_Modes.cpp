#include <cmath> // sqrt
#include <cstring> // memcpy, memset
#include "DataSet_Modes.h"
#include "CpptrajStdio.h"
#include "FrameBuffer.h"
//#incl ude "vectormath.h" // printModes

#ifndef NO_MATHLIB
// Definition of Fortran subroutines called from this class
extern "C" {
  // LAPACK
  void dspev_(char&, char&, int&, double*, double*, double*, int&, double*, int&);
  // ARPACK
  void dsaupd_(int&, char&, int&, char*, int&, double&, double*,
               int&, double*, int&, int*, int*, double*, double*,
               int&, int&);
  void dseupd_(int&, char&, int*, double*, double*, int&, double&,
               char&, int&, char*, int&, double&, double*,
               int&, double*, int&, int*, int*, double*, double*,
               int&, int&);
}
#endif

// CONSTRUCTOR
DataSet_Modes::DataSet_Modes() :
  DataSet(MODES, 10, 5, 0), // 0 dim indicates DataSet-specific write
  evalues_(0),
  evectors_(0),
  nmodes_(0),
  vecsize_(0),
  reduced_(false)
{}

// DESTRUCTOR
DataSet_Modes::~DataSet_Modes() {
  if (evalues_!=0) delete[] evalues_;
  if (evectors_!=0) delete[] evectors_;
}

// DataSet_Modes::SetAvgCoords()
void DataSet_Modes::SetAvgCoords(int ncoords, const double* Xin) {
  const double* Xend = Xin + ncoords;
  for (const double* Xptr = Xin; Xptr < Xend; Xptr += 3)
    avg_.AddXYZ( Xptr );
}

/** Get eigenvectors and eigenvalues. */
int DataSet_Modes::CalcEigen(DataSet_Matrix& mIn, int n_to_calc) {
  bool eigenvaluesOnly;
  int info = 0;
  if (mIn.Nrows() > 0) {
    mprinterr("Error: DataSet_Modes: Eigenvector/value calc only for symmetric matrices.\n");
    return 1;
  }
  // If number to calc is 0, assume we want eigenvalues only
  if (n_to_calc < 1) {
    eigenvaluesOnly = true;
    nmodes_ = mIn.Ncols();
  } else {
    eigenvaluesOnly = false;
    nmodes_ = n_to_calc;
  }
  if (nmodes_ > mIn.Ncols()) {
    mprintf("Warning: Specified # of eigenvalues to calc (%i) > matrix dimension (%i).\n",
            nmodes_, mIn.Ncols());
    nmodes_ = mIn.Ncols();
    mprintf("Warning: Only calculating %i eigenvalues.\n", nmodes_);
  }
  // -----------------------------------------------------------------
  if (nmodes_ == mIn.Ncols()) {
    // Calculate all eigenvalues (and optionally eigenvectors). 
    char jobz = 'V'; // Default: Calc both eigenvectors and eigenvalues
    vecsize_ = mIn.Ncols();
    // Check if only calculating eigenvalues
    if (eigenvaluesOnly) {
      jobz = 'N';
      vecsize_ = 1;
    }
    // Set up space to hold eigenvectors
    if (evectors_ != 0) delete[] evectors_;
    if (!eigenvaluesOnly)
      evectors_ = new double[ nmodes_ * vecsize_ ];
    else
      evectors_ = 0;
    // Set up space to hold eigenvalues
    if (evalues_ != 0) delete[] evalues_;
    evalues_ = new double[ nmodes_ ];
    // Create copy of matrix since call to dspev destroys it
    double* mat = new double[ mIn.Size() ];
    memcpy(mat, mIn.MatrixPtr(), mIn.Size() * sizeof(double));
    // Lower triangle; not upper since fortran array layout is inverted w.r.t. C/C++
    char uplo = 'L'; 
    // Allocate temporary workspace
    double* work = new double[ 3 * nmodes_ ];
    // NOTE: The call to dspev is supposed to store eigenvectors in columns. 
    //       However as mentioned above fortran array layout is inverted
    //       w.r.t. C/C++ so eigenvectors end up in rows.
    // NOTE: Eigenvalues/vectors are returned in ascending order.
    dspev_(jobz, uplo, nmodes_, mat, evalues_, evectors_, vecsize_, work, info);
    delete[] work;
    delete[] mat;
    if (info != 0) {
      if (info < 0) {
        mprinterr("Internal Error: from dspev: Argument %i had illegal value.\n", -info);
        mprinterr("Args: %c %c %i matrix %x %x %i work %i\n", jobz, uplo, nmodes_,  
                  evalues_, evectors_, vecsize_, info);
      } else { // info > 0
        mprinterr("Internal Error: from dspev: The algorithm failed to converge.\n");
        mprinterr("%i off-diagonal elements of an intermediate tridiagonal form\n", info);
        mprinterr("did not converge to zero.\n");
      }
      return 1;
    }
  // -----------------------------------------------------------------
  } else {
    // Calculate up to n-1 eigenvalues/vectors using the Implicitly Restarted
    // Arnoldi iteration.
    // FIXME: Eigenvectors obtained with this method appear to have signs
    //        flipped compared to full method - is dot product wrong?
    int nelem = mIn.Ncols(); // Dimension of input matrix (N)
    // Allocate memory to store eigenvectors
    vecsize_ = mIn.Ncols();
    int ncv; // # of columns of the matrix V (evectors_), <= N (mIn.Ncols())
    if (evectors_!=0) delete[] evectors_;
    if (nmodes_*2 <= nelem) 
      ncv = nmodes_*2;
    else 
      ncv = nelem;
    evectors_ = new double[ ncv * nelem ];
    // Allocate memory to store eigenvalues
    if ( evalues_ != 0) delete[] evalues_;
    evalues_ = new double[ nelem ] ; // NOTE: Should this be nmodes?
    // Allocate workspace
    double* workd = new double[ 3 * nelem ];
    int lworkl = ncv * (ncv+8); // NOTE: should this be ncv^2 * (ncv+8)
    double* workl = new double[ lworkl ];
    double* resid = new double[ nelem ];
    // Set parameters for dsaupd (Arnolid)
    int ido = 0; // Reverse comm. flag; 0 = first call
    // The iparam array is used to set parameters for the calc.
    int iparam[11];
    memset(iparam, 0, 11 * sizeof(int));
    iparam[0] = 1;   // Method for selecting implicit shifts; 1 = exact
    iparam[2] = 300; // Max # of iterations allowed
    iparam[3] = 1;   // blocksize to be used in the recurrence (code works with only 1).
    iparam[6] = 1;   // Type of eigenproblem being solved; 1: A*x = lambda*x
    double tol = 0;  // Stopping criterion (tolerance); 0 = arpack default 
    char bmat = 'I'; // Type of matrix B that defines semi-inner product; I = identity
    char which[2];   // Which of the Ritz values of OP to compute;
    which[0] = 'L';  // 'LA' = compute the NEV largest eigenvalues
    which[1] = 'A';
    // The ipntr array will hold starting locations in workd and workl arrays
    // for matrices/vectors used by the Lanczos iteration.
    int ipntr[11];
    memset(ipntr, 0, 11 * sizeof(int));
    // Create copy of matrix since it will be modified 
    double* mat = new double[ mIn.Size() ];
    memcpy(mat, mIn.MatrixPtr(), mIn.Size() * sizeof(double));
    // LOOP
    bool loop = false;
    do {
      if (loop) {
        // Dot products
        double* target = workd + (ipntr[1] - 1); // -1 since fortran indexing starts at 1
        double* vec    = workd + (ipntr[0] - 1);
        memset(target, 0, nelem*sizeof(double));
        for(int i = 0; i < nelem; i++) {
          for(int j = i; j < nelem; j++) {
            int ind = nelem * i + j - (i * (i + 1)) / 2;
            target[i] += mat[ind] * vec[j];
            if(i != j)
              target[j] += mat[ind] * vec[i];
          }
        }
      }

      dsaupd_(ido, bmat, nelem, which, nmodes_, tol, resid,
              ncv, evectors_, nelem, iparam, ipntr, workd, workl,
              lworkl, info);
      loop = (ido == -1 || ido == 1);
    } while ( loop ); // END LOOP

    if (info != 0) {
      mprinterr("Error: DataSet_Modes: dsaupd returned %i\n",info);
    } else {
      int rvec = 1;
      char howmny = 'A';
      double sigma = 0.0;
      int* select = new int[ ncv ];
      dseupd_(rvec, howmny, select, evalues_, evectors_, nelem, sigma,
              bmat, nelem, which, nmodes_, tol, resid,
              ncv, evectors_, nelem, iparam, ipntr, workd, workl,
              lworkl, info);
      delete[] select;
    } 
    delete[] mat;
    delete[] workl;
    delete[] workd;
    delete[] resid;
    if (info != 0) { 
      mprinterr("Error: DataSet_Modes: dseupd returned %i\n",info);
      return 1;
    }
  }
  return 0;
}

// DataSet_Modes::PrintModes()
void DataSet_Modes::PrintModes() {
  mprintf("%s: %i modes.\n",Legend().c_str(),nmodes_);
  for (int mode = 0; mode < nmodes_; ++mode) {
    mprintf("Mode %i: Eigenvalue= %f\n", mode, evalues_[mode]);
    if (evectors_!=0) {
      mprintf("\tEigenvector={");
      const double* Vec = Eigenvector(mode);
      for (int veci = 0; veci < vecsize_; ++veci) 
        mprintf(" %f", Vec[veci]);
      mprintf(" }\n");
    }
  }
  //printMatrix("Eigenvectors (Rows):", evectors_, nmodes_, vecsize_);
}

/** Write eigenvalues and if present eigenvectors/avg coords to file
  * in PTRAJ-compatible format.
  */
void DataSet_Modes::WriteToFile(CpptrajFile& outfile) {
  if (!outfile.IsOpen()) {
    mprinterr("Internal Error: DataSet_Modes: File %s is not open for write.\n",
              outfile.FullFileStr());
    return;
  }
  if (reduced_)
    outfile.Printf(" Reduced Eigenvector file: ");
  else
    outfile.Printf(" Eigenvector file: ");
  outfile.Printf("%s\n", DataSet_Matrix::MatrixOutputString[type_]);
  // First number is # avg coords, second is size of each vector
  outfile.Printf(" %4i %4i\n", avg_.size(), vecsize_);
  // Set up framebuffer, default 7 columns
  // Since data format has leading space, actual width is width + 1
  int colwidth = width_ + 1;
  int bufsize;
  if (avg_.size() > vecsize_)
    bufsize = avg_.size();
  else
    bufsize = vecsize_;
  FrameBuffer fbuffer(bufsize, colwidth, 7, outfile.IsDos());
  // Print average coords
  fbuffer.DoubleToBuffer( avg_.xAddress(), avg_.size(), data_format_, colwidth, 7);
  outfile.Write( fbuffer.Buffer(), fbuffer.CurrentSize() );
  // Eigenvectors and eigenvalues
  // TODO: Reverse order of eigenvalues prior to this call.
  int imode = 1;
  for (int mode = nmodes_ - 1; mode >= 0; --mode) {
    outfile.Printf(" ****\n %4i ", imode++);
    outfile.Printf(data_format_, evalues_[mode]);
    outfile.Printf("\n");
    if (evectors_ != 0) {
      const double* Vec = Eigenvector(mode);
      fbuffer.BufferBegin();
      fbuffer.DoubleToBuffer( Vec, vecsize_, data_format_, colwidth, 7 );
      outfile.Write( fbuffer.Buffer(), fbuffer.CurrentSize() ); 
    }
  }
}

/** Convert eigenvalues to cm^-1 */
int DataSet_Modes::EigvalToFreq() {
  mprintf("\tConverting eigenvalues to frequencies.\n");
  for (int i = 0; i < nmodes_; ++i) {
    // "0.6" is conversion of kT for 300K into kcal/mol(?)
    if (evalues_[i] > 0)
      evalues_[i] =  108.587 * sqrt( 0.6 / evalues_[i]);
    else if (evalues_[i] < 0.0)
      evalues_[i] = -108.587 * sqrt(-0.6 / evalues_[i]);
    else {
      mprinterr("Error: DataSet_Modes: bad eigenvalue %i = %f\n", i, evalues_[i]);
      return 1;
    }
  }
  return 0;
}

/** Mass-weght Eigenvectors. Currently only works when vector size
  * is a multiple of 3 (i.e. COVAR-type matrix. Size of massIn
  * must be == number of modes (TODO: Make std::vector). The
  * ith xyz elements of each eigenvector is multiplied by mass i.
  */
int DataSet_Modes::MassWtEigvect(const double* massIn) {
  if (massIn == 0) return 1;
  if (evectors_ == 0) return 0;
  mprintf("\tMass-weighting %i eigenvectors\n", nmodes_);
  int vend = nmodes_ * vecsize_; // == size of evectors array
  const double* mptr = massIn;
  for (int vi = 0; vi < vecsize_; vi += 3) {
    double mass = 1.0 / sqrt( *(mptr++) );
    for (int modev = vi; modev < vend; modev += vecsize_) {
//      mprinterr("evectors[%i] *= %f\nevectors[%i] *= %f\nevectors[%i] *= %f\n", // DEBUG
//                modev,mass,modev+1,mass,modev+2,mass); // DEBUG
      evectors_[modev  ] *= mass;
      evectors_[modev+1] *= mass;
      evectors_[modev+2] *= mass;
    }
  }
  return 0;
}

/** Reduce covariance eigenvectors. Each eigenvector is assumed to have
  * X, Y, and Z components. Each eigenvector element is reduced via
  * Ei = Eix^2 + Eiy^2 + Eiz^2. See Abseher & Nilges, JMB 1998, 279, 911-920.
  */
int DataSet_Modes::ReduceCovar() {
  if (evectors_ == 0) {
    mprinterr("Error: reduce: No eigenvectors present.\n");
    return 1;
  }
  int newvecsize = vecsize_ / 3;
  mprintf("\tReducing size of %i eigenvectors from %i to %i\n",nmodes_,vecsize_,newvecsize);
  double* newEvectors = new double[ nmodes_ * newvecsize ];
  for (int mode = 0; mode < nmodes_; ++mode) {
    const double* Vec = Eigenvector(mode);
    double* newVec = newEvectors + (mode * newvecsize);
    for (int vi = 0; vi < vecsize_; vi += 3) { 
      //mprinterr("newVec[%u]=%f*%f + %f*%f + %f*%f\n",newVec-(newEvectors + (mode * newvecsize)),
      //          Vec[vi],Vec[vi],Vec[vi+1],Vec[vi+1],Vec[vi+2],Vec[vi+2]); // DEBUG
      *(newVec++) = Vec[vi]*Vec[vi] + Vec[vi+1]*Vec[vi+1] + Vec[vi+2]*Vec[vi+2];
    }
  }
  delete[] evectors_;
  evectors_ = newEvectors;
  vecsize_ = newvecsize;
  reduced_ = true;
  return 0;
}

/** Reduce distance covariance eigenvectors. Each eigenvector element 
  * corresponds to a different atom pair. E.g., for 4 atoms there are 
  * 6 possible pairs: {0[0,1], 1[0,2], 2[0,3], 3[1,2], 4[1,3], 5[2,3]}
  * which in a symmetric half-matrix (no diagonal) looks like:
  *   X 0 1 2
  *   0 X 3 4
  *   1 3 X 5
  *   2 4 5 X
  * Eigenvectors are reduced by taking the sum of the squares of each row:
  * 0[0^2 + 1^2 + 2^2], 1[0^2 + 3^2 + 4^2], 2[1^2 + 3^2 + 5^2], etc
  */
int DataSet_Modes::ReduceDistCovar(int nelts) {
  int i, j;
  if (evectors_ == 0) {
    mprinterr("Error: reduce: No eigenvectors present.\n");
    return 1;
  }
  // TODO: Check that nelts * (nelts-1) / 2 == vecsize
  int newvecsize = nelts;
  mprintf("\tReducing size of %i eigenvectors from %i to %i\n",nmodes_,vecsize_,newvecsize);
  double* newEvectors = new double[ nmodes_ * newvecsize ];
  double* newVec = newEvectors;
  for (int mode = 0; mode < nmodes_; ++mode) {
    const double* Vec = Eigenvector(mode);
    for (int row = 0; row < nelts; ++row) {
      *newVec = 0.0;
      for (int col = 0; col < nelts; ++col) {
        if (row != col) {
          // Calculate distance index into half-matrix w.o. diagonal,
          // see TriangleMatrix::calcIndex
          if (row > col) {
            j = row;
            i = col;
          } else {
            i = row;
            j = col;
          }
          int i1 = i + 1;
          double v = Vec[ ( (nelts * i) - ((i1 * i) / 2) ) + j - i1 ];
          *newVec += (v * v);
        }
      }
      ++newVec;
    }
  }
  delete[] evectors_;
  evectors_ = newEvectors;
  vecsize_ = newvecsize;
  reduced_ = true;
  return 0;
}  
