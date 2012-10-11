#include <cmath> // sqrt
#include <cstring> // memcpy, memset
#include "DataSet_Modes.h"
#include "CpptrajStdio.h"
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
  evalues_(0),
  evectors_(0),
  nmodes_(0),
  vecsize_(0)
{}

// DESTRUCTOR
DataSet_Modes::~DataSet_Modes() {
  if (evalues_!=0) delete[] evalues_;
  if (evectors_!=0) delete[] evectors_;
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
    //       w.r.t. C/C++ so eigenvectors are in rows.
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

/** Convert eigencalues to cm^-1 */
int DataSet_Modes::EigvalToFreq() {
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
  * must be == number of modes (TODO: Make std::vector).
  */
int DataSet_Modes::MassWtEigvect(const double* massIn) {
  if (massIn == 0) return 1;
  if (evectors_ == 0) return 0;
  int vend = nmodes_ * vecsize_;
  const double* mptr = massIn;
  for (int i = 0; i < nmodes_; ++i) {
    double mass = 1.0 / sqrt( *(mptr++) );
    for (int v = i * 3; v < vend; v += vecsize_) {
      evectors_[v  ] *= mass; 
      evectors_[v+1] *= mass; 
      evectors_[v+2] *= mass;
    }
  }
  return 0;
}
 
