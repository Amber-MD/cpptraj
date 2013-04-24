#include <cmath> // sqrt
#include <cstdio> // sscanf
#include <cstring> // memcpy, memset
#include "DataSet_Modes.h"
#include "CpptrajStdio.h"
#include "BufferedFrame.h"

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
  avgcrd_(0),
  evalues_(0),
  evectors_(0),
  nmodes_(0),
  vecsize_(0),
  navgcrd_(0),
  reduced_(false)
{}

// DESTRUCTOR
DataSet_Modes::~DataSet_Modes() {
  if (avgcrd_ != 0) delete[] avgcrd_;
  if (evalues_!=0) delete[] evalues_;
  if (evectors_!=0) delete[] evectors_;
}

// DataSet_Modes::SetAvgCoords()
void DataSet_Modes::SetAvgCoords(int ncoords, const double* Xin) {
  if (avgcrd_!=0) delete[] avgcrd_;
  navgcrd_ = ncoords;
  if (navgcrd_ > 0) {
    avgcrd_ = new double[ navgcrd_ ];
    memcpy( avgcrd_, Xin, navgcrd_ * sizeof(double));
  } else
    avgcrd_ = 0;
}

/** Get eigenvectors and eigenvalues. They will be stored in descending 
  * order (largest eigenvalue first).
  */
int DataSet_Modes::CalcEigen(DataSet_Matrix& mIn, int n_to_calc) {
#ifdef NO_MATHLIB
  mprinterr("Error: modes: Compiled without ARPACK/LAPACK/BLAS routines.\n");
  return 1;
#else
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
    // If no eigenvectors calcd set vecsize to 0
    if (evectors_==0)
      vecsize_ = 0;
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
  // Eigenvalues and eigenvectors are in ascending order. Resort so that
  // they are in descending order (i.e. largest eigenvalue first).
  int pivot = nmodes_ / 2;
  int nmode = nmodes_ - 1;
  double* vtmp = 0;
  if (evectors_ != 0) 
    vtmp = new double[ vecsize_ ];
  for (int mode = 0; mode < pivot; ++mode) {
    // Swap eigenvalue
    double eval = evalues_[mode];
    evalues_[mode] = evalues_[nmode];
    evalues_[nmode] = eval;
    // Swap eigenvector
    if (vtmp != 0) {
      double* Vec0 = evectors_ + (mode  * vecsize_);
      double* Vec1 = evectors_ + (nmode * vecsize_);
      memcpy( vtmp, Vec0, vecsize_ * sizeof(double) );
      memcpy( Vec0, Vec1, vecsize_ * sizeof(double) );
      memcpy( Vec1, vtmp, vecsize_ * sizeof(double) );
    }
    --nmode;
  }
  if (vtmp != 0) delete[] vtmp;

  return 0;
#endif
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
int DataSet_Modes::WriteToFile(std::string const& fname) {
  if (fname.empty()) {
    mprinterr("Internal Error: DataSet_Modes::WriteToFile: No filename given.\n");
    return 1;
  }
  BufferedFrame outfile;
  if (outfile.OpenWrite( fname )) {
    mprinterr("Error: Could not open %s for writing.\n", fname.c_str());
    return 1;
  }
  if (reduced_)
    outfile.Printf(" Reduced Eigenvector file: ");
  else
    outfile.Printf(" Eigenvector file: ");
  outfile.Printf("%s", DataSet_Matrix::MatrixOutputString[type_]);
  // Write out # of modes on title line to not break compat. with older modes files
  outfile.Printf(" nmodes %i", nmodes_);
  // Write out col width on title line to not break compat. with older modes files
  // Since data format has leading space, actual width is width + 1
  int colwidth = width_ + 1;
  outfile.Printf(" width %i\n", colwidth);
  // First number is # avg coords, second is size of each vector
  outfile.Printf(" %4i %4i\n", navgcrd_, vecsize_);
  // Set up framebuffer, default 7 columns
  int bufsize;
  if (navgcrd_ > vecsize_)
    bufsize = navgcrd_;
  else
    bufsize = vecsize_;
  outfile.SetupFrameBuffer( bufsize, colwidth, 7 );
  // Print average coords
  outfile.DoubleToBuffer( avgcrd_, navgcrd_, data_format_ );
  outfile.WriteFrame();
  // Eigenvectors and eigenvalues
  for (int mode = 0; mode < nmodes_; ++mode) {
    outfile.Printf(" ****\n %4i ", mode+1);
    outfile.Printf(data_format_, evalues_[mode]);
    outfile.Printf("\n");
    if (evectors_ != 0) {
      const double* Vec = Eigenvector(mode);
      outfile.BufferBegin();
      outfile.DoubleToBuffer( Vec, vecsize_, data_format_ );
      outfile.WriteFrame(); 
    }
  }
  outfile.CloseFile();
  return 0;
}

// DataSet_Modes::ReadEvecFile()
int DataSet_Modes::ReadEvecFile(std::string const& modesfile, int ibeg, int iend) {
  int modesToRead = iend - ibeg + 1;
  if (modesToRead < 1) {
    mprinterr("Error: Specified # of modes to read (%i) must be > 0\n",modesToRead);
    return 1;
  }

  BufferedFrame infile;
  if (infile.OpenRead( modesfile)) return 1;
  // Read title line, convert to arg list
  const char* buffer = 0;
  if ( (buffer = infile.NextLine())==0 ) {
    mprinterr("Error: ReadEvecFile(): error while reading title (%s)\n",infile.Filename().full());
    return 1;
  }
  ArgList title(buffer);
  // Check if reduced
  reduced_ = title.hasKey("Reduced");
  // Determine modes file type
  type_ = DataSet_Matrix::NO_OP;
  for (int matidx = (int)DataSet_Matrix::NO_OP + 1; 
           matidx != (int)DataSet_Matrix::NMAT; ++matidx)
  {
    if ( title.hasKey( DataSet_Matrix::MatrixOutputString[matidx] ))
    {
      type_ = (DataSet_Matrix::MatrixType)matidx;
      break;
    }
  }
  // For compatibility with quasih and nmode output
  if (type_ == DataSet_Matrix::NO_OP) {
    mprintf("Warning: ReadEvecFile(): Unrecognized type [%s]\n", title.ArgLine());
    mprintf("         Assuming MWCOVAR.\n");
    type_ = DataSet_Matrix::MWCOVAR;
  }
  // For newer modesfiles, get # of modes in file.
  int modesInFile = title.getKeyInt("nmodes",-1);
  if (modesInFile == -1) {
    modesInFile = modesToRead; 
    mprintf("Warning: Older modes file, # of modes not known.\n");
    mprintf("Warning: Will try to read at least %i modes.\n", modesToRead);
  } else {
    mprintf("\tFile contains %i modes.\n", modesInFile);
    if (modesToRead > modesInFile) {
      mprintf("Warning: # modes to read (%i) > modes in file. Only reading %i modes.\n",
              modesToRead, modesInFile);
      modesToRead = modesInFile;
    }
  }
  // For newer modesfiles, get width of data elts
  int colwidth = title.getKeyInt("width", -1);
  if (colwidth == -1) 
    colwidth = 11; // Default, 10 + 1 space
  width_ = colwidth - 1;
  SetDataSetFormat(false);
  // Read number of elements in avg coords and eigenvectors
  if ( (buffer = infile.NextLine())==0 ) {
    mprinterr("Error: ReadEvecFile(): error while reading number of atoms (%s)\n",
              infile.Filename().full());
    return 1;
  }
  int nvals = sscanf(buffer, "%i %i", &navgcrd_, &vecsize_);
  if (nvals == 0) {
    mprinterr("Error: ReadEvecFile(): sscanf on coords failed (%s)\n",infile.Filename().full());
    return 1;
  } else if (nvals == 1) {
    mprintf("Warning: ReadEvecFile(): No value for eigenvector size found in %s,\n",
            infile.Filename().full());
    mprintf("         assuming it is equal to #average elements (%i)\n",navgcrd_);
    vecsize_ = navgcrd_;
  }
  // Allocate FrameBuffer
  int bufsize;
  if (navgcrd_ > vecsize_)
    bufsize = navgcrd_;
  else
    bufsize = vecsize_;
  infile.SetupFrameBuffer( bufsize, colwidth, 7 );
  // Allocate memory for avg coords and read in
  if (avgcrd_ != 0) delete[] avgcrd_;
  if (navgcrd_ > 0) {
    avgcrd_ = new double[ navgcrd_ ];
    infile.ReadFrame( );
    infile.BufferToDouble( avgcrd_, navgcrd_ );
    infile.BufferBegin(); // Reset buffer position
  }
  // Allocate memory for eigenvalues and eigenvectors
  if (evalues_!=0) delete[] evalues_;
  evalues_ = 0;
  if (evectors_!=0) delete[] evectors_;
  evectors_ = 0;
  evalues_ = new double[ modesToRead ];
  if (vecsize_ > 0) 
    evectors_ = new double[ modesToRead * vecsize_ ];
  nmodes_ = 0;
  int currentMode = 0;
  int nno = 0;
  bool firstRead = true;
  while ( (buffer = infile.NextLine())!=0 ) { // This should read in ' ****'
    if (strncmp(buffer," ****", 5)!=0) {
      mprinterr("Error: ReadEvecFile(): When reading eigenvector %i, expected ' ****',\n",
                currentMode+1);
      mprinterr("Error: got %s [%s]\n", buffer, infile.Filename().full());
      return 1;
    }
    // Read number and eigenvalue 
    if ( (buffer = infile.NextLine())==0 ) {
      mprinterr("Error: ReadEvecFile(): error while reading number and eigenvalue (%s)\n",
                infile.Filename().full());
      return 1;
    }
    if (sscanf(buffer, "%i%lf", &nno, evalues_ + nmodes_) != 2) {
      mprinterr("Error: ReadEvecFile(): error while scanning number and eigenvalue (%s)\n",
                infile.Filename().full());
      return 1;
    }
    if (vecsize_ > 0) {
      // Read eigenvector
      // Older modesfiles could have vecsize > 0 but no eigenvectors, only 
      // blanks. If this is the case set vecsize to -1 to indicate a blank
      // read is needed after reading eigenvalues.
      double* Vec = evectors_ + (nmodes_ * vecsize_);
      int vi = 0;
      while (vi < vecsize_) {
        buffer = infile.NextLine();
        if (firstRead && (buffer[0] == '\n' || buffer[0] == '\r')) {
          mprintf("Warning: Old modes file with vecsize > 0 but no eigenvectors.\n");
          vecsize_ = -1;
          delete[] evectors_;
          evectors_ = 0;
          break;
        }
        double tmpval[7];
        int nvals = sscanf(buffer, "%lf %lf %lf %lf %lf %lf %lf", tmpval, 
                           tmpval+1, tmpval+2, tmpval+3, tmpval+4, tmpval+5, tmpval+6);
        for (int ti = 0; ti < nvals; ++ti)
          Vec[vi++] = tmpval[ti];
      }
      // Check if mode read was between ibeg and iend (which start from 1).
      // If so, increment number of modes.
      if (currentMode+1 >= ibeg && currentMode < iend) ++nmodes_;
      if (nmodes_ == modesToRead) break; 
      ++currentMode;
    } else if (vecsize_ == -1) {
      // Blank read past empty eigenvector
      buffer = infile.NextLine();
    }
    firstRead = false;
  }
  if (nmodes_ != modesToRead) 
    mprintf("Warning: Number of read modes is %i, requested %i\n", nmodes_, modesToRead);
  return 0;
}

/** Convert eigenvalues to cm^-1. 
  * Frequency = sqrt( Ene / (mass * MSF)) = sqrt( Ene / Eigval )
  */
int DataSet_Modes::EigvalToFreq() {
  mprintf("\tConverting eigenvalues to frequencies (cm^-1).\n");
  for (int i = 0; i < nmodes_; ++i) {
    // "0.6" is conversion of kT for 300K into kcal/mol
    // "108.597" is conversion to units of cm^-1
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
