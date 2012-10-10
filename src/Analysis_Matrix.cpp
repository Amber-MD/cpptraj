#include <cmath> // sqrt
#include "Analysis_Matrix.h"
#include "CpptrajStdio.h"
#include "Thermo.h"

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
Analysis_Matrix::Analysis_Matrix() :
  minfo_(0),
  thermopt_(false),
  nevec_(0),
  reduce_(false),
  modinfo_(0),
  freeWork_(true),
  eigval_(0),
  vout_(0)
{}

// DESTRUCTOR
Analysis_Matrix::~Analysis_Matrix() {
  if (freeWork_) {
    if (eigval_!=0) delete[] eigval_;
    if (vout_!=0) delete[] vout_;
  }
}

// Analysis_Matrix::WorkspaceStored()
/** This is called if no errors in analysis and workspace information is 
  * successfully stored in the modinfo structure. The modinfo structure
  * is always placed inside the master DataSetList so it becomes responsible
  * for freeing the workspace.
  */
void Analysis_Matrix::WorkspaceStored() {
  freeWork_ = false;
}

// Analysis_Matrix::Setup()
/** analyze matrix <matrixname> [out <filename>] [name <name>] [thermo | order] 
  *                [vecs <vecs>] [reduce]
  */
int Analysis_Matrix::Setup(DataSetList* DSLin) {
#ifdef NO_MATHLIB
  mprinterr("Error: analyze matrix: Compiled without LAPACK routines.\n");
  return 1;
#else
  // Ensure first 2 args (should be 'analyze' 'matrix') are marked.
  analyzeArgs_.MarkArg(0);
  analyzeArgs_.MarkArg(1);
  // Get matrix name
  std::string mname = analyzeArgs_.GetStringNext();
  if (mname.empty()) {
    mprinterr("Error: analyze matrix: missing matrix name (first argument).\n");
    return 1;
  }
  // Find matrix in DataSetList.
  minfo_ = (DataSet_Matrix*)DSLin->FindSetOfType( mname, DataSet::MATRIX );
  if (minfo_ == 0) {
    mprinterr("Error: analyze matrix: Could not find matrix named %s\n",mname.c_str());
    return 1;
  }
  // Filenames
  outfilename_ = analyzeArgs_.GetStringKey("out");
  // Thermo flag
  thermopt_ = analyzeArgs_.hasKey("thermo");
  if (thermopt_)
    outthermo_ = analyzeArgs_.GetStringKey("outthermo");
  if (thermopt_ && minfo_->Mass()==0) {
    mprinterr("Error: analyze matrix: parameter 'thermo' only works for\n");
    mprinterr("       mass-weighted covariance matrix ('mwcovar').\n");
    return 1;
  }
  // Number of eigenvectors; allow "0" only in case of 'thermo'
  nevec_ = analyzeArgs_.getKeyInt("vecs",0);
  if (thermopt_) {
    if (nevec_ < 0) nevec_ = 0;
  } else if (nevec_ <= 0) {
    nevec_ = 1;
  }
  // Reduce flag
  reduce_ = analyzeArgs_.hasKey("reduce");
  /*if ( reduce_ && minfo_->Type() != MatrixType::MATRIX_MWCOVAR &&
                  minfo_->Type() != MatrixType::MATRIX_COVAR &&
                  minfo_->Type() != MatrixType::MATRIX_DISTCOVAR)*/
  if ( reduce_ && minfo_->Ncols() == minfo_->Nelts() )
  {
    mprinterr("Error: analyze matrix: parameter 'reduce only works for\n");
    mprinterr("       covariance and distance covariance matrices.\n");
    return 1;
  }
  // Set up ModesInfo
  std::string modesname = analyzeArgs_.GetStringKey("name");
  // Unlike PTRAJ, always put ModesInfo into the DataSetList; this way
  // DataSetList always handles the destruction.
  modinfo_ = (ModesInfo*)DSLin->AddSet( DataSet::MODES, modesname, "Modes" );
  // TODO: Pass in MatrixType Nelt
  if (modinfo_==0) return 1;

  // Print Status
  mprintf("    ANALYZE MATRIX: Analyzing matrix %s",minfo_->Legend().c_str());
  if (!outfilename_.empty())
    mprintf(" and dumping results to %s\n", outfilename_.c_str());
  else
    mprintf(" and printing results to STDOUT\n");
  mprintf("      Calculating %i eigenvectors", nevec_);
  if (thermopt_)
    mprintf(" and thermodynamic data");
  mprintf("\n");
  if (nevec_>0 && reduce_)
    mprintf("      Eigenvectors will be reduced\n");
  if (!modesname.empty())
    mprintf("      Storing modes with name: %s\n",modesname.c_str());
  return 0;
#endif
}

#ifndef NO_MATHLIB
// dotprod()
// NOTE: Only needed when using dsaupd
static void dotprod(int nelem, double *mat, double *vec, double *target){

  int i, j, ind;
  for(i = 0; i < nelem; i++){
    target[i] = 0.0;
  }

  for(i = 0; i < nelem; i++){
    for(j = i; j < nelem; j++){
      ind = nelem * i + j - (i * (i + 1)) / 2;
      target[i] += mat[ind] * vec[j];
      if(i != j)
        target[j] += mat[ind] * vec[i];
    }
  }
}

// calcIndex()
/** Assuming an upper-right half matrix (no diagonal) calculate the proper
  * index for a 2d array. SHOULD NEVER BE CALLED WITH iIn == jIn!
  * Adapted from TriangleMatrix::calcIndex
  */
static int calcIndex(int nrows, int iIn, int jIn) {
  int i, j, i1;

  if (iIn > jIn) {
    j = iIn;
    i = jIn;
  } else {
    i = iIn;
    j = jIn;
  }

  i1 = i + 1;
  return ( ( (nrows * i) - ((i1 * i) / 2) ) + j - i1 );
}
#endif

// Analysis_Matrix::Analyze()
int Analysis_Matrix::Analyze() {
#ifdef NO_MATHLIB
  mprinterr("Error: Compiled without LAPACK routines.\n");
  return 1;
#else
  // Find eigenvalues and eigenvectors
  int neval = 0;
  int info = 0;
  double* workd = 0;

  const double* vect = minfo_->Vect();
  // TODO: Just pass in matrix Nelt?
  int nelem = minfo_->Ncols();
  //modinfo_->SetNavgElem( nelem ); // NOTE: Now done below when avg crds are set in ModesInfo
  //mprintf("CDBG: nevec=%i  nelem=%i\n",nevec_, nelem);
  if (nevec_ > nelem) {
    nevec_ = nelem;
    mprintf("Warning: NEVEC > NELEM: Number of calculated evecs were reduced to %i\n",nevec_);
  }
  if (nevec_ > minfo_->Nsnap()) {
    nevec_ = minfo_->Nsnap();
    mprintf("Warning: NEVEC > #Snapshots: Number of calculated evecs were reduced to %i\n",nevec_);
  }

  if (nevec_ == 0 || nelem == nevec_) {
    char uplo, jobz;
    int ldz = 0;

    neval = nelem;
    eigval_ = new double[ nelem ];
    // scratch for the function dspev
    workd = new double[ 3 * nelem ];
    // lower triangle of matrix is expected as input for dspev
    uplo = 'L';

    if (nevec_ == nelem) {
      // Get all eigenvectors
      vout_ = new double[ nelem * nelem ];
      // calculate eigenvalues and eigenvectors
      jobz = 'V';
      // Dimension of vout
      ldz = nelem;
    } else if (nevec_ == 0) {
      // Get only eigenvalues, only possible if thermo flag is set, otherwise
      // nevec is set to 1 as default.
      vout_ = new double[ nelem ];
      // only calculate eigenvalues
      jobz = 'N';
      // Dimension of vout
      ldz = 1;
    }

    dspev_(jobz, uplo, nelem, minfo_->MatrixPtr(), eigval_, vout_, ldz, workd, info);

    delete[] workd;
    if (info!=0) {
      mprinterr("Error: analyze matrix: dspev returned %i\n",info);
      return 1;
    }
  } else {
    // get up to n-1 eigenvectors
    int ncv;
    neval = nevec_;
    if (nevec_*2 <= nelem) {
      ncv = nevec_ * 2;
      vout_ = new double[ nelem * ncv ];
    } else {
      ncv = nelem;
      if (reduce_)
        vout_ = new double[ nelem * 2 * ncv ];
      else
        vout_ = new double[ nelem * ncv ];
    }
    eigval_ = new double[ nelem ];

    workd = new double[ 3 * nelem ];
    int lworkl = ncv * (ncv+8);
    double* workl = new double[ lworkl ];
    double* resid = new double[ nelem ];
    int* select = new int[ ncv ];
    
    int ido = 0;
    info = 0;
    int ipntr[11];
    int iparam[11];
    for (int i = 0; i < 11; ++i) {
      ipntr[i] = 0;
      iparam[i] = 0;
    }
    iparam[0] = 1;
    iparam[2] = 300;
    iparam[3] = 1;
    iparam[6] = 1;
    double tol = 0;
    char bmat = 'I';
    char which[2];
    which[0] = 'L';
    which[1] = 'A';
    bool loop = false;

    //mprintf("CDBG: At do loop, nelem=%i\n", nelem);
    do {
      if (loop)
        dotprod(nelem, minfo_->MatrixPtr(), workd + (ipntr[0] - 1), workd + (ipntr[1] - 1));

      dsaupd_(ido, bmat, nelem, which, nevec_, tol, resid,
              ncv, vout_, nelem, iparam, ipntr, workd, workl,
              lworkl, info);
     
      loop = (ido == -1 || ido == 1); 
    } while ( loop );

    if (info != 0) {
      mprinterr("Error: analyze matrix: dsaupd returned %i\n",info);
    } else {
      int rvec = 1;
      char howmny = 'A';
      double sigma = 0.0;

      dseupd_(rvec, howmny, select, eigval_, vout_, nelem, sigma,
              bmat, nelem, which, nevec_, tol, resid,
              ncv, vout_, nelem, iparam, ipntr, workd, workl,
              lworkl, info);
      if (info!=0)
        mprinterr("Error: analyze matrix: dseupd returned %i\n",info);
    }
    delete[] workl;
    delete[] workd;
    delete[] resid;
    delete[] select;
    if (info != 0) {
      return 1;
    }
  }

  //if (minfo_->Type() == MatrixType::MATRIX_MWCOVAR) {
  if (minfo_->Mass() != 0) {
    // Convert eigenvalues to cm^-1
    for (int i = 0; i < neval; ++i) {
      if (eigval_[i] > 0)
        // "0.6" is conversion of kT for 300K into kcal/mol(?)
        eigval_[i] =  108.587 * sqrt( 0.6 / eigval_[i]);
      else if (eigval_[i] < 0.0)
        eigval_[i] = -108.587 * sqrt(-0.6 / eigval_[i]);
      else {
        mprinterr("Error: analyze matrix: bad eigenvalue %i = %lf\n", i, eigval_[i]);
        return 1;
      }
    }
    // Mass-weight eigenvectors 
    //int crow = 0;
    const double* mptr = minfo_->Mass();
    for ( int crow = 0; crow < minfo_->Ncols(); crow += 3) { // This will ONLY work for COVAR
    /*for (AtomMask::const_iterator atomi = minfo_->Mask1().begin(); 
                                  atomi != minfo_->Mask1().end(); ++atomi)
    {*/
      double mass = 1.0 / sqrt( *(mptr++) );
      for (int i = 0; i < nevec_; ++i) {
        int idx = i*nelem + crow;
        vout_[idx  ] *= mass;
        vout_[idx+1] *= mass;
        vout_[idx+2] *= mass;
      }
    }
    // Calc thermo chemistry if specified
    if (thermopt_) {
      double* eigvali = new double[ neval ];
      double* vibn = new double[ 4 * neval ];
      //int natoms = minfo_->Mask1Tot();
      int natoms = minfo_->Ncols() / 3; // This will ONLY work for COVAR
      //double* masses = new double[ natoms ];
      for (int i = 0; i < neval; ++i)
        eigvali[i] = eigval_[neval - 1 - i];
      /*crow = 0;
      for (AtomMask::const_iterator atomi = minfo_->Mask1().begin();
                                  atomi != minfo_->Mask1().end(); ++atomi)
        masses[crow++] = (*(minfo_->Parm()))[*atomi].Mass();*/
      CpptrajFile outfile;
      outfile.OpenWrite(outthermo_);
      thermo( outfile, natoms, neval, 1, vect, minfo_->Mass(), eigvali, 
              vibn, vibn + 1*neval, vibn + 2*neval, vibn + 3*neval,
              298.15, 1.0 );
      outfile.CloseFile();
      delete[] eigvali;
      delete[] vibn;
      //delete[] masses; 
    }
  }

  // Reduction of eigenvectors (s. Abseher & Nilges, JMB 1998, 279, 911-920.)
  double* voutput = 0;
  int nout = 0;
  if (reduce_) {
    voutput = vout_ + nevec_ * nelem;
    nout = minfo_->Nelts();
    mprinterr("CDBG: reduce: nout = %i\n", nout);
    if ( nout * 3 == minfo_->Ncols() ) { // 3 cols per elt, assume COVAR
    /*if (minfo_->Type()==MatrixType::MATRIX_COVAR ||
        minfo_->Type()==MatrixType::MATRIX_MWCOVAR)
    {*/
      for (int i = 0; i < nevec_; ++i) {
        for (int j = 0; j < nout; ++j) {
          int idx = i*nelem + j*3;
          voutput[i*nout + j] = vout_[idx  ] * vout_[idx  ] +
                                vout_[idx+1] * vout_[idx+1] +
                                vout_[idx+2] * vout_[idx+2];
        }
      }
    } else { // Otherwise assume DISTCOVAR - could be dangerous...
    //if (minfo_->Type()==MatrixType::MATRIX_DISTCOVAR) {
      for (int i = 0; i < nevec_; ++i) {
        for (int j = 0; j < nout; ++j) {
          int vidx = i*nout + j;
          voutput[vidx] = 0;
          for (int k = 0; k < nout; ++k) {
            if (k != j) {
              int idx = calcIndex( nout, j, k );
              double vecelem = vout_[i*nelem + idx];
              voutput[vidx] += (vecelem * vecelem);
            }
          }
        }
      }
    }
  } else {
    voutput = vout_;
    nout = nelem;
  }

  // ---------- Output average coordinates, eigenvectors, and eigenvalues
  if (!outfilename_.empty()) {
    CpptrajFile outfile;
    if (outfile.SetupWrite(outfilename_, debug_)) {
      mprinterr("Error opening matrix analysis output file %s\n",outfilename_.c_str());
      return 1;
    }
    outfile.OpenFile();
    if (reduce_)
      outfile.Printf(" Reduced Eigenvector file: ");
    else
      outfile.Printf(" Eigenvector file: ");
    //outfile.Printf("%s\n", MatrixType::MatrixOutput[minfo_->Type()]);
    outfile.Printf("TEMP\n");
    outfile.Printf(" %4i %4i\n", nelem, nout);
    //mprintf("CDBG: neval = %i  nevec = %i  nout = %i\n", neval,nevec_, nout);
    //if (minfo_->Type()!=MatrixType::MATRIX_DIST) {
      for (int i = 0; i < nelem; ++i) {
        outfile.Printf(" %10.5f", vect[i]);
        if ((i+1)%7 == 0)
          outfile.Printf("\n");
      }
      if (nelem%7 != 0)
        outfile.Printf("\n");
    //}
    // Eigenvectors and eigenvalues
    for (int i = neval - 1; i >= 0; i--) {
      outfile.Printf(" ****\n %4i  %10.5f\n", neval - i, eigval_[i]);
      if (nevec_ > 0) {
        for (int j = 0; j < nout; ++j) {
          outfile.Printf(" %10.5f", voutput[i * nout + j]);
          if((j+1) % 7 == 0)
            outfile.Printf("\n");
        }
      }
      if(nout % 7 != 0)
        outfile.Printf("\n");
    }
    outfile.CloseFile();
  }

  // Get same ordering of arrays as in output if storing on stack requested
  if (modinfo_->Source() == ModesInfo::MS_STACK) {
    int k = (int) (neval / 2);
    for(int i = 0; i < k; ++i){
      int idx = neval - 1 - i;
      double tmp = eigval_[i];
      eigval_[i] = eigval_[idx];
      eigval_[idx] = tmp;
      if(reduce_){
        // Shift values from voutput(=vout+nevec*nelem) to vout and invert on the fly
        if(nevec_ > 0)
          for(int j = 0; j < nout; ++j)
            vout_[i * nout + j] = voutput[idx * nout + j];
      } else {
        // Invert within vout(=voutput)
        if(nevec_ > 0){
          for(int j = 0; j < nout; ++j) {
            tmp = vout_[i * nout + j];
            vout_[i * nout + j] = vout_[idx * nout + j];
            vout_[idx * nout + j] = tmp;
          }
        }
      }
    }
  }

  // Store average coordinates, eigenvectors, and eigenvalues in modesInfo
  modinfo_->SetAvg( minfo_->Ncols(), vect );
  if (nevec_ > 0)
    modinfo_->SetNvect( neval );
  else
    modinfo_->SetNvect( 0 );
  modinfo_->SetNvectElem( nout );
  modinfo_->SetFreq( eigval_ );
  modinfo_->SetEvec( vout_ );
  WorkspaceStored();

  return 0;
#endif
}

