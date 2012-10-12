#include "Analysis_Matrix.h"
#include "CpptrajStdio.h"
#include "Thermo.h"

// CONSTRUCTOR
Analysis_Matrix::Analysis_Matrix() :
  matrix_(0),
  modes_(0),
  nevec_(0),
  thermopt_(false),
  reduce_(false),
  eigenvaluesOnly_(false)
{}

// Analysis_Matrix::Setup()
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
  matrix_ = (DataSet_Matrix*)DSLin->FindSetOfType( mname, DataSet::MATRIX );
  if (matrix_ == 0) {
    mprinterr("Error: analyze matrix: Could not find matrix named %s\n",mname.c_str());
    return 1;
  }
  // Check that matrix is symmetric (half-matrix incl. diagonal).
  if (matrix_->Nrows() > 0) {
    mprinterr("Error: analyze matrix: Only works for symmetric matrices (i.e. no mask2)\n");
    return 1;
  }
  // Filenames
  outfilename_ = analyzeArgs_.GetStringKey("out");
  // Thermo flag
  thermopt_ = analyzeArgs_.hasKey("thermo");
  if (thermopt_)
    outthermo_ = analyzeArgs_.GetStringKey("outthermo");
  if (thermopt_ && matrix_->Mass()==0) {
    mprinterr("Error: analyze matrix: parameter 'thermo' only works for\n");
    mprinterr("       mass-weighted covariance matrix ('mwcovar').\n");
    return 1;
  }
  // Number of eigenvectors; allow "0" only in case of 'thermo'
  nevec_ = analyzeArgs_.getKeyInt("vecs",0);
  if (thermopt_) {
    if (nevec_ < 0) nevec_ = 0;
  } else if (nevec_ <= 0) {
    mprintf("Warning: # of eigenvectors specified is < 1 (%i) and 'thermopt' not specified.\n",
            nevec_);
    mprintf("Warning: Specify # eigenvectors with 'vecs <#>'. Setting to 1.\n");
    nevec_ = 1;
  }
  // Reduce flag
  reduce_ = analyzeArgs_.hasKey("reduce");
  if ( reduce_ && matrix_->Type() != DataSet_Matrix::MWCOVAR &&
                  matrix_->Type() != DataSet_Matrix::COVAR   &&
                  matrix_->Type() != DataSet_Matrix::DISTCOVAR  )
  {
    mprinterr("Error: analyze matrix: reduce not supported for %s\n", 
              DataSet_Matrix::MatrixTypeString[matrix_->Type()]);
    mprinterr("Error: reduce only works for covariance and distance covariance matrices.\n");
    return 1;
  }
  // Set up DataSet_Modes
  std::string modesname = analyzeArgs_.GetStringKey("name");
  modes_ = (DataSet_Modes*)DSLin->AddSet( DataSet::MODES, modesname, "Modes" );
  if (modes_==0) return 1;
  // Output string for writing modes file.
  modes_->SetType( matrix_->Type() );
  modes_->SetAvgCoords( matrix_->VectSize(), matrix_->Vect() );

  // Print Status
  mprintf("    ANALYZE MATRIX: Analyzing matrix %s",matrix_->Legend().c_str());
  if (!outfilename_.empty())
    mprintf(" and dumping results to %s\n", outfilename_.c_str());
  else
    mprintf(" and printing results to STDOUT\n");
  mprintf("      Calculating %i eigenvectors", nevec_);
  if (thermopt_) {
    mprintf(" and thermodynamic data, output to");
    if (!outthermo_.empty())
      mprintf(" %s", outthermo_.c_str());
    else
      mprintf(" STDOUT");
  }
  mprintf("\n");
  if (nevec_>0 && reduce_)
    mprintf("      Eigenvectors will be reduced\n");
  if (!modesname.empty())
    mprintf("      Storing modes with name: %s\n",modesname.c_str());
  return 0;
#endif
}

int Analysis_Matrix::Analyze() {
  // Calculate eigenvalues / eigenvectors
  if (modes_->CalcEigen( *matrix_, nevec_ )) return 1;
  if (matrix_->Type() == DataSet_Matrix::MWCOVAR) {
    // Convert eigenvalues to cm^-1
    if (modes_->EigvalToFreq()) return 1;
    // Mass-wt eigenvectors
    if (modes_->MassWtEigvect( matrix_->Mass() )) return 1;
    // Calc thermo-chemistry if specified
    if (thermopt_) {
      // Temp storage for sorting eigenvalues in descending order
      double* eigvali = new double[ modes_->Nmodes() ];
      int j = modes_->Nmodes() - 1;
      for (int i = 0; i < modes_->Nmodes(); ++i)
        eigvali[j--] = modes_->Eigenvalue(i);
      // # of atoms - currently assuming COVAR (i.e. 3 matrix elts / coord)
      int natoms = matrix_->Nelts();
      CpptrajFile outfile;
      outfile.OpenWrite(outthermo_);
      thermo( outfile, natoms, modes_->Nmodes(), 1, matrix_->Vect(), 
              matrix_->Mass(), eigvali, 298.15, 1.0 );
      outfile.CloseFile();
      delete[] eigvali;
    }
  }
  if (reduce_) {
    if ( matrix_->Type() == DataSet_Matrix::COVAR ||
         matrix_->Type() == DataSet_Matrix::MWCOVAR )
      modes_->ReduceCovar();
    else if ( matrix_->Type() == DataSet_Matrix::DISTCOVAR )
      modes_->ReduceDistCovar( matrix_->Nelts() );
  }
  //modes_->PrintModes(); // DEBUG
  return 0;
}

void Analysis_Matrix::Print(DataFileList* DFLin) {
  if (!outfilename_.empty()) {
    CpptrajFile outfile;
    outfile.OpenWrite( outfilename_ );
    modes_->WriteToFile( outfile );
    outfile.CloseFile();
  }
}
