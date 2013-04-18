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

void Analysis_Matrix::Help() {
  mprintf("\t<name> [out <filename>] [thermo outthermo <filename>]\n");
  mprintf("\t[vecs <#>] [name <modesname>] [reduce]\n");
  mprintf("\tDiagonalize given symmetric matrix to obtain eigenvectors\n");
  mprintf("\tand eigenvalues.\n");
}

// Analysis_Matrix::Setup()
Analysis::RetType Analysis_Matrix::Setup(ArgList& analyzeArgs, DataSetList* DSLin,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
#ifdef NO_MATHLIB
  mprinterr("Error: analyze matrix: Compiled without LAPACK routines.\n");
  return Analysis::ERR;
#else
  // Get matrix name
  std::string mname = analyzeArgs.GetStringNext();
  if (mname.empty()) {
    mprinterr("Error: analyze matrix: missing matrix name (first argument).\n");
    return Analysis::ERR;
  }
  // Find matrix in DataSetList.
  matrix_ = (DataSet_Matrix*)DSLin->FindSetOfType( mname, DataSet::MATRIX );
  if (matrix_ == 0) {
    mprinterr("Error: analyze matrix: Could not find matrix named %s\n",mname.c_str());
    return Analysis::ERR;
  }
  // Check that matrix is symmetric (half-matrix incl. diagonal).
  if (matrix_->Nrows() > 0) {
    mprinterr("Error: analyze matrix: Only works for symmetric matrices (i.e. no mask2)\n");
    return Analysis::ERR;
  }
  // Filenames
  outfilename_ = analyzeArgs.GetStringKey("out");
  // Thermo flag
  thermopt_ = analyzeArgs.hasKey("thermo");
  if (thermopt_)
    outthermo_ = analyzeArgs.GetStringKey("outthermo");
  if (thermopt_ && matrix_->Type()!=DataSet_Matrix::MWCOVAR) {
    mprinterr("Error: analyze matrix: parameter 'thermo' only works for\n");
    mprinterr("       mass-weighted covariance matrix ('mwcovar').\n");
    return Analysis::ERR;
  }
  // Number of eigenvectors; allow "0" only in case of 'thermo'
  nevec_ = analyzeArgs.getKeyInt("vecs",0);
  if (thermopt_) {
    if (nevec_ < 0) nevec_ = 0;
  } else if (nevec_ <= 0) {
    mprintf("Warning: # of eigenvectors specified is < 1 (%i) and 'thermo' not specified.\n",
            nevec_);
    mprintf("Warning: Specify # eigenvectors with 'vecs <#>'. Setting to 1.\n");
    nevec_ = 1;
  }
  // Reduce flag
  reduce_ = analyzeArgs.hasKey("reduce");
  if ( reduce_ && matrix_->Type() != DataSet_Matrix::MWCOVAR &&
                  matrix_->Type() != DataSet_Matrix::COVAR   &&
                  matrix_->Type() != DataSet_Matrix::DISTCOVAR  )
  {
    mprinterr("Error: analyze matrix: reduce not supported for %s\n", 
              DataSet_Matrix::MatrixTypeString[matrix_->Type()]);
    mprinterr("Error: reduce only works for covariance and distance covariance matrices.\n");
    return Analysis::ERR;
  }
  // Set up DataSet_Modes
  std::string modesname = analyzeArgs.GetStringKey("name");
  modes_ = (DataSet_Modes*)DSLin->AddSet( DataSet::MODES, modesname, "Modes" );
  if (modes_==0) return Analysis::ERR;
  // Output string for writing modes file.
  modes_->SetType( matrix_->Type() );

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
  return Analysis::OK;
#endif
}

Analysis::RetType Analysis_Matrix::Analyze() {
  modes_->SetAvgCoords( matrix_->VectSize(), matrix_->Vect() );
  // Calculate eigenvalues / eigenvectors
  if (modes_->CalcEigen( *matrix_, nevec_ )) return Analysis::ERR;
  if (matrix_->Type() == DataSet_Matrix::MWCOVAR) {
    if ( matrix_->Mass() == 0) {
      mprinterr("Error: MWCOVAR Matrix %s does not have mass info.\n", matrix_->Legend().c_str());
      return Analysis::ERR;
    }
    // Convert eigenvalues to cm^-1
    if (modes_->EigvalToFreq()) return Analysis::ERR;
    // Mass-wt eigenvectors
    if (modes_->MassWtEigvect( matrix_->Mass() )) return Analysis::ERR;
    // Calc thermo-chemistry if specified
    if (thermopt_) {
      // # of atoms - currently assuming COVAR (i.e. 3 matrix elts / coord)
      int natoms = matrix_->Nelts();
      CpptrajFile outfile;
      outfile.OpenWrite(outthermo_);
      thermo( outfile, natoms, modes_->Nmodes(), 1, matrix_->Vect(), 
              matrix_->Mass(), modes_->Eigenvalues(), 298.15, 1.0 );
      outfile.CloseFile();
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
  if (!outfilename_.empty())
    modes_->WriteToFile( outfilename_ ); 

  return Analysis::OK;
}
