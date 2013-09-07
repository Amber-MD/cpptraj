#include "Analysis_Matrix.h"
#include "CpptrajStdio.h"
#include "DataSet_MatrixDbl.h"

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
  matrix_ = (DataSet_2D*)DSLin->FindSetOfType( mname, DataSet::MATRIX_DBL );
  if (matrix_ == 0)
    matrix_ = (DataSet_2D*)DSLin->FindSetOfType( mname, DataSet::MATRIX_FLT );
  if (matrix_ == 0) {
    mprinterr("Error: analyze matrix: Could not find matrix named %s\n",mname.c_str());
    return Analysis::ERR;
  }
  // Check that matrix is symmetric (half-matrix incl. diagonal).
  if (matrix_->Kind() != DataSet_2D::HALF) {
    mprinterr("Error: analyze matrix: Only works for symmetric matrices (i.e. no mask2)\n");
    return Analysis::ERR;
  }
  // Filenames
  outfilename_ = analyzeArgs.GetStringKey("out");
  // Thermo flag
  thermopt_ = analyzeArgs.hasKey("thermo");
  if (thermopt_)
    outthermo_ = analyzeArgs.GetStringKey("outthermo");
  if (thermopt_ && matrix_->Type()!=DataSet_2D::MWCOVAR) {
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

// Analysis_Matrix::Analyze()
Analysis::RetType Analysis_Matrix::Analyze() {
  modes_->SetAvgCoords( *matrix_ );
  // Calculate eigenvalues / eigenvectors
  if (modes_->CalcEigen( *matrix_, nevec_ )) return Analysis::ERR;
  if (matrix_->Type() == DataSet_2D::MWCOVAR) {
    DataSet_MatrixDbl* Dmatrix = static_cast<DataSet_MatrixDbl*>( matrix_ );
    if ( Dmatrix->Mass().empty() ) {
      mprinterr("Error: MWCOVAR Matrix %s does not have mass info.\n", matrix_->Legend().c_str());
      return Analysis::ERR;
    }
    // Convert eigenvalues to cm^-1
    if (modes_->EigvalToFreq()) return Analysis::ERR;
    // Mass-wt eigenvectors
    if (modes_->MassWtEigvect( Dmatrix->Mass() )) return Analysis::ERR;
    // Calc thermo-chemistry if specified
    if (thermopt_) {
      CpptrajFile outfile;
      outfile.OpenWrite(outthermo_);
      modes_->Thermo( outfile, 1, 298.15, 1.0 );
      outfile.CloseFile();
    }
  }
  if (reduce_) {
    if (modes_->Reduce()) return Analysis::ERR;
  }
  //modes_->PrintModes(); // DEBUG
  if (!outfilename_.empty())
    modes_->WriteToFile( outfilename_ ); 

  return Analysis::OK;
}
