#include "Analysis_Matrix.h"
#include "CpptrajStdio.h"
#include "DataSet_MatrixDbl.h"

// CONSTRUCTOR
Analysis_Matrix::Analysis_Matrix() :
  matrix_(0),
  modes_(0),
  outthermo_(0),
  thermo_temp_(298.15),
  nevec_(0),
  thermopt_(false),
  reduce_(false),
  eigenvaluesOnly_(false),
  nmwizopt_(false),
  nmwizvecs_(0),
  nmwizfile_(0)
{}

void Analysis_Matrix::Help() const {
  mprintf("\t<name> [out <filename>] [thermo [outthermo <filename>] [temp <T>]]\n"
          "\t[vecs <#>] [name <modesname>] [reduce]\n"
          "\t[ nmwiz [nmwizvecs <n>] [nmwizfile <file>] %s\n"
          "\t  nmwizmask <mask> ]\n"
          "  Diagonalize given symmetric matrix to obtain eigenvectors\n"
          "  and eigenvalues.\n", DataSetList::TopArgs);
}

// Analysis_Matrix::Setup()
Analysis::RetType Analysis_Matrix::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
#ifdef NO_MATHLIB
  mprinterr("Error: Compiled without LAPACK routines.\n");
  return Analysis::ERR;
#else
  // Get matrix name
  std::string mname = analyzeArgs.GetStringNext();
  if (mname.empty()) {
    mprinterr("Error: Missing matrix name (first argument).\n");
    return Analysis::ERR;
  }
  // Find matrix in DataSetList.
  matrix_ = (DataSet_2D*)setup.DSL().FindSetOfType( mname, DataSet::MATRIX_DBL );
  if (matrix_ == 0)
    matrix_ = (DataSet_2D*)setup.DSL().FindSetOfType( mname, DataSet::MATRIX_FLT );
  if (matrix_ == 0) {
    mprinterr("Error: Could not find matrix named %s\n",mname.c_str());
    return Analysis::ERR;
  }
  // Check that matrix is symmetric (half-matrix incl. diagonal).
  if (matrix_->MatrixKind() != DataSet_2D::HALF) {
    mprinterr("Error: Only works for symmetric matrices (i.e. no mask2)\n");
    return Analysis::ERR;
  }
  
  // nmwiz flag
  nmwizopt_ = analyzeArgs.hasKey("nmwiz");
  if (nmwizopt_) { 
    nmwizvecs_ = analyzeArgs.getKeyInt("nmwizvecs", 20);
    if (nmwizvecs_ < 1) {
      mprinterr("Error: nmwizvecs must be >= 1\n");
      return Analysis::ERR;
    }
    nmwizfile_ = setup.DFL().AddCpptrajFile(analyzeArgs.GetStringKey("nmwizfile"), "NMwiz output",
                                       DataFileList::TEXT, true);
    Topology* parmIn = setup.DSL().GetTopology( analyzeArgs); // TODO: Include with matrix
    if (parmIn == 0) {
      mprinterr("Error: nmwiz: No topology specified.\n");
      return Analysis::ERR;
    }
    AtomMask nmwizMask( analyzeArgs.GetStringKey("nmwizmask") );
    if (parmIn->SetupIntegerMask( nmwizMask )) return Analysis::ERR;
    nmwizMask.MaskInfo();
    Topology* nparm = parmIn->partialModifyStateByMask( nmwizMask );
    if (nparm == 0) return Analysis::ERR;
    nmwizParm_ = *nparm;
    delete nparm;
    nmwizParm_.Brief("nmwiz topology");
  }
  
  // Filenames
  DataFile* outfile = setup.DFL().AddDataFile( analyzeArgs.GetStringKey("out"), analyzeArgs);
  // Thermo flag
  thermopt_ = analyzeArgs.hasKey("thermo");
  if (thermopt_) {
    outthermo_ = setup.DFL().AddCpptrajFile(analyzeArgs.GetStringKey("outthermo"), "'thermo' output",
                                       DataFileList::TEXT, true);
    if (outthermo_ == 0) return Analysis::ERR;
  }
  thermo_temp_ = analyzeArgs.getKeyDouble("temp", 298.15);
  if (thermopt_ && matrix_->Meta().ScalarType() != MetaData::MWCOVAR) {
    mprinterr("Error: Parameter 'thermo' only works for mass-weighted covariance matrix ('mwcovar').\n");
    return Analysis::ERR;
  }
  // Number of eigenvectors; allow "0" only in case of 'thermo'. -1 means 'All'
  nevec_ = analyzeArgs.getKeyInt("vecs",-1);
  if (nevec_ == 0 && !thermopt_) {
    mprintf("Warning: # of eigenvectors specified is 0 and 'thermo' not specified.\n");
    mprintf("Warning: Specify # eigenvectors with 'vecs <#>'. Setting to All.\n");
    nevec_ = -1;
  }
  // Reduce flag
  reduce_ = analyzeArgs.hasKey("reduce");
  // Set up DataSet_Modes. Set Modes DataSet type to be same as input matrix. 
  MetaData md( analyzeArgs.GetStringKey("name") );
  md.SetScalarType( matrix_->Meta().ScalarType() );
  modes_ = (DataSet_Modes*)setup.DSL().AddSet( DataSet::MODES, md, "Modes" );
  if (modes_==0) return Analysis::ERR;
  if (outfile != 0) outfile->AddDataSet( modes_ );

  // Print Status
  mprintf("    DIAGMATRIX: Diagonalizing matrix %s",matrix_->legend());
  if (outfile != 0)
    mprintf(" and writing modes to %s", outfile->DataFilename().full());
  if (nevec_ > 0)
    mprintf("\n\tCalculating %i eigenvectors.\n", nevec_);
  else if (nevec_ == 0)
    mprintf("\n\tNot calculating eigenvectors.\n");
  else
    mprintf("\n\tCalculating all eigenvectors.\n");
  if (thermopt_)
    mprintf("\tCalculating thermodynamic data at %.2f K, output to %s\n",
            thermo_temp_, outthermo_->Filename().full());
  if (nmwizopt_)
    mprintf("\tWriting %i modes to NMWiz file %s", nmwizvecs_, nmwizfile_->Filename().full());
  if (nevec_>0 && reduce_)
    mprintf("\tEigenvectors will be reduced\n");
  mprintf("\tStoring modes with name: %s\n", modes_->Meta().Name().c_str());
  return Analysis::OK;
#endif
}

// Analysis_Matrix::Analyze()
Analysis::RetType Analysis_Matrix::Analyze() {
  // Set the averaged coordinates and masses from matrix.
  if (modes_->SetAvgCoords( *matrix_ )) return Analysis::ERR;
  mprintf("\tEigenmode calculation for '%s'\n", matrix_->legend());
  // Check that the matrix was generated with enough snapshots.
  if (matrix_->Type() == DataSet::MATRIX_DBL) {
    DataSet_MatrixDbl const& Dmatrix = static_cast<DataSet_MatrixDbl const&>( *matrix_ );
    if (Dmatrix.Nsnapshots() < Dmatrix.Ncols())
      mprintf("Warning: In matrix '%s', # of frames %u is less than # of columns %zu.\n"
              "Warning: The max # of non-zero eigenvalues will be %u\n", Dmatrix.legend(),
              Dmatrix.Nsnapshots(), Dmatrix.Ncols(), Dmatrix.Nsnapshots());
  }
  // Calculate eigenvalues / eigenvectors
  if (modes_->CalcEigen( *matrix_, nevec_ )) return Analysis::ERR;
  // If mass-weighted covariance, mass-weight the resulting eigenvectors.
  if (matrix_->Meta().ScalarType() == MetaData::MWCOVAR) {
    mprintf("Info: Converting eigenvalues to cm^-1 and mass-weighting eigenvectors.\n");
    // Convert eigenvalues to cm^-1
    if (modes_->EigvalToFreq(thermo_temp_)) return Analysis::ERR;
    // Mass-wt eigenvectors
    if (modes_->MassWtEigvect()) return Analysis::ERR;
    // Calc thermo-chemistry if specified
    if (thermopt_)
      modes_->Thermo( *outthermo_, 1, thermo_temp_, 1.0 );
  }
  // Print nmwiz file if specified
  if (nmwizopt_) NMWizOutput(); 

  if (reduce_) {
    if (modes_->ReduceVectors()) return Analysis::ERR;
  }
  //modes_->PrintModes(); // DEBUG

  return Analysis::OK;
}

int Analysis_Matrix::NMWizOutput() const {
  // Check # vecs
  int nvecs;
  if (nmwizvecs_ <= modes_->Nmodes())
    nvecs = nmwizvecs_;
  else {
    mprintf("Warning: nmwizvecs > # eigenvectors, only writing %i vecs.\n",
            modes_->Nmodes());
    nvecs = modes_->Nmodes();
  }
  // Check # atoms
  if (nmwizParm_.Natom() * 3 != modes_->VectorSize()) {
    mprinterr("Error: nmwiz topology size %i does not match eigenvector size %i.\n",
              nmwizParm_.Natom() * 3, modes_->VectorSize());
    return 1;
  }
  
  nmwizfile_->Printf("nmwiz_load %s\n", nmwizfile_->Filename().full());

  nmwizfile_->Printf("name default_name\n");  //TODO: get from optionally provided pdb file

  nmwizfile_->Printf("atomnames ");
  for (Topology::atom_iterator atom = nmwizParm_.begin(); atom != nmwizParm_.end(); ++atom)
        nmwizfile_->Printf("%s ", atom->c_str());
  nmwizfile_->Printf("\n");

  nmwizfile_->Printf("resnames ");
  for (Topology::atom_iterator atom = nmwizParm_.begin(); atom != nmwizParm_.end(); ++atom)
    nmwizfile_->Printf("%s ", nmwizParm_.Res(atom->ResNum()).c_str());
  nmwizfile_->Printf("\n");

  nmwizfile_->Printf("resids ");
  for (Topology::atom_iterator atom = nmwizParm_.begin(); atom != nmwizParm_.end(); ++atom)
    nmwizfile_->Printf("%d ", atom->ResNum()+1);
  nmwizfile_->Printf("\n");

  nmwizfile_->Printf("chainids \n");    //TODO: get from optionally provided pdb file

  nmwizfile_->Printf("bfactors \n");    //TODO: get from optionally provided pdb file

  nmwizfile_->Printf("coordinates ");
  for (int i = 0; i < modes_->NavgCrd(); ++i)
          nmwizfile_->Printf("%8.3f ", modes_->AvgCrd()[i]);
  nmwizfile_->Printf("\n");

  for (int vec = 0; vec < nvecs; ++vec){
    nmwizfile_->Printf("mode %i %12.10f ", vec+1, 1/modes_->Eigenvalue(vec));
    const double* Vec = modes_->Eigenvector(vec);
    for (int i = 0 ; i < modes_->VectorSize(); ++i)
      nmwizfile_->Printf("%12.5f ", Vec[i]);
    nmwizfile_->Printf("\n");
  }
  return 0;
}
