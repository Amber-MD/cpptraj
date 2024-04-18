#include "Analysis_TICA.h"
#include "CoordCovarMatrix_Full.h"
#include "CoordCovarMatrix_Half.h"
#include "CpptrajStdio.h"
#include "DataSet_1D.h"
#include "DataSet_double.h"
#include "DataSet_MatrixDbl.h" // TODO remove?
#include "DataSet_Modes.h"
#include <cmath> //sqrt
#include <algorithm> // std::max

/** CONSTRUCTOR */
Analysis_TICA::Analysis_TICA() :
  TgtTraj_(0),
  lag_(0),
  useMass_(false),
  debugC0_(0),
  debugCT_(0),
  evectorScale_(NO_SCALING),
  ticaModes_(0),
  cumulativeVariance_(0)
{
  SetHidden(true);
}

// Analysis_TICA::Help()
void Analysis_TICA::Help() const {
  mprintf("\t{crdset <COORDS set name>|data <input set arg1> ...} [lag <time lag>]\n"
          "\t[mask <mask>] [mass] [map {kinetic|commute|none}]\n"
          "\t[name <output set name>] [out <file>] [cumvarout <file>]\n"
         );
          
}

// Analysis_TICA::Setup()
Analysis::RetType Analysis_TICA::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  // Attempt to get coords dataset or datasets from datasetlist
  TgtTraj_ = 0;
  sets_.clear();
  std::string setname = analyzeArgs.GetStringKey("crdset");
  std::string dataarg = analyzeArgs.GetStringKey("data");
  if (!setname.empty() && !dataarg.empty()) {
    mprinterr("Error: Specify either 'crdset' or 'data', not both.\n");
    return Analysis::ERR;
  }
  if (!setname.empty()) {
    TgtTraj_ = (DataSet_Coords*)setup.DSL().FindCoordsSet( setname );
    if (TgtTraj_ == 0) {
      mprinterr("Error: Could not locate COORDS set corresponding to %s\n",
                setname.c_str());
      Help();
      return Analysis::ERR;
    }
  } else if (!dataarg.empty()) {
    while (!dataarg.empty()) {
      if (sets_.AppendSetsFromArgs( ArgList(dataarg), setup.DSL() )) {
        mprinterr("Error: Could not add data sets using argument '%s'\n", dataarg.c_str());
        return Analysis::ERR;
      }
      dataarg = analyzeArgs.GetStringKey("data");
    }
  } else {
    mprinterr("Error: Must specify either 'crdset' or 'data'.\n");
    return Analysis::ERR;
  }
  // Other keywords
  std::string maparg = analyzeArgs.GetStringKey("map");
  if (!maparg.empty()) {
    if (maparg == "kinetic")
      evectorScale_ = KINETIC_MAP;
    else if (maparg == "commute")
      evectorScale_ = COMMUTE_MAP;
    else if (maparg == "none")
      evectorScale_ = NO_SCALING;
    else {
      mprinterr("Error: Unrecognized keyword for 'map' argument: %s\n", maparg.c_str());
      return Analysis::ERR;
    }
  } else
    evectorScale_ = KINETIC_MAP;
  lag_ = analyzeArgs.getKeyInt("lag", 1);
  std::string maskstr = analyzeArgs.GetStringKey("mask");
  if (mask1_.SetMaskString( maskstr )) {
    mprinterr("Error: Could not set atom mask string '%s'\n", maskstr.c_str());
    return Analysis::ERR;
  }
  maskstr = analyzeArgs.GetStringKey("mask2");
  if (!maskstr.empty()) {
    mprintf("DEBUG: Second mask detected.\n");
    if (mask2_.SetMaskString( maskstr )) {
      mprinterr("Error: Could not set second atom mask string '%s'\n", maskstr.c_str());
      return Analysis::ERR;
    }
  }

  useMass_ = analyzeArgs.hasKey("mass");

  // Allocate data sets/data file
  std::string dsname = analyzeArgs.GetStringKey("name");
  ticaModes_ = (DataSet_Modes*)setup.DSL().AddSet( DataSet::MODES, dsname, "TICA" );
  if (ticaModes_ == 0) {
    mprinterr("Error: Could not allocate output TICA modes.\n");
    return Analysis::ERR;
  }
  DataFile* outfile = setup.DFL().AddDataFile( analyzeArgs.GetStringKey("out"), analyzeArgs);
  if (outfile != 0)
    outfile->AddDataSet( ticaModes_ );

  cumulativeVariance_ = (DataSet_1D*)setup.DSL().AddSet(DataSet::DOUBLE, MetaData(ticaModes_->Meta().Name(), "cumvar"));
  if (cumulativeVariance_ == 0) {
    mprinterr("Error: Could not allocate output cumulative variance.\n");
    return Analysis::ERR;
  }
  // Change default precision
  cumulativeVariance_->SetupFormat().SetFormatWidthPrecision(12,8);
  DataFile* cvoutfile = setup.DFL().AddDataFile( analyzeArgs.GetStringKey("cumvarout"), analyzeArgs);
  if (cvoutfile != 0)
    cvoutfile->AddDataSet( cumulativeVariance_ );

  debugC0_ = setup.DFL().AddCpptrajFile(analyzeArgs.GetStringKey("debugc0"), "TICA C0 debug",
                                                                   DataFileList::TEXT, true);
  if (debugC0_ == 0) {
    mprinterr("Error: Could not open C0 debug file.\n");
    return Analysis::ERR;
  }
  debugCT_ = setup.DFL().AddCpptrajFile(analyzeArgs.GetStringKey("debugct"), "TICA CT debug",
                                                                   DataFileList::TEXT, true);
  if (debugCT_ == 0) {
    mprinterr("Error: Could not open CT debug file.\n");
    return Analysis::ERR;
  }

  // Print analysis info
  mprintf("    TICA: Time independent correlation analysis.\n");
  if (TgtTraj_ != 0) {
    mprintf("\tUsing coordinates from set '%s'\n", TgtTraj_->legend());
    mprintf("\tUsing atoms selected by mask '%s'\n", mask1_.MaskString());
    if (useMass_)
      mprintf("\tMass-weighted.\n");
    else
      mprintf("\tNot mass-weighted.\n");
  }
  if (!sets_.empty()) {
    mprintf("\tUsing %zu data sets:", sets_.size());
    for (Array1D::const_iterator it = sets_.begin(); it != sets_.end(); ++it)
      mprintf(" %s", (*it)->legend());
    mprintf("\n");
  }
  mprintf("\tTime lag: %i frames.\n", lag_);
  if (debugC0_ != 0)
    mprintf("\tDebug C0 output to %s\n", debugC0_->Filename().full());
  if (debugCT_ != 0)
    mprintf("\tDebug CT output to %s\n", debugCT_->Filename().full());
  switch (evectorScale_) {
    case NO_SCALING  : mprintf("\tNot scaling eigenvectors.\n"); break;
    case KINETIC_MAP : mprintf("\tScaling eigenvectors by eigenvalues.\n"); break;
    case COMMUTE_MAP : mprintf("\tScaling eigenvectors by regularized time scales.\n"); break;
  }
  mprintf("\tTICA data set: %s\n", ticaModes_->legend());
  if (outfile != 0)
    mprintf("\tTICA data set written to '%s'\n", outfile->DataFilename().full());
  mprintf("\tCumulative variance set: %s\n", cumulativeVariance_->legend());
  if (cvoutfile != 0)
    mprintf("\tCumulative variance written to '%s'\n", cvoutfile->DataFilename().full());

  return Analysis::OK;
}

// Analysis_TICA::Analyze()
Analysis::RetType Analysis_TICA::Analyze() {
  if (calcMatrices()) return Analysis::ERR;

  if (evectorScale_ == KINETIC_MAP) {
    // Weight eigenvectors by eigenvalues
    for (int ii = 0; ii < ticaModes_->Nmodes(); ii++)
      ticaModes_->MultiplyEvecByFac( ii, ticaModes_->Eigenvalue(ii) );
  } else if (evectorScale_ == COMMUTE_MAP) {
    // Weight eigenvectors by regularized time scales
    mprinterr("Internal Error: Commute_map not yet implemented.\n");
    return Analysis::ERR;
  }
  // Calculate cumulative variance
  cumulativeVariance_->Allocate(DataSet::SizeArray(1, ticaModes_->Nmodes()));
  double cumulativeSum = 0;
  for (int ii = 0; ii < ticaModes_->Nmodes(); ii++) {
    double dval = fabs(ticaModes_->Eigenvalue(ii));
    dval *= dval;
    cumulativeSum += dval;
    cumulativeVariance_->Add(ii, &cumulativeSum);
  }
  // Normalize the cumulative variance
  DataSet_double& cvset = static_cast<DataSet_double&>( *cumulativeVariance_ );
  for (int ii = 0; ii < ticaModes_->Nmodes(); ii++)
    cvset[ii] /= cumulativeSum;

  //if (TgtTraj_ != 0)
  //  return analyze_crdset();
  //else if (!sets_.empty())
  //  return analyze_datasets();

  return Analysis::OK;
}

// -----------------------------------------------------------------------------
/// FOR DEBUG
static inline void printDarray(const char* desc, std::vector<double> const& arrayIn, const char* fmt)
{
  //static const char* fmt = "%15.8f";
  mprintf("DEBUG: %s:  [", desc);
  int col = 0;
  for (std::vector<double>::const_iterator it = arrayIn.begin(); it != arrayIn.end(); ++it)
  {
    mprintf(fmt, *it);
    if ( (it+1) == arrayIn.end() ) mprintf("]");
    col++;
    if (col==4) {
      mprintf("\n");
      col = 0;
    }
  }
  if (col != 0)
    mprintf("\n");
}

/// For debugging - print eigenvalues/eigenvectors to stdout
/*static void printEigen(DataSet_Modes const& C0_Modes, const char* desc) {
  // DEBUG - print eigenvalues
  std::vector<double> tmpevals;
  for (int ii = 0; ii < C0_Modes.Nmodes(); ii++)
    tmpevals.push_back( C0_Modes.Eigenvalue(ii) );
  printDarray(desc, tmpevals, "%16.8e");
  // DEBUG - print first 3 values of each eigenvector
  for (int ii = 0; ii < C0_Modes.Nmodes(); ii++) {
    const double* evec = C0_Modes.Eigenvector(ii);
    for (int jj = 0; jj < 3; jj++)
      mprintf("%12.8f", evec[jj]);
    mprintf("\n");
  }
}*/

/// For debugging - print eigenvalues to file
static void printEvals(DataSet_Modes const& modes, const char* fname) {
  CpptrajFile outfile;
  if (outfile.OpenWrite(fname)) {
    mprinterr("Internal Error: printEvals: could not open file %s\n", fname);
    return;
  }
  for (int ii = 0; ii < modes.Nmodes(); ii++)
    outfile.Printf("%12.8f\n", modes.Eigenvalue(ii));
  outfile.CloseFile();
}

/// For debugging - print eigenvectors to file
static void printEvecs(DataSet_Modes const& modes, const char* fname) {
  CpptrajFile outfile;
  if (outfile.OpenWrite(fname)) {
    mprinterr("Internal Error: printEvals: could not open file %s\n", fname);
    return;
  }
  for (int ii = 0; ii < modes.Nmodes(); ii++) {
    const double* evec = modes.Eigenvector(ii);
    for (int jj = 0; jj < modes.VectorSize(); jj++)
      outfile.Printf("%12.8f", evec[jj]);
    outfile.Printf("\n");
  }
  outfile.CloseFile();
}

/// For debugging. Print matrix to file
static void printMatrix(const char* fname, DataSet_2D& matR, ArgList& tmpArgs,
                        int width, int precision, TextFormat::FmtType fmt)
{
  if (width > 0) {
    matR.SetupFormat().SetFormatWidthPrecision( width, precision );
    matR.SetupFormat().SetFormatType( fmt );
  }
  DataFile outfile8;
  outfile8.SetupDatafile(fname, tmpArgs, 0);
  outfile8.AddDataSet( &matR );
  outfile8.WriteDataOut();
  tmpArgs.SetAllUnmarked();
}

/// For debugging. Print matrix to file
static void printMatrix(const char* fname, DataSet_2D& matR, ArgList& tmpArgs)
{
  printMatrix(fname, matR, tmpArgs, -1, -1, TextFormat::DOUBLE);
}

// -----------------------------------------------------------------------------
/** Calculate combined effective sums for tau=0 and tau=lag for 1D data sets. */
void Analysis_TICA::calc_sums_from1Dsets(Darray& means, unsigned int end1) const {
  means.assign( sets_.size(), 0 );
  for (unsigned int idx = 0; idx != sets_.size(); ++idx) {
    DataSet_1D* seti = sets_[idx];
    //unsigned int end1 = seti->Size() - lag_;
    unsigned int k2 = lag_;
    for (unsigned int k1 = 0; k1 < end1; k1++, k2++) {
      means[idx] += (seti->Dval(k1) + seti->Dval(k2));
    }
  }
}

/** Calculate combined effective sums for tau=0 and tau=lag for periodic 1D data sets. */
void Analysis_TICA::calc_sums_fromPeriodicSets(Darray& means, unsigned int end1) const {
  means.assign( sets_.size() * 2, 0 );
  unsigned int jdx = 0;
  for (unsigned int idx = 0; idx != sets_.size(); ++idx, jdx += 2) {
    DataSet_1D* seti = sets_[idx];
    //unsigned int end1 = seti->Size() - lag_;
    unsigned int k2 = lag_;
    for (unsigned int k1 = 0; k1 < end1; k1++, k2++) {
      double theta1 = seti->Dval(k1) * Constants::DEGRAD;
      double theta2 = seti->Dval(k2) * Constants::DEGRAD;
      double x1 = cos( theta1 );
      double y1 = sin( theta1 );
      double x2 = cos( theta2 );
      double y2 = sin( theta2 );
      means[jdx  ] += (x1 + x2);
      means[jdx+1] += (y1 + y2);
    }
  }
}

/** Calculate combined effective sums for tau=0 and tau=lag for COORDS set. */
void Analysis_TICA::calc_sums_fromCoordsSet(Darray& means, unsigned int end1) const {
  means.assign( mask1_.Nselected() * 3, 0 );
  unsigned int k2 = lag;
  for (unsigned int k1 = 0; k1 < end1; k1++, k2++) {

/** Calculate total weight */
double Analysis_TICA::calc_total_weight(Darray const& weights, unsigned int end1) {
  // Total weight is times 2 because symmetric
  double total_weight = 0;
  if (!weights.empty()) {
    for (Darray::const_iterator it = weights.begin(); it != weights.end(); ++it)
      total_weight += *it;
    total_weight *= 2;
  } else
    total_weight = end1 * 2;
  return total_weight;
}

/** Create XXYY and XYYX matrices. */
void Analysis_TICA::create_matrices_from1Dsets( DataSet_2D* matXXYY,
                                                DataSet_2D* matXYYX,
                                                Darray const& means,
                                                unsigned int end1 )
const
{
  unsigned int Ncols = sets_.size(); // FIXME
  unsigned int Nrows = Ncols;
  unsigned int xxyyIdx = 0;
  unsigned int xyyxIdx = 0;

  matXXYY->AllocateHalf( Ncols );
  //matXYYX->Allocate2D( Ncols, Nrows );
  matXYYX->AllocateHalf( Ncols );

  for (unsigned int row = 0; row < Nrows; row++) {
    
    //for (unsigned int col = 0; col < Ncols; col++) {
    for (unsigned int col = row; col < Ncols; col++) {

      // XYYX
      DataSet_1D* seti = sets_[row];
      double offi = means[row];
      DataSet_1D* setj = sets_[col];
      double offj = means[col];
      //unsigned int end1 = seti->Size() - lag_;
      double sum = 0;
      double sumxx = 0;
      unsigned int k2 = lag_;
      for (unsigned int k1 = 0; k1 < end1; k1++, k2++) {
        double dvali1 = seti->Dval(k1) - offi;
        double dvalj2 = setj->Dval(k2) - offj;
        double dvali2 = seti->Dval(k2) - offi;
        double dvalj1 = setj->Dval(k1) - offj;
        //if (row == 0 && col == 0) mprintf("DBG1: %u %u %u %g %g\n", row, col, k1, dvali, dvalj);
        sum += (dvali1 * dvalj2);
        sum += (dvali2 * dvalj1);
        // XXYY
        sumxx += (dvali1 * dvalj1);
        sumxx += (dvalj2 * dvali2);
      }
      matXYYX->SetElement(xyyxIdx++, sum);
      matXXYY->SetElement(xxyyIdx++, sumxx);
      // XXYY
/*      if ( col >= row ) {
        double sumxx = 0;
        unsigned int k2 = lag_;
        for (unsigned int k1 = 0; k1 < end1; k1++, k2++) {
          double dvali1 = seti->Dval(k1) - offi;
          double dvalj2 = setj->Dval(k2) - offj;
          double dvali2 = seti->Dval(k2) - offi;
          double dvalj1 = setj->Dval(k1) - offj;
          sumxx += (dvali1 * dvalj1);
          sumxx += (dvalj2 * dvali2);
        }
        matXXYY->SetElement(xxyyIdx++, sumxx);
      }*/
    }
  }
}

// -----------------------------------------------------------------------------
/** Calculate XXYY and XYYX matrices. */
int Analysis_TICA::calcMatrices() const {
  // Check that sets have same size
  unsigned int maxFrames = 0;
  if (!sets_.empty()) {
    maxFrames = sets_.Array().front()->Size();
    for (DSarray::const_iterator it = sets_.begin(); it != sets_.end(); ++it)
    {
      if ((*it)->Size() != maxFrames) {
        mprinterr("Error: Set '%s' does not have same size (%zu) as first set (%u)\n",
                  (*it)->legend(), (*it)->Size(), maxFrames);
        return 1;
      }
    }
  } else {
    mprinterr("Internal Error: Finish calcMatrices().\n");
    return 1;
  }
  // Calculate start and end times for C0 and CT
  if ( (unsigned int)lag_ >= maxFrames ) {
    mprinterr("Error: lag %i >= max frames %u\n", lag_, maxFrames);
    return 1;
  }
  unsigned int c0end = maxFrames - (unsigned int)lag_;
  mprintf("DEBUG: C0 start = %u end = %u\n", 0, c0end);

  // Calculate means 
  Darray means;
  calc_sums_from1Dsets( means, c0end );
  double total_weight = calc_total_weight( Darray(), c0end );
  for (Darray::iterator it = means.begin(); it != means.end(); ++it)
    *it /= total_weight;
  CpptrajFile meanout; // DEBUG
  if (meanout.OpenWrite("test.mean.dat")) {
    mprinterr("Error: Could not open test.mean.dat\n");
    return 1;
  }
  for (Darray::const_iterator it = means.begin(); it != means.end(); ++it) {
    meanout.Printf("%12.6f\n", *it); // DEBUG
  }
  meanout.CloseFile();

  // Calculate matrices
  DataSet_MatrixDbl matXXYY, matXYYX;
  create_matrices_from1Dsets(static_cast<DataSet_2D*>(&matXXYY), static_cast<DataSet_2D*>(&matXYYX), means, c0end);
  ArgList tmpArgs("square2d noheader");
  printMatrix("test.xxyy.dat", matXXYY, tmpArgs, 12, 6, TextFormat::DOUBLE);
  printMatrix("test.xyyx.dat", matXYYX, tmpArgs, 12, 6, TextFormat::DOUBLE);
  // Normalize
  matXXYY.Normalize( 1.0 / total_weight );
  matXYYX.Normalize( 1.0 / total_weight );
  printMatrix("test.xxyy.norm.dat", matXXYY, tmpArgs);
  printMatrix("test.xyyx.norm.dat", matXYYX, tmpArgs);

  if (calculateTICA( means, matXXYY, matXYYX)) {
    mprinterr("Error: Caclulation of TICA modes failed.\n");
    return 1;
  }

  return 0;
}

// -----------------------------------------------------------------------------
int Analysis_TICA::calculateTICA(Darray const& meanX, DataSet_2D const& CXXYY, DataSet_2D const& Cxy) const {
  ArgList tmpArgs("square2d noheader");
  // ---------------------------------------------
  // Get C0 eigenvectors/eigenvalues
  DataSet_Modes C0_Modes;
  C0_Modes.SetAvgCoords( meanX );
  // Want all eigenvectors
  if (C0_Modes.CalcEigen( CXXYY, CXXYY.Ncols() )) {
    mprinterr("Error: Could not calculate eigenvectors and eigenvales for C0 matrix.\n");
    return 1;
  }
  // Eigenvec/vals should come out of CalcEigen sorted from largest Eigenval to smallest.
  // Determine cutoff; C0 is symmetric positive-definite so select
  // cutoff so that negative eigenvalues vanish.
  double min_eval = C0_Modes.Eigenvalue(C0_Modes.Nmodes()-1);
/*  double min_eval = C0_Modes.Eigenvalue(0);
  for (int ii = 1; ii < C0_Modes.Nmodes(); ii++) {
    if (C0_Modes.EigenValue(ii) < min_eval)
      min_eval = C0_Modes.Eigenvalue(ii);
  }*/
  mprintf("Min C0 eigenvalue= %g\n", min_eval);
  double epsilon = 1e-6; // FIXME TODO make user-definable
  if (min_eval < 0) {
    epsilon = std::max(epsilon, -min_eval + 1e-16);
  }
  // Get the absolute value of each eigenvector
  Darray abs_c0_evals;
  abs_c0_evals.reserve( C0_Modes.Nmodes() );
  for (int ii = 0; ii < C0_Modes.Nmodes(); ii++)
    abs_c0_evals.push_back( fabs( C0_Modes.Eigenvalue(ii) ) );
  // Find the index of the first abs(eigenvalue) smaller than epsilon
  int idx_first_smaller = 0;
  for (; idx_first_smaller < C0_Modes.Nmodes(); idx_first_smaller++) {
    if ( abs_c0_evals[idx_first_smaller] < epsilon ) {
      mprintf("%g < %g\n", abs_c0_evals[idx_first_smaller], epsilon );
      break;
    }
  }
  mprintf("DEBUG: Index of eigenvalue smaller than %g %i\n", epsilon, idx_first_smaller);
  C0_Modes.ResizeModes( idx_first_smaller );
  // Enforce canonical eigenvector signs
  C0_Modes.SetCanonicalEvecSigns();
/*  for (int ii = 0; ii < C0_Modes.Nmodes(); ii++) {
    // Find the maximum absolute value of the eigenvector
    double abs_maxval = 0;
    int abs_maxidx = 0;
    const double* evec = C0_Modes.Eigenvector(ii);
    for (int jj = 0; jj < C0_Modes.VectorSize(); jj++) {
      double dval = fabs( evec[jj] );
      if ( dval > abs_maxval ) {
        abs_maxval = dval;
        abs_maxidx = jj;
      }
    }
    mprintf("argmax %i %i %g\n", ii, abs_maxidx, evec[abs_maxidx]);
    double sign;
    if ( evec[abs_maxidx] < 0 )
      sign = -1.0;
    else
      sign =  1.0;
    // Multiply all elements of the eigenvector by the sign of the max abs element
    C0_Modes.MultiplyEvecByFac( ii, sign );
  }*/
  
  // DEBUG - print eigenvalues
  //printEigen( C0_Modes, "C0evals" );
  printEvals(C0_Modes, "sm.dat");
  printEvecs(C0_Modes, "Vm.dat");
/*  Darray tmpevals;
  for (int ii = 0; ii < C0_Modes.Nmodes(); ii++)
    tmpevals.push_back( C0_Modes.Eigenvalue(ii) );
  printDarray("C0evals", tmpevals, "%16.8e");
  // DEBUG - print first 3 values of each eigenvector
  for (int ii = 0; ii < C0_Modes.Nmodes(); ii++) {
    const double* evec = C0_Modes.Eigenvector(ii);
    for (int jj = 0; jj < 3; jj++)
      mprintf("%12.8f", evec[jj]);
    mprintf("\n");
  }*/

  // Create matrix Ltrans, where rows are eigenvectors of C0 times eigenvalues of C0
  DataSet_MatrixDbl matLtrans;
  matLtrans.Allocate2D( C0_Modes.VectorSize(), C0_Modes.Nmodes() );
  unsigned int idx = 0;
  for (int ii = 0; ii < C0_Modes.Nmodes(); ii++) {
    double fac = 1.0 / sqrt(C0_Modes.Eigenvalue(ii));
    const double* evec = C0_Modes.Eigenvector(ii);
    for (int jj = 0; jj < C0_Modes.VectorSize(); ++jj)
      matLtrans.SetElement(idx++, evec[jj] * fac);
  }
  // DEBUG - write unnormalized matrix
  printMatrix("matLtrans.dat", matLtrans, tmpArgs, 12, 8, TextFormat::DOUBLE);

  // L is transposed already (eigenvectors are in rows)
  // Calculate L^T * Ct
  DataSet_MatrixDbl matLtransCt;
  DataSet_2D::RetType ret = matLtransCt.Multiply(matLtrans, Cxy);
  if (ret != DataSet_2D::OK) {
    mprinterr("Error: Could not multiply L^T x Ct\n");
    return 1;
  }
  // DEBUG - write unnormalized matrix
  printMatrix("matLCt.dat", matLtransCt, tmpArgs, 15, 8, TextFormat::SCIENTIFIC);

  // Calculate (L^T * Ct) * L
  // Need to use transpose of Ltrans to get L
  DataSet_MatrixDbl Ct_trans;
  ret = Ct_trans.Multiply_M2transpose(matLtransCt, matLtrans);
  if (ret != DataSet_2D::OK) {
    mprinterr("Error: Could not multiply (L^T x Ct) x L\n");
    return 1;
  }
  // DEBUG - write unnormalized matrix
  printMatrix("Ct_trans.dat", Ct_trans, tmpArgs, 12, 8, TextFormat::DOUBLE);

  DataSet_Modes Ct_Modes;
  Ct_Modes.SetAvgCoords( meanX );
  // Want all eigenvectors
  if (Cxy.IsSymmetric()) {
    mprintf("\tCt is symmetric.\n");
    if (Ct_Modes.CalcEigen( Ct_trans, Ct_trans.Ncols() )) {
      mprinterr("Error: Could not calculate eigenvectors and eigenvalues for Ct matrix.\n");
      return 1;
    }
  } else {
    mprintf("\tCt is not symmetric.\n");
    if (Ct_Modes.CalcEigen_General( Ct_trans )) {
      mprinterr("Error: Could not calculate eigenvectors and eigenvalues for Ct matrix (non-symmetric).\n");
      return 1;
    }
  }
  // Sort by abs(eigenvalue)
  Ct_Modes.SortByAbsEigenvalue();
  //printEigen( Ct_Modes, "Ctevals" );
  printEvals(Ct_Modes, "ctvals.dat");
  printEvecs(Ct_Modes, "ctvecs.dat");

  // Calculate R^T * L^T to get (LR)^T
  DataSet_MatrixDbl ctm; // FIXME DataSet_Modes eigenvectors should be in a DataSet_MatrixDbl?
  ctm.Allocate2D( Ct_Modes.VectorSize(), Ct_Modes.Nmodes() );
  std::copy( Ct_Modes.Eigenvectors(),
             Ct_Modes.Eigenvectors() + (Ct_Modes.Nmodes() * Ct_Modes.VectorSize()),
             (double*)ctm.MatrixPtr() );
  DataSet_MatrixDbl matRt;
  matRt.Multiply( ctm, matLtrans );
  printMatrix("matRt.dat", matRt, tmpArgs, 12, 8, TextFormat::DOUBLE);
  //DataSet_Modes matRmodes;
  ticaModes_->SetAvgCoords( Ct_Modes.AvgCrd() );
  ticaModes_->SetModes(false, Ct_Modes.Nmodes(), Ct_Modes.VectorSize(),
                       Ct_Modes.EigenvaluePtr(), (const double*)matRt.MatrixPtr());
  ticaModes_->SetCanonicalEvecSigns();
  //DataFile outfile9;
  //outfile9.SetupDatafile("matRmodes.dat", 0);
  //outfile9.AddDataSet( ticaModes_ );
  //outfile9.WriteDataOut();
/*
  // The following is to test the math in the same exact order as L x Ct_Modes^T
  DataSet_MatrixDbl matL;
  ret = matL.TransposeOf( matLtrans );
  if (ret != DataSet_2D::OK) {
    mprinterr("Error: Could not get transpose of L^T\n");
    if ( ret == DataSet_2D::ALLOC_ERR) mprinterr("Error: Allocation error.\n");
    return 1;
  }
  // DEBUG - write unnormalized matrix
  printMatrix("matL.dat", matL, tmpArgs, 12, 8, TextFormat::DOUBLE);

  // Put eigenvectors in columns
  DataSet_MatrixDbl R_trans;
  if (R_trans.Allocate2D( Ct_Modes.Nmodes(), Ct_Modes.VectorSize() )) {
    mprinterr("Error: Could not allocate matrix for Ct eigenvectors.\n");
    return 1;
  }
  idx = 0;
  for (unsigned int row = 0; row < R_trans.Nrows(); row++) {
    for (unsigned int col = 0; col < R_trans.Ncols(); col++) {
      const double* evec = Ct_Modes.Eigenvector(col);
      R_trans[idx++] = evec[row];
    }
  }
  // DEBUG - write unnormalized matrix
  printMatrix("R_trans.dat", R_trans, tmpArgs, 12, 8, TextFormat::DOUBLE);

  // Calculate L * R^T
  // TODO calculate R^T * L^T to get (LR)^T instead?
  DataSet_MatrixDbl matR;
  matR.Multiply( matL, R_trans );

  // Change eigenvector (in columns) signs
  for (unsigned int col = 0; col < matR.Ncols(); col++) {
    // Find the maximum absolute value of the eigenvector (in column)
    double abs_maxval = 0;
    unsigned int abs_maxidx = 0;
    for (unsigned int row = 0; row < matR.Nrows(); row++) {
      double dval = fabs( matR.GetElement( col, row ) );
      if (dval > abs_maxval) {
        abs_maxval = dval;
        abs_maxidx = row;
      }
    }
    double elt = matR.GetElement(col, abs_maxidx);
    mprintf("argmax %u %u %g\n", col, abs_maxidx, elt);
    double sign;
    if (elt < 0)
      sign = -1.0;
    else
      sign = 1.0;
    // Multiply all elements of eigenvector (in column) by sign of max abs element
    for (unsigned int row = 0; row < matR.Nrows(); row++) {
      double dval = matR.GetElement( col, row ) * sign;
      matR.SetElement( col, row, dval );
    }
  }
  // DEBUG - write unnormalized matrix
  printMatrix("matR.dat", matR, tmpArgs, 12, 8, TextFormat::DOUBLE);
*/
  return 0;
}

// =============================================================================
/// Calculate sum over each set TODO weights
static inline void calculate_sum(std::vector<double>& sumX,
                                 std::vector<DataSet_1D*> const& sets, 
                                 unsigned int nelt_,
                                 unsigned int startFrame, unsigned int endFrame)
{
  sumX.assign( sets.size() * nelt_, 0 );

  if (nelt_ == 2) {
    mprinterr("Internal Error: CoordCovarMatrix_Half::AddDataToMatrix_C0CT(): Not implemented.\n");
    return;
  } else if (nelt_ == 1) {
    for (unsigned int jdx = 0; jdx < sets.size(); jdx++) {
      for (unsigned int idx = startFrame; idx < endFrame; idx++) {
        sumX[jdx] += sets[jdx]->Dval(idx);
      }
    }
  }
}

/// Subtract the set mean from every element of the set
static inline void subtract_mean(DataSet_double& out,
                                 DataSet_1D const* in,
                                 unsigned int nelt_,
                                 double total_weight,
                                 double sumX,
                                 unsigned int startFrame, unsigned int endFrame)
{
//  DataSet_double& out = static_cast<DataSet_double&>( *outPtr );
  //out.clear();
  out.Allocate(DataSet::SizeArray(1, endFrame - startFrame));
  //out.reserve(endFrame - startFrame);
  double mean = sumX / total_weight;
  for (unsigned int idx = startFrame; idx != endFrame; idx++) {
    out.AddElement( in->Dval(idx) - mean );
    //out.push_back( in->Dval(idx) - mean );
  }
}

/// Multiply transpose of matrix by matrix TODO since the resulting matrix is supposed to be symmetric we can speed this up
static void matT_times_mat( DataSet_2D* out,
                            std::vector<DataSet_double> const& M1,
                            std::vector<DataSet_double> const& M2)
{
  unsigned int Nrows = M1.size();
  unsigned int Ncols = M2.size();
  out->Allocate2D( Ncols, Nrows );
  unsigned int idx = 0;
  for (unsigned int row = 0; row < Nrows; row++) {
    DataSet_double const& seti = M1[row];
    for (unsigned int col = 0; col < Ncols; col++) {
      DataSet_double const& setj = M2[col];
      // seti->Size() must equal setj->Size()
      double sum = 0;
      for (unsigned int k = 0; k < seti.Size(); k++) {
        //if (row == 0 && col == 0) mprintf("DBG0: %u %u %u %g %g\n", row, col, k, seti.Dval(k), setj.Dval(k));
        sum += (seti.Dval(k) * setj.Dval(k));
      }
      out->SetElement(idx++, sum);
      //mprintf(" %12.4f", sum);
    }
    //mprintf("\n");
  }
}

/// For eaach matrix, Multiply transpose of matrix by matrix, then sum
static void matT_times_mat_symmetric( DataSet_2D* out,
                                      std::vector<DataSet_double> const& M1,
                                      std::vector<DataSet_double> const& M2 )
{
  if (M1.size() != M2.size()) {
    mprinterr("Internal Error: matT_times_mat_symmetric: Different # of sets.\n");
    return;
  }
  out->AllocateHalf( M1.size() );
  unsigned int Nrows = M1.size();
  unsigned int Ncols = Nrows;
  unsigned int idx = 0;
  for (unsigned int row = 0; row < Nrows; row++) {
    DataSet_double const& seti1 = M1[row];
    DataSet_double const& seti2 = M2[row];
    for (unsigned int col = row; col < Ncols; col++) {
      DataSet_double const& setj1 = M1[col];
      DataSet_double const& setj2 = M2[col];
      double sum = 0;
      // ALL SETS MUST HAVE SAME SIZE
      unsigned int nframes = seti1.Size();
      for (unsigned int k = 0; k < nframes; k++) {
        sum += (seti1.Dval(k) * setj1.Dval(k));
        sum += (seti2.Dval(k) * setj2.Dval(k));
      }
      out->SetElement(idx++, sum);
      //mprintf(" %12.4f", sum);
    }
    //mprintf("\n");
  }
}


/** Calculate instantaneous covariance and lagged covariance arrays */
int Analysis_TICA::calculateCovariance_C0CT(DSarray const& sets)
const
{
  static unsigned int nelt_ = 1; // FIXME
  // Check that sets have same size
  unsigned int maxFrames = sets.front()->Size();
  for (DSarray::const_iterator it = sets.begin(); it != sets.end(); ++it)
  {
    if ((*it)->Size() != maxFrames) {
      mprinterr("Error: Set '%s' does not have same size (%zu) as first set (%u)\n",
                (*it)->legend(), (*it)->Size(), maxFrames);
      return 1;
    }
  }
  // Calculate start and end times for C0 and CT
  if ( (unsigned int)lag_ >= maxFrames ) {
    mprinterr("Error: lag %i >= max frames %u\n", lag_, maxFrames);
    return 1;
  }
  unsigned int c0end = maxFrames - (unsigned int)lag_;
  mprintf("DEBUG: C0 start = %u end = %u\n", 0, c0end);
  unsigned int ctstart = (unsigned int)lag_;
  mprintf("DEBUG: CT start = %u end = %u\n", ctstart, maxFrames);
  // Calculate sum over each set
  Darray sumX;
  calculate_sum(sumX, sets, nelt_, 0, c0end);
  printDarray("sx_raw", sumX, "%15.8f");

  Darray sumY;
  calculate_sum(sumY, sets, nelt_, ctstart, maxFrames);
  printDarray("sy_raw", sumY, "%15.8f");

  // Sanity check
  if (sumX.size() != sumY.size()) {
    mprinterr("Internal Error: Analysis_TICA::calculateCovariance_C0CT(): sumX size != sumY size\n");
    return 1;
  }

  // Calculate effective sum (symmetric)
  Darray sx;
  sx.reserve( sumX.size() );
  for (unsigned int jdx = 0; jdx != sumX.size(); jdx++)
    sx.push_back( sumX[jdx] + sumY[jdx] );

  // Total weight is times 2 because symmetric
  //double total_weight = 0;
  //for (Darray::const_iterator it = weights.begin(); it != weights.end(); ++it)
  //  total_weight += *it;
  //total_weight *= 2;
  double total_weight = c0end * 2;
  mprintf("DEBUG: Total weight= %f\n", total_weight);
  // DEBUG - print mean
  CpptrajFile meanout;
  if (meanout.OpenWrite("debug.mean.dat")) {
    mprinterr("Error: Could not open debug.mean.dat\n");
    return 1;
  }
  Darray meanX;
  meanX.reserve( sx.size() );
  for (Darray::const_iterator it = sx.begin(); it != sx.end(); ++it) {
    meanX.push_back( *it / total_weight);
    meanout.Printf("%12.6f\n", meanX.back()); // DEBUG
  }
  meanout.CloseFile();


  // Center TODO sx_centered and sy_centered may not need to be calced
  Darray sx_centered, sy_centered;
  sx_centered.reserve( sumX.size() );
  sy_centered.reserve( sumY.size() );
  for (unsigned int idx = 0; idx != sumX.size(); idx++)
  {
    sx_centered.push_back( sumX[idx] - (0.5 * sx[idx]) );
    sy_centered.push_back( sumY[idx] - (0.5 * sx[idx]) );
  }
  printDarray("sx_raw_centered", sx_centered, "%16.8e");
  printDarray("sy_raw_centered", sy_centered, "%16.8e");

  // Remove mean
  typedef std::vector<DataSet_double> DDArray;
  DDArray CenteredX(sets.size());
  DDArray CenteredY(sets.size());

  Darray tmpx, tmpy; // DEBUG FIXME
  // Because symmetric, sy = sx //TODO use meanX set
  Darray const& sy = sx;
  for (unsigned int jdx = 0; jdx != sets.size(); jdx++) {
    subtract_mean(CenteredX[jdx], sets[jdx], nelt_, total_weight, sx[jdx], 0, c0end);
    subtract_mean(CenteredY[jdx], sets[jdx], nelt_, total_weight, sy[jdx], ctstart, maxFrames);
    tmpx.push_back( CenteredX[jdx].Dval(0) ); // DEBUG
    tmpy.push_back( CenteredY[jdx].Dval(0) ); // DEBUG
  }
  printDarray("X0", tmpx, "%16.8e");
  printDarray("Y0", tmpy, "%16.8e");

  // ---------------------------------------------
  ArgList tmpArgs("square2d noheader");
  // Calculate Cxxyy
  DataSet_MatrixDbl CXXYY;// = (DataSet_2D*)new DataSet_MatrixDbl();
  matT_times_mat_symmetric(static_cast<DataSet_2D*>(&CXXYY), CenteredX, CenteredY);
  // DEBUG - write unnormalized matrix
  printMatrix("cxxyy.dat", CXXYY, tmpArgs, 12, 6, TextFormat::DOUBLE);

  // Calculate Cxyyx
  DataSet_MatrixDbl Cxy, Cyx;
  matT_times_mat(static_cast<DataSet_2D*>(&Cxy), CenteredX, CenteredY);
  matT_times_mat(static_cast<DataSet_2D*>(&Cyx), CenteredY, CenteredX);
  for (unsigned int idx = 0; idx != Cxy.Size(); idx++) {
    double sum = Cxy.GetElement(idx) + Cyx.GetElement(idx);
    Cxy.SetElement(idx, sum);
  }
  // DEBUG - write unnormalized matrix
  printMatrix("cxyyx.dat", Cxy, tmpArgs, 12, 6, TextFormat::DOUBLE);

  // Normalize
  CXXYY.Normalize( 1.0 / total_weight );
  Cxy.Normalize( 1.0 / total_weight );
  // DEBUG - write normalized matrices
  printMatrix("cxxyy.norm.dat", CXXYY, tmpArgs);
  printMatrix("cxyyx.norm.dat", Cxy, tmpArgs);

  // ---------------------------------------------
  if (calculateTICA( meanX, CXXYY, Cxy)) {
    mprinterr("Error: Caclulation of TICA modes failed.\n");
    return 1;
  }

  return 0;
}

/** Analyze multiple 1D data sets. */
Analysis::RetType Analysis_TICA::analyze_datasets() {
  if (calculateCovariance_C0CT( sets_.Array() ) || ticaModes_ == 0) {
    mprinterr("Error: Could not calculate C0/CT matrices.\n");
    return Analysis::ERR;
  }
  
/*  // Matrix - half
  CoordCovarMatrix_Half covarMatrix;
  if (covarMatrix.SetupMatrix( sets_.Array() )) {
    mprinterr("Error: Could not set up C0 matrix for data sets.\n");
    return Analysis::ERR;
  }
  covarMatrix.AddDataToMatrix( sets_.Array() );
  // Normalize
  if (covarMatrix.FinishMatrix()) {
    mprinterr("Error: Could not normalize coordinate covariance matrix for C0.\n");
    return Analysis::ERR;
  }
  // DEBUG PRINT
  covarMatrix.DebugPrint("C0", *debugC0_, " %12.6f");*/

  return Analysis::OK;
}

/** Analyze coordinates data set. */
Analysis::RetType Analysis_TICA::analyze_crdset() {
  unsigned int Nframes = TgtTraj_->Size();
  if (Nframes < 1) {
    mprinterr("Error: No frames to analyze.\n");
    return Analysis::ERR;
  }
  // Evaluate mask
  if ( TgtTraj_->Top().SetupIntegerMask( mask1_ ) ) {
    mprinterr("Error: Could not evaluate atom mask '%s'\n", mask1_.MaskString());
    return Analysis::ERR;
  }
  mask1_.MaskInfo();
  if (mask1_.None()) {
    mprinterr("Error: No atoms selected by mask '%s'\n", mask1_.MaskString());
    return Analysis::ERR;
  }
  // DEBUG
  if (mask2_.MaskStringSet()) {
    if ( TgtTraj_->Top().SetupIntegerMask( mask2_ ) ) {
      mprinterr("Error: Could not evaluate second atom mask '%s'\n", mask2_.MaskString());
      return Analysis::ERR;
    }
    mask2_.MaskInfo();
    if (mask2_.None()) {
      mprinterr("Error: No atoms selected by second mask '%s'\n", mask2_.MaskString());
      return Analysis::ERR;
    }
  }
  // Allocate frames
  Frame coords0 = TgtTraj_->AllocateFrame();
  //coords0.SetupFrameFromMask( mask1_, TgtTraj_->Top().Atoms(), TgtTraj_->CoordsInfo() );
  //Frame coords1 = coords0;
  // Matrix - half
  CoordCovarMatrix_Half covarMatrix;
  if (covarMatrix.SetupMatrix(TgtTraj_->Top().Atoms(), mask1_, useMass_ )) {
    mprinterr("Error: Could not set up C0 matrix for '%s'.\n", TgtTraj_->legend());
    return Analysis::ERR;
  }
  // Loop over frames
  for (unsigned int frm0 = 0; frm0 < Nframes; frm0++) {
    //mprintf("DEBUG: Frame %i\n", frm0);
    //TgtTraj_->GetFrame(frm0, coords0, mask1_);
    TgtTraj_->GetFrame(frm0, coords0);
    // Covariance
    covarMatrix.AddFrameToMatrix( coords0, mask1_ );
  } // END loop over frames

  // Normalize
  if (covarMatrix.FinishMatrix()) {
    mprinterr("Error: Could not normalize coordinate covariance matrix for C0.\n");
    return Analysis::ERR;
  }

  // DEBUG PRINT
  covarMatrix.DebugPrint("C0", *debugC0_);

  if (mask2_.MaskStringSet()) {
    // DEBUG
    // Matrix - full
    CoordCovarMatrix_Full CT;
    CT.SetupMatrix(TgtTraj_->Top().Atoms(), mask1_,
                   TgtTraj_->Top().Atoms(), mask2_, useMass_);
    // Loop over frames
    for (unsigned int frm0 = 0; frm0 < Nframes; frm0++) {
      //mprintf("DEBUG: Frame %i\n", frm0);
      //TgtTraj_->GetFrame(frm0, coords0, mask1_);
      TgtTraj_->GetFrame(frm0, coords0);
      // Covariance
      CT.AddFrameToMatrix( coords0, mask1_, coords0, mask2_ );
    } // END loop over frames

    // Normalize
    if (CT.FinishMatrix()) {
      mprinterr("Error: Could not normalize coordinate covariance matrix for CT.\n");
      return Analysis::ERR;
    }

    // DEBUG PRINT
    CT.DebugPrint("CT", *debugCT_);
  }

  return Analysis::OK;
}
