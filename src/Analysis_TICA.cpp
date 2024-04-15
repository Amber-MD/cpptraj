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
  evectorScale_(NO_SCALING)
{
  SetHidden(true);
}

// Analysis_TICA::Help()
void Analysis_TICA::Help() const {
  mprintf("{crdset <set name>|data <set arg1> ...} [lag <time lag>]\n"
          "[mask <mask>] [mass] [map {kinetic|commute|none}]\n");
          
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
      if (sets_.AddSetsFromArgs( ArgList(dataarg), setup.DSL() )) {
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

  return Analysis::OK;
}

// Analysis_TICA::Analyze()
Analysis::RetType Analysis_TICA::Analyze() {
  if (TgtTraj_ != 0)
    return analyze_crdset();
  else if (!sets_.empty())
    return analyze_datasets();

  return Analysis::ERR;
}

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

/// For debugging - print eigenvalues/eigenvectors to stdout
static void printEigen(DataSet_Modes const& C0_Modes, const char* desc) {
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
}

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
  //for (unsigned int jdx = 0; jdx != sets.size(); jdx++) {
  //  CenteredX[jdx] = (DataSet_1D*)new DataSet_double();
  //  CenteredY[jdx] = (DataSet_1D*)new DataSet_double();
  //}
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
  CXXYY.SetupFormat().SetFormatWidthPrecision(12,6); // DEBUG
  CXXYY.SetupFormat().SetFormatType(TextFormat::DOUBLE); // DEBUG
  mprintf("CXXYY\n");
  matT_times_mat_symmetric(static_cast<DataSet_2D*>(&CXXYY), CenteredX, CenteredY);
  // DEBUG - write unnormalized matrix
  DataFile outfile1;
  outfile1.SetupDatafile("cxxyy.dat", tmpArgs, 0);
  outfile1.AddDataSet( &CXXYY );
  outfile1.WriteDataOut();
  tmpArgs.SetAllUnmarked();

  // Calculate Cxyyx
  DataSet_MatrixDbl Cxy, Cyx;
  Cxy.SetupFormat().SetFormatWidthPrecision(12,6); // DEBUG
  Cxy.SetupFormat().SetFormatType(TextFormat::DOUBLE); // DEBUG
  matT_times_mat(static_cast<DataSet_2D*>(&Cxy), CenteredX, CenteredY);
  matT_times_mat(static_cast<DataSet_2D*>(&Cyx), CenteredY, CenteredX);
  for (unsigned int idx = 0; idx != Cxy.Size(); idx++) {
    double sum = Cxy.GetElement(idx) + Cyx.GetElement(idx);
    Cxy.SetElement(idx, sum);
  }
  // DEBUG - write unnormalized matrix
  DataFile outfile2;
  outfile2.SetupDatafile("cxyyx.dat", tmpArgs, 0);
  outfile2.AddDataSet( &Cxy );
  outfile2.WriteDataOut();
  tmpArgs.SetAllUnmarked();

  // Normalize
  CXXYY.Normalize( 1.0 / total_weight );
  Cxy.Normalize( 1.0 / total_weight );
  // DEBUG - write normalized matrices
  outfile1.SetupDatafile("cxxyy.norm.dat", tmpArgs, 0);
  outfile1.WriteDataOut();
  tmpArgs.SetAllUnmarked();
  outfile2.SetupDatafile("cxyyx.norm.dat", tmpArgs, 0);
  outfile2.WriteDataOut();
  tmpArgs.SetAllUnmarked();

  // Free memory
  //for (unsigned int jdx = 0; jdx != sets.size(); jdx++) {
  //  delete CenteredX[jdx];
  //  delete CenteredY[jdx];
  //}

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
  for (int ii = 0; ii < C0_Modes.Nmodes(); ii++) {
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
  }
  
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
  matLtrans.SetupFormat().SetFormatWidthPrecision(12,8); // DEBUG
  matLtrans.SetupFormat().SetFormatType(TextFormat::DOUBLE); // DEBUG
  matLtrans.Allocate2D( C0_Modes.VectorSize(), C0_Modes.Nmodes() );
  unsigned int idx = 0;
  for (int ii = 0; ii < C0_Modes.Nmodes(); ii++) {
    double fac = 1.0 / sqrt(C0_Modes.Eigenvalue(ii));
    const double* evec = C0_Modes.Eigenvector(ii);
    for (int jj = 0; jj < C0_Modes.VectorSize(); ++jj)
      matLtrans.SetElement(idx++, evec[jj] * fac);
  }
  // DEBUG - write unnormalized matrix
  DataFile outfile3;
  outfile3.SetupDatafile("matLtrans.dat", tmpArgs, 0);
  outfile3.AddDataSet( &matLtrans );
  outfile3.WriteDataOut();
  tmpArgs.SetAllUnmarked();

  // L is transposed already (eigenvectors are in rows)
  // Calculate L^T * Ct
  DataSet_MatrixDbl matLtransCt;
  matLtransCt.SetupFormat().SetFormatWidthPrecision(15,8); // DEBUG
  matLtransCt.SetupFormat().SetFormatType(TextFormat::SCIENTIFIC); // DEBUG
  DataSet_2D::RetType ret = matLtransCt.Multiply(matLtrans, Cxy);
  if (ret != DataSet_2D::OK) {
    mprinterr("Error: Could not multiply L^T x Ct\n");
    return 1;
  }
  // DEBUG - write unnormalized matrix
  DataFile outfile4;
  outfile4.SetupDatafile("matLCt.dat", tmpArgs, 0);
  outfile4.AddDataSet( &matLtransCt );
  outfile4.WriteDataOut();
  tmpArgs.SetAllUnmarked();

  // Calculate (L^T * Ct) * L
  // Need to use transpose of Ltrans to get L
  DataSet_MatrixDbl Ct_trans;
  Ct_trans.SetupFormat().SetFormatWidthPrecision(12,8); // DEBUG
  Ct_trans.SetupFormat().SetFormatType(TextFormat::DOUBLE); // DEBUG
  ret = Ct_trans.Multiply_M2transpose(matLtransCt, matLtrans);
  if (ret != DataSet_2D::OK) {
    mprinterr("Error: Could not multiply (L^T x Ct) x L\n");
    return 1;
  }
  // DEBUG - write unnormalized matrix
  DataFile outfile5;
  outfile5.SetupDatafile("Ct_trans.dat", tmpArgs, 0);
  outfile5.AddDataSet( &Ct_trans );
  outfile5.WriteDataOut();
  tmpArgs.SetAllUnmarked();

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

  // The following is to test the math in the same exact order as L x Ct_Modes^T
  DataSet_MatrixDbl matL;
  ret = matL.TransposeOf( matLtrans );
  if (ret != DataSet_2D::OK) {
    mprinterr("Error: Could not get transpose of L^T\n");
    if ( ret == DataSet_2D::ALLOC_ERR) mprinterr("Error: Allocation error.\n");
    return 1;
  }
  // DEBUG - write unnormalized matrix
  DataFile outfile6;
  outfile6.SetupDatafile("matL.dat", tmpArgs, 0);
  outfile6.AddDataSet( &matL );
  outfile6.WriteDataOut();
  tmpArgs.SetAllUnmarked();
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
  R_trans.SetupFormat().SetFormatWidthPrecision(12,8); // DEBUG
  R_trans.SetupFormat().SetFormatType(TextFormat::DOUBLE); // DEBUG
  DataFile outfile7;
  outfile7.SetupDatafile("R_trans.dat", tmpArgs, 0);
  outfile7.AddDataSet( &R_trans );
  outfile7.WriteDataOut();
  tmpArgs.SetAllUnmarked();

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
  matR.SetupFormat().SetFormatWidthPrecision(12,8); // DEBUG
  matR.SetupFormat().SetFormatType(TextFormat::DOUBLE); // DEBUG
  DataFile outfile8;
  outfile8.SetupDatafile("matR.dat", tmpArgs, 0);
  outfile8.AddDataSet( &matR );
  outfile8.WriteDataOut();
  tmpArgs.SetAllUnmarked();

  return 0;
}

/** Analyze multiple 1D data sets. */
Analysis::RetType Analysis_TICA::analyze_datasets() {
  calculateCovariance_C0CT( sets_.Array() );

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
