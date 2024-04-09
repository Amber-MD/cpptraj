#include "Analysis_TICA.h"
#include "CoordCovarMatrix_Full.h"
#include "CoordCovarMatrix_Half.h"
#include "CpptrajStdio.h"
#include "DataSet_1D.h"
#include "DataSet_double.h"
#include "DataSet_MatrixDbl.h" // TODO remove?

/** CONSTRUCTOR */
Analysis_TICA::Analysis_TICA() :
  TgtTraj_(0),
  lag_(0),
  useMass_(false),
  debugC0_(0),
  debugCT_(0)
{
  SetHidden(true);
}

// Analysis_TICA::Help()
void Analysis_TICA::Help() const {
  mprintf("[crdset <set name>] [lag <time lag>] [mask <mask>] [mass]\n");
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
static inline void subtract_mean(DataSet_1D* outPtr,
                                 DataSet_1D const* in,
                                 unsigned int nelt_,
                                 double total_weight,
                                 double sumX,
                                 unsigned int startFrame, unsigned int endFrame)
{
  DataSet_double& out = static_cast<DataSet_double&>( *outPtr );
  //out.clear();
  out.Allocate(DataSet::SizeArray(1, endFrame - startFrame));
  //out.reserve(endFrame - startFrame);
  double mean = sumX / total_weight;
  for (unsigned int idx = startFrame; idx != endFrame; idx++) {
    out.AddElement( in->Dval(idx) - mean );
    //out.push_back( in->Dval(idx) - mean );
  }
}

/// Multiply transpose of matrix by matrix
static void matT_times_mat( DataSet_2D* out,
                            std::vector<DataSet_1D*> const& M1,
                            std::vector<DataSet_1D*> const& M2)
{
  unsigned int Nrows = M1.size();
  unsigned int Ncols = M2.size();
  out->Allocate2D( Ncols, Nrows );
  unsigned int idx = 0;
  for (unsigned int row = 0; row < Nrows; row++) {
    DataSet_1D* seti = M1[row];
    for (unsigned int col = 0; col < Ncols; col++) {
      DataSet_1D* setj = M2[col];
      // seti->Size() must equal setj->Size()
      double sum = 0;
      for (unsigned int k = 0; k < seti->Size(); k++) {
        sum += (seti->Dval(k) * setj->Dval(k));
      }
      out->SetElement(idx++, sum);
      //mprintf(" %12.4f", sum);
    }
    //mprintf("\n");
  }
}

/// For eaach matrix, Multiply transpose of matrix by matrix, then sum
static void matT_times_mat_symmetric( DataSet_2D* out,
                                      std::vector<DataSet_1D*> const& M1,
                                      std::vector<DataSet_1D*> const& M2 )
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
    DataSet_1D* seti1 = M1[row];
    DataSet_1D* seti2 = M2[row];
    for (unsigned int col = row; col < Ncols; col++) {
      DataSet_1D* setj1 = M1[col];
      DataSet_1D* setj2 = M2[col];
      double sum = 0;
      // ALL SETS MUST HAVE SAME SIZE
      unsigned int nframes = seti1->Size();
      for (unsigned int k = 0; k < nframes; k++) {
        sum += (seti1->Dval(k) * setj1->Dval(k));
        sum += (seti2->Dval(k) * setj2->Dval(k));
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

  // Center
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

  // Remove mean // FIXME do without new
  typedef std::vector<DataSet_1D*> DDArray;
  DDArray CenteredX(sets.size());
  DDArray CenteredY(sets.size());
  for (unsigned int jdx = 0; jdx != sets.size(); jdx++) {
    CenteredX[jdx] = (DataSet_1D*)new DataSet_double();
    CenteredY[jdx] = (DataSet_1D*)new DataSet_double();
  }
  Darray tmpx, tmpy; // DEBUG FIXME
  // Because symmetric, sy = sx
  Darray const& sy = sx;
  for (unsigned int jdx = 0; jdx != sets.size(); jdx++) {
    subtract_mean(CenteredX[jdx], sets[jdx], nelt_, total_weight, sx[jdx], 0, c0end);
    subtract_mean(CenteredY[jdx], sets[jdx], nelt_, total_weight, sy[jdx], ctstart, maxFrames);
    tmpx.push_back( CenteredX[jdx]->Dval(0) );
    tmpy.push_back( CenteredY[jdx]->Dval(0) );
    //tmpx.push_back( CenteredX[jdx][0] ); // DEBUG FIXME
    //tmpy.push_back( CenteredY[jdx][0] ); // DEBUG FIXME
  }
  printDarray("X0", tmpx, "%16.8e");
  printDarray("Y0", tmpy, "%16.8e");

  // Calculate Cxxyy
  DataSet_MatrixDbl CXXYY;// = (DataSet_2D*)new DataSet_MatrixDbl();
  mprintf("CXXYY\n");
  matT_times_mat_symmetric(static_cast<DataSet_2D*>(&CXXYY), CenteredX, CenteredY);
  DataFile outfile1;
  outfile1.SetupDatafile("cxxyy.dat", 0);
  outfile1.AddDataSet( &CXXYY );
  outfile1.WriteDataOut();

  // Calculate Cxyyx
  DataSet_MatrixDbl Cxy, Cyx;
  matT_times_mat(static_cast<DataSet_2D*>(&Cxy), CenteredX, CenteredY);
  matT_times_mat(static_cast<DataSet_2D*>(&Cyx), CenteredY, CenteredX);
  for (unsigned int idx = 0; idx != Cxy.Size(); idx++) {
    double sum = Cxy.GetElement(idx) + Cyx.GetElement(idx);
    Cxy.SetElement(idx, sum);
  }
  DataFile outfile2;
  outfile2.SetupDatafile("cxyyx.dat", 0);
  outfile2.AddDataSet( &Cxy );
  outfile2.WriteDataOut();

  //matT_times_mat(CenteredX, CenteredX);
  //mprintf("CYY\n");
  //matT_times_mat(CenteredY);
  //matT_times_mat(CenteredY, CenteredY);
/*  CoordCovarMatrix_Half Cxx;
  if (Cxx.SetupMatrix( CenteredX )) {
    mprinterr("Error: Could not set up Cxx matrix.\n");
    return 1;
  }
  Cxx.AddDataToMatrix( CenteredX );
  CpptrajFile cxxfile;
  cxxfile.OpenWrite("cxx.dat");
  Cxx.DebugPrint("Cxx", cxxfile, " %12.6f");
  CoordCovarMatrix_Half Cyy;
  if (Cyy.SetupMatrix( CenteredY )) {
    mprinterr("Error: Could not set up Cyy matrix.\n");
    return 1;
  }
  Cyy.AddDataToMatrix( CenteredY );
  CpptrajFile cyyfile;
  cyyfile.OpenWrite("cyy.dat");
  Cyy.DebugPrint("Cyy", cyyfile, " %12.6f");*/


  // Free memory
  for (unsigned int jdx = 0; jdx != sets.size(); jdx++) {
    delete CenteredX[jdx];
    delete CenteredY[jdx];
  }

  return 0;
}

/** Analyze multiple 1D data sets. */
Analysis::RetType Analysis_TICA::analyze_datasets() {
  calculateCovariance_C0CT( sets_.Array() );

  // Matrix - half
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
  covarMatrix.DebugPrint("C0", *debugC0_, " %12.6f");

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
