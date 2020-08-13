#include <cmath> // sqrt
#include "Analysis_Timecorr.h"
#include "CpptrajStdio.h"
#include "DataSet_double.h"
#include "DataSet_Vector.h"

/// Strings corresponding to modes, used in output.
const char* Analysis_Timecorr::ModeString[] = { 
  "Auto-correlation function", "Cross-correlation function" };

// CONSTRUCTOR
Analysis_Timecorr::Analysis_Timecorr() :
  tstep_(1.0),
  tcorr_(10000.0),
  order_(2),
  mode_(AUTOCORR),
  dplr_(false),
  norm_(false),
  drct_(false),
  ptrajformat_(false),
  vinfo1_(0),
  vinfo2_(0),
  tc_c_(0),
  tc_p_(0),
  tc_r3r3_(0),
  outfile_(0)
{}

const char* Analysis_Timecorr::Plegend_[] = {"<P0>", "<P1>", "<P2>"};

// Analysis_TimeCorr::CalculateAverages()
std::vector<double> Analysis_Timecorr::CalculateAverages(DataSet_Vector const& vIn, 
                                                         AvgResults& avgOut) 
{
  std::vector<double> R3i;
  R3i.reserve(vIn.Size());
  Vec3 avg(0.0, 0.0, 0.0);
  avgOut.rave_ = 0.0;
  avgOut.r3iave_ = 0.0;
  avgOut.r6iave_ = 0.0;
  // Loop over all vectors
  for (DataSet_Vector::const_iterator vec = vIn.begin();
                                      vec != vIn.end(); ++vec)
  {
    const Vec3& VXYZ = *vec;
    // Calc vector length
    double len = sqrt( VXYZ.Magnitude2() );
    // Update avgcrd, rave, r3iave, r6iave
    avg += VXYZ;
    avgOut.rave_ += len;
    double r3i = 1.0 / (len*len*len);
    avgOut.r3iave_ += r3i;
    avgOut.r6iave_ += r3i*r3i;
    R3i.push_back( r3i );
  }
  // Normalize averages
  double dnorm = 1.0 / (double)vIn.Size();
  avgOut.rave_ *= dnorm;
  avgOut.r3iave_ *= dnorm;
  avgOut.r6iave_ *= dnorm;
  avgOut.avgr_ = sqrt(avg.Magnitude2()) * dnorm;
  return R3i;
}

// Analysis_Timecorr::CalcCorr()
void Analysis_Timecorr::CalcCorr(int frame) {
  if (drct_) {
    // Calc correlation function using direct approach
    if (mode_ == AUTOCORR)
      corfdir_.AutoCorr(data1_);
    else
      corfdir_.CrossCorr(data1_, data2_);
  } else {
    // Pad with zero's at the end
    data1_.PadWithZero(frame);
    if (mode_ == AUTOCORR)
      pubfft_.AutoCorr(data1_);
    else {
      data2_.PadWithZero(frame);
      pubfft_.CrossCorr(data1_, data2_);
    }
  }
}

void Analysis_Timecorr::Help() const {
  mprintf("\tvec1 <vecname1> [vec2 <vecname2>] out <filename> [name <dsname>]\n"
          "\t[order <order>] [tstep <tstep>] [tcorr <tcorr>]\n"
          "\t[dplr] [norm] [drct] [dplrout <dplrfile>] [ptrajformat]\n"
          "  Calculate auto/cross-correlation functions for specified vectors.\n");
}

// Analysis_Timecorr::Setup()
Analysis::RetType Analysis_Timecorr::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  // Get Vectors
  std::string vec1name = analyzeArgs.GetStringKey("vec1");
  if (vec1name.empty()) {
    mprinterr("Error: no vec1 given, ignoring command\n");
    return Analysis::ERR;
  }
  vinfo1_ = (DataSet_Vector*)setup.DSL().FindSetOfGroup( vec1name, DataSet::VECTOR_1D );
  if (vinfo1_==0) {
    mprinterr("Error: vec1: no vector with name %s found.\n", 
              vec1name.c_str());
    return Analysis::ERR;
  }
  std::string vec2name = analyzeArgs.GetStringKey("vec2");
  if (!vec2name.empty()) {
    vinfo2_ = (DataSet_Vector*)setup.DSL().FindSetOfGroup( vec2name, DataSet::VECTOR_1D );
    if (vinfo2_==0) {
      mprinterr("Error: vec2: no vector with name %s found.\n", 
                vec2name.c_str());
      return Analysis::ERR;
    }
  } else
    vinfo2_ = 0;
  // Get output DataSet name
  std::string setname = analyzeArgs.GetStringKey("name");
  if (setname.empty())
    setname = setup.DSL().GenerateDefaultName("TC");
  // Determine auto or cross correlation 
  if (vinfo2_ == 0)
    mode_ = AUTOCORR;
  else
    mode_ = CROSSCORR;
  // Get dplr, norm, drct
  dplr_ = analyzeArgs.hasKey("dplr");
  norm_ = analyzeArgs.hasKey("norm");
  drct_ = analyzeArgs.hasKey("drct");
  std::string dplrname = analyzeArgs.GetStringKey("dplrout");
  // Get order for Legendre polynomial, tstep, and tcorr
  order_ = analyzeArgs.getKeyInt("order",2);
  if (order_ < 0 || order_ > 2) {
    mprintf("Warning: vector order out of bounds (should be 0, 1, or 2), resetting to 2.\n");
    order_ = 2;
  }
  tstep_ = analyzeArgs.getKeyDouble("tstep", 1.0);
  tcorr_ = analyzeArgs.getKeyDouble("tcorr", 10000.0);
  // File output. For ptrajformat, time correlation functions and dipolar
  // are output to file specified by 'out'. Otherwise time correlation
  // functions are written to file specified by 'out' using DataFile 
  // framework and dipolar output to 'dplrname'.
  ptrajformat_ = analyzeArgs.hasKey("ptrajformat");
  std::string filename = analyzeArgs.GetStringKey("out");
  if (ptrajformat_ && filename.empty()) {
    mprinterr("Error: No output file name given ('out <filename>'). Required for 'ptrajformat'.\n");
    return Analysis::ERR;
  }
  DataFile* dataout = 0;
  if (!ptrajformat_) {
    dataout = setup.DFL().AddDataFile( filename, analyzeArgs );
    if (dplr_) {
      if (!dplrname.empty() && dplrname == filename) {
        mprinterr("Error: 'dplrname' cannot be the same file as 'out' when 'ptrajformat' not specified.\n");
        return Analysis::ERR;
      }
      outfile_ = setup.DFL().AddCpptrajFile( dplrname, "Timecorr dipolar", DataFileList::TEXT, true );
      if (outfile_ == 0) return Analysis::ERR;
    }
  } else {
    outfile_ = setup.DFL().AddCpptrajFile( filename, "Timecorr output" );
    if (outfile_ == 0) return Analysis::ERR;
  }
  // Set up output DataSets
  tc_p_ = setup.DSL().AddSet( DataSet::DOUBLE, MetaData(setname, "P"));
  if (tc_p_ == 0) return Analysis::ERR;
  tc_p_->SetLegend( Plegend_[order_] );
  if (dataout != 0) dataout->AddDataSet( tc_p_ );
  if (dplr_) {
    tc_c_ = setup.DSL().AddSet( DataSet::DOUBLE, MetaData(setname, "C"));
    tc_r3r3_ = setup.DSL().AddSet( DataSet::DOUBLE, MetaData(setname, "R3R3"));
    if (tc_c_ == 0 || tc_r3r3_ == 0) return Analysis::ERR;
    tc_c_->SetLegend("<C>");
    tc_r3r3_->SetLegend( "<1/(r^3*r^3)>" );
    if (dataout != 0) {
      dataout->AddDataSet( tc_c_ );
      dataout->AddDataSet( tc_r3r3_ );
    }
  }

  // Print Status
  mprintf("    TIMECORR: Calculating %s", ModeString[mode_]);
  if (mode_ == AUTOCORR)
    mprintf(" of vector %s\n", vinfo1_->legend());
  else // CROSSCORR
    mprintf(" of vectors %s and %s\n", vinfo1_->legend(),
            vinfo2_->legend());
  mprintf("\tCorrelation time %f, time step %f, order %i\n", tcorr_, tstep_, order_);
  mprintf("\tCorr. func. are");
  if (dplr_)
    mprintf(" for dipolar interactions and");
  if (norm_)
    mprintf(" normalized.\n");
  else
    mprintf(" not normalized.\n");
  mprintf("\tCorr. func. are calculated using the");
  if (drct_)
    mprintf(" direct approach.\n");
  else
    mprintf(" FFT approach.\n");
  if (ptrajformat_)
    mprintf("\tResults are written to %s\n", outfile_->Filename().full());
  else {
    if (dataout != 0) mprintf("\tTime correlation functions written to %s\n", dataout->DataFilename().full());
    if (outfile_ != 0) mprintf("\tDipolar results written to %s\n", outfile_->Filename().full()); 
  }
  return Analysis::OK;
}

// Analysis_Timecorr::Analyze()
Analysis::RetType Analysis_Timecorr::Analyze() {
  // If 2 vectors, ensure they have the same # of frames
  if (vinfo2_!=0) {
    if (vinfo1_->Size() != vinfo2_->Size()) {
      mprinterr("Error: # Frames in vec %s (%zu) != # Frames in vec %s (%zu)\n",
                vinfo1_->legend(), vinfo1_->Size(), 
                vinfo2_->legend(), vinfo2_->Size());
      return Analysis::ERR;
    }
  }
  // Determine sizes
  int frame = vinfo1_->Size();
  int time = (int)(tcorr_ / tstep_) + 1;
  // nsteps
  int nsteps = 0;
  if (time > frame)
    nsteps = frame;
  else
    nsteps = time;
  // Allocate memory to hold complex numbers for direct or FFT
  if (drct_) {
    data1_.Allocate( frame );
    if (mode_ == CROSSCORR)
      data2_.Allocate( frame ); 
    corfdir_.CorrSetup( nsteps ); 
  } else {
    // Initialize FFT
    pubfft_.CorrSetup( frame );
    data1_ = pubfft_.Array();
    if (mode_ == CROSSCORR)
      data2_ = data1_;
  }
  // ----- Calculate spherical harmonics ---------
  // Real + Img. for each -order <= m <= order
  if (vinfo1_->CalcSphericalHarmonics(order_)) return Analysis::ERR;
  if (vinfo2_ != 0) {
    if (vinfo2_->CalcSphericalHarmonics(order_)) return Analysis::ERR;
  }
  // ----- Initialize PN output array memory -----
  DataSet_double& pncf_ = static_cast<DataSet_double&>( *tc_p_ );
  pncf_.Resize( nsteps );
  Dimension Xdim(0.0, tstep_, "Time");
  pncf_.SetDim(Dimension::X, Xdim);
  // ----- Calculate PN --------------------------
  for (int midx = -order_; midx <= order_; ++midx) {
    data1_.Assign( vinfo1_->SphericalHarmonics( midx ) );
    if (vinfo2_ != 0)
      data2_.Assign( vinfo2_->SphericalHarmonics( midx ) );
    CalcCorr( frame );
    for (int k = 0; k < nsteps; ++k)
      pncf_[k] += data1_[2 * k];
  }
  // ----- Dipolar Calc. -------------------------
  AvgResults Avg1, Avg2;
  if (dplr_) {
    DataSet_double& cf_ = static_cast<DataSet_double&>( *tc_c_ );
    cf_.Resize( nsteps );
    cf_.SetDim(Dimension::X, Xdim);
    DataSet_double& rcf_ = static_cast<DataSet_double&>( *tc_r3r3_ );
    rcf_.Resize( nsteps );
    rcf_.SetDim(Dimension::X, Xdim);
    // Calculate averages
    std::vector<double> R3i_1 = CalculateAverages(*vinfo1_, Avg1);
    std::vector<double> R3i_2;
    if (vinfo2_ != 0)
      R3i_2 = CalculateAverages(*vinfo2_, Avg2);
    // C
    for (int midx = -order_; midx <= order_; ++midx) {
      data1_.Assign( vinfo1_->SphericalHarmonics( midx ) );
      if (vinfo2_ != 0)
        data2_.Assign( vinfo2_->SphericalHarmonics( midx ) );
      for (int i = 0, i2 = 0; i < frame; ++i, i2 += 2) {
        double r3i = R3i_1[ i ]; 
        data1_[i2  ] *= r3i;
        data1_[i2+1] *= r3i;
        if ( vinfo2_ != 0 ) {
          r3i = R3i_2[ i ];
          data2_[i2  ] *= r3i;
          data2_[i2+1] *= r3i;
        }
      }
      CalcCorr( frame );
      for (int k = 0; k < nsteps; ++k) 
        cf_[k] += data1_[2 * k];
    }
    // 1 / R^6
    for (int i = 0, i2 = 0; i < frame; ++i, i2 += 2) {
      data1_[i2  ] = R3i_1[ i ];
      data1_[i2+1] = 0.0;
      if ( vinfo2_ != 0 ) {
        data2_[i2  ] = R3i_2[ i ];
        data2_[i2+1] = 0.0;
      }
    }
    CalcCorr( frame );
    for (int k = 0; k < nsteps; ++k)
      rcf_[k] = data1_[2 * k];
  }
  // ----- NORMALIZATION -------------------------
  // 4*PI / ((2*order)+1) due to spherical harmonics addition theorem
  double KN = DataSet_Vector::SphericalHarmonicsNorm( order_ );
  Normalize( tc_p_,    frame, KN );
  if (dplr_) {
    Normalize( tc_c_,    frame, KN );
    Normalize( tc_r3r3_, frame, 1.0 );
  }
  // ----- PRINT PTRAJ FORMAT --------------------
  if (outfile_ != 0) { 
    outfile_->Printf("%ss, normal type\n",ModeString[mode_]);
    if (dplr_) {
      outfile_->Printf("***** Vector length *****\n");
      outfile_->Printf("%10s %10s %10s %10s\n", "<r>", "<rrig>", "<1/r^3>", "<1/r^6>");
      outfile_->Printf("%10.4f %10.4f %10.4f %10.4f\n",
                     Avg1.rave_, Avg1.avgr_, Avg1.r3iave_, Avg1.r6iave_);
      if (mode_ == CROSSCORR)
        outfile_->Printf("%10.4f %10.4f %10.4f %10.4f\n",
                       Avg2.rave_, Avg2.avgr_, Avg2.r3iave_, Avg2.r6iave_);
    }
    if (ptrajformat_) {
      outfile_->Printf("\n***** Correlation functions *****\n");
      if (dplr_) {
        DataSet_double& cf_ = static_cast<DataSet_double&>( *tc_c_ );
        DataSet_double& rcf_ = static_cast<DataSet_double&>( *tc_r3r3_ );
        outfile_->Printf("%10s %10s %10s %10s\n", "Time", "<C>", Plegend_[order_], "<1/(r^3*r^3)>");
        for (int i = 0; i < nsteps; ++i)
          outfile_->Printf("%10.3f %10.4f %10.4f %10.4f\n", (double)i * tstep_,
                         cf_[i], pncf_[i], rcf_[i]);
      } else {
        outfile_->Printf("%10s %10s\n", "Time", Plegend_[order_]);
        for (int i = 0; i < nsteps; ++i)
          outfile_->Printf("%10.3f %10.4f\n", (double)i * tstep_, pncf_[i]);
      }
    }
  }
  return Analysis::OK;
}

// Analysis_Timecorr::Normalize()
void Analysis_Timecorr::Normalize(DataSet* ds, int frame, double Kin) {
  if (ds == 0) return;
  DataSet_double& data = static_cast<DataSet_double&>( *ds );
  double Kn;
  if (norm_)
    Kn = (double)frame / data[0];
  else
    Kn = Kin;
  int nsteps = (int)data.Size();
  for (int i = 0; i < nsteps; ++i) {
    data[i] *= (Kn / (double)(frame - i));
  }
}
