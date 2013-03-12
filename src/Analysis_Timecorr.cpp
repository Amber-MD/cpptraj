#include <cmath> // sqrt
#include "Analysis_Timecorr.h"
#include "CpptrajStdio.h"
#include "Constants.h" // PI

/// Strings corresponding to modes, used in output.
const char* Analysis_Timecorr::ModeString[] = { "Auto", "Cross" };

// CONSTRUCTOR
Analysis_Timecorr::Analysis_Timecorr() :
  tstep_(1),
  tcorr_(10000),
  order_(2),
  mode_(AUTO),
  dplr_(false),
  norm_(false),
  drct_(false),
  vinfo1_(0),
  vinfo2_(0)
{}

void Analysis_Timecorr::Help() {
  mprintf("\tvec1 <vecname1> [vec2 <vecname2>] out <filename>\n");
  mprintf("\t[order <order>] tstep <tstep> tcorr <tcorr>\n");
  mprintf("\t[dplr] [norm] [drct]\n");
  mprintf("\tCalculate auto or cross-correlation function for specified vectors.\n");
}

// Analysis_TimeCorr::CalculateAverages()
std::vector<double> Analysis_Timecorr::CalculateAverages(DataSet_Vector& vIn, 
                                                         AvgResults& avgOut) 
{
  std::vector<double> R3i;
  R3i.reserve(vIn.Size());
  double avgx = 0.0;
  double avgy = 0.0;
  double avgz = 0.0;
  avgOut.rave_ = 0.0;
  avgOut.r3iave_ = 0.0;
  avgOut.r6iave_ = 0.0;
  // Loop over all vectors
  for (DataSet_Vector::iterator vec = vIn.begin();
                                vec != vIn.end(); ++vec)
  {
    const double* VXYZ = *vec;
    // Calc vector length
    double len = sqrt(VXYZ[0]*VXYZ[0] + VXYZ[1]*VXYZ[1] + VXYZ[2]*VXYZ[2]);
    // Update avgcrd, rave, r3iave, r6iave
    avgx += VXYZ[0];
    avgy += VXYZ[1];
    avgz += VXYZ[2];
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
  avgOut.avgr_ = sqrt(avgx*avgx + avgy*avgy + avgz*avgz) * dnorm;
  return R3i;
}

// Analysis_Timecorr::CalcCorr()
void Analysis_Timecorr::CalcCorr(int frame) {
  if (drct_) {
    // Calc correlation function using direct approach
    if (mode_ == AUTO)
      corfdir_.AutoCorr(data1_);
    else
      corfdir_.CrossCorr(data1_, data2_);
  } else {
    // Pad with zero's at the end
    data1_.PadWithZero(frame);
    if (mode_ == AUTO)
      pubfft_.AutoCorr(data1_);
    else {
      data2_.PadWithZero(frame);
      pubfft_.CrossCorr(data1_, data2_);
    }
  }
}

// Analysis_Timecorr::Setup()
Analysis::RetType Analysis_Timecorr::Setup(ArgList& analyzeArgs, DataSetList* DSLin,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  // Get Vectors
  std::string vec1name = analyzeArgs.GetStringKey("vec1");
  if (vec1name.empty()) {
    mprinterr("Error: analyze timecorr: no vec1 given, ignoring command\n");
    return Analysis::ERR;
  }
  vinfo1_ = (DataSet_Vector*)DSLin->FindSetOfType( vec1name, DataSet::VECTOR );
  if (vinfo1_==0) {
    mprinterr("Error: analyze timecorr: vec1: no vector with name %s found.\n", 
              vec1name.c_str());
    return Analysis::ERR;
  }
  std::string vec2name = analyzeArgs.GetStringKey("vec2");
  if (!vec2name.empty()) {
    vinfo2_ = (DataSet_Vector*)DSLin->FindSetOfType( vec2name, DataSet::VECTOR );
    if (vinfo2_==0) {
      mprinterr("Error: analyze timecorr: vec2: no vector with name %s found.\n", 
                vec2name.c_str());
      return Analysis::ERR;
    }
  }
  // Determine auto or cross correlation 
  if (vinfo2_ == 0)
    mode_ = AUTO;
  //else if (vinfo1_ == vinfo2_) {
  //  mode_ = AUTO; 
  //  vinfo2_ = 0;}
  else
    mode_ = CROSS;
  // Get order for Legendre polynomial
  order_ = analyzeArgs.getKeyInt("order",2);
  if (order_ < 0 || order_ > 2) {
    mprintf("Warning: vector order out of bounds (<0 or >2), resetting to 2.\n");
    order_ = 2;
  }

  // Get tstep, tcorr, filename
  tstep_ = analyzeArgs.getKeyDouble("tstep", 1.0);
  tcorr_ = analyzeArgs.getKeyDouble("tcorr", 10000.0);
  filename_ = analyzeArgs.GetStringKey("out");
  if (filename_.empty()) {
    mprinterr("Error: analyze timecorr: No outfile given ('out <filename>').\n");
    return Analysis::ERR;
  }

  // Get dplr, norm, drct
  dplr_ = analyzeArgs.hasKey("dplr");
  norm_ = analyzeArgs.hasKey("norm");
  drct_ = analyzeArgs.hasKey("drct");

  // Print Status
  mprintf("    ANALYZE TIMECORR: Calculating");
  if (mode_ == AUTO)
    mprintf(" auto-correlation function of vector %s", vinfo1_->Legend().c_str());
  else if (mode_ == CROSS)
    mprintf(" cross-correlation function of vectors %s and %s", vinfo1_->Legend().c_str(),
            vinfo2_->Legend().c_str());
  mprintf("\n\t\tCorrelation time %lf, time step %lf\n", tcorr_, tstep_);
  mprintf("\t\tCorr. func. are");
  if (dplr_)
    mprintf(" for dipolar interactions and");
  if (norm_)
    mprintf(" normalized.\n");
  else
    mprintf(" not normalized.\n");
  mprintf("\t\tCorr. func. are calculated using the");
  if (drct_)
    mprintf(" direct approach.\n");
  else
    mprintf(" FFT approach.\n");
  mprintf("\t\tResults are written to %s\n", filename_.c_str());

  return Analysis::OK;
}

// Analysis_Timecorr::Analyze()
Analysis::RetType Analysis_Timecorr::Analyze() {
  // If 2 vectors, ensure they have the same # of frames
  if (vinfo2_!=0) {
    if (vinfo1_->Size() != vinfo2_->Size()) {
      mprinterr("Error: # Frames in vec %s (%i) != # Frames in vec %s (%i)\n",
                vinfo1_->Legend().c_str(), vinfo1_->Size(), 
                vinfo2_->Legend().c_str(), vinfo2_->Size());
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
    if (mode_ == CROSS)
      data2_.Allocate( frame ); 
    corfdir_.Allocate( nsteps ); 
  } else {
    // Initialize FFT
    pubfft_.Allocate( frame );
    data1_ = pubfft_.Array();
    if (mode_ == CROSS)
      data2_ = data1_;
  }

  // ---------------------------------------------------------------------------
  // Real + Img. for each -order <= m <= order, spherical Harmonics for each frame
  vinfo1_->CalcSphericalHarmonics(order_);
  if (vinfo2_ != 0)
    vinfo2_->CalcSphericalHarmonics(order_);
  // ---------------------------------------------------------------------------

  // Initialize output array memory
  std::vector<double> p2cf_(nsteps, 0.0);
  std::vector<double> cf_;
  std::vector<double> rcf_;
  if (dplr_) {
    cf_.assign(nsteps, 0.0);
    rcf_.assign(nsteps, 0.0);
  }

  // P2
  for (int midx = -order_; midx <= order_; ++midx) {
    vinfo1_->FillData( data1_, midx );
    if (vinfo2_ != 0)
      vinfo2_->FillData( data2_, midx);
    CalcCorr( frame );
    for (int k = 0; k < nsteps; ++k)
      p2cf_[k] += data1_[2 * k];
  }
  // Only needed if dplr
  AvgResults Avg1, Avg2;
  if (dplr_) {
    // Calculate averages
    std::vector<double> R3i_1 = CalculateAverages(*vinfo1_, Avg1);
    std::vector<double> R3i_2;
    if (vinfo2_ != 0)
      R3i_2 = CalculateAverages(*vinfo2_, Avg2);
    // C
    for (int midx = -order_; midx <= order_; ++midx) {
      vinfo1_->FillData( data1_, midx );
      if (vinfo2_ != 0)
        vinfo2_->FillData( data2_, midx);
      int i2 = 0;
      for (int i = 0; i < frame; ++i) {
        double r3i = R3i_1[ i ]; 
        data1_[i2  ] *= r3i;
        data1_[i2+1] *= r3i;
        if ( vinfo2_ != 0 ) {
          r3i = R3i_2[ i ];
          data2_[i2  ] *= r3i;
          data2_[i2+1] *= r3i;
        }
        i2 += 2;
      }
      CalcCorr( frame );
      for (int k = 0; k < nsteps; ++k) 
        cf_[k] += data1_[2 * k];
    }
    // 1 / R^6
    for (int i = 0; i < frame; ++i) {
      int i2 = i * 2;
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
    
  // ----- PRINT NORMAL -----
  CpptrajFile outfile;
  if (outfile.OpenWrite(filename_)) return Analysis::ERR;
  outfile.Printf("%s-correlation functions, normal type\n",ModeString[mode_]);
  if (dplr_) {
    outfile.Printf("***** Vector length *****\n");
    outfile.Printf("%10s %10s %10s %10s\n", "<r>", "<rrig>", "<1/r^3>", "<1/r^6>");
    outfile.Printf("%10.4f %10.4f %10.4f %10.4f\n",
                   Avg1.rave_, Avg1.avgr_, Avg1.r3iave_, Avg1.r6iave_);
    if (mode_ == CROSS)
      outfile.Printf("%10.4f %10.4f %10.4f %10.4f\n",
                     Avg2.rave_, Avg2.avgr_, Avg2.r3iave_, Avg2.r6iave_);
  }
  outfile.Printf("\n***** Correlation functions *****\n");
  if (dplr_) {
    outfile.Printf("%10s %10s %10s %10s\n", "Time", "<C>", "<P2>", "<1/(r^3*r^3)>");
    if (norm_) {
      for (int i = 0; i < nsteps; ++i)
        outfile.Printf("%10.3f %10.4f %10.4f %10.4f\n",
                       (double)i * tstep_,
                       cf_[i]   * frame / (cf_[0]   * (frame - i)),
                       p2cf_[i] * frame / (p2cf_[0] * (frame - i)),
                       rcf_[i]  * frame / (rcf_[0]  * (frame - i)));  
    } else {
      // 4/5*PI due to spherical harmonics addition theorem
      for (int i = 0; i < nsteps; ++i)
        outfile.Printf("%10.3f %10.4f %10.4f %10.4f\n",
                       (double)i * tstep_,
                       FOURFIFTHSPI * cf_[i]   / (frame - i),
                       FOURFIFTHSPI * p2cf_[i] / (frame - i),
                       rcf_[i]  / (frame - i));
    }
  } else {
    outfile.Printf("%10s %10s\n", "Time", "<P2>");
    if (norm_) {
      for (int i = 0; i < nsteps; ++i)
        outfile.Printf("%10.3f %10.4f\n",
                       (double)i * tstep_,
                       p2cf_[i] * frame / (p2cf_[0] * (frame - i)));
    } else {
      // 4/5*PI due to spherical harmonics addition theorem
      for (int i = 0; i < nsteps; ++i)
        outfile.Printf("%10.3f %10.4f\n",
                       (double)i * tstep_,
                       FOURFIFTHSPI * p2cf_[i] / (frame - i));
    }
  }
  outfile.CloseFile();

  return Analysis::OK;
}
