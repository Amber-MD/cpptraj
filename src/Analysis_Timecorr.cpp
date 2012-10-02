#include "Analysis_Timecorr.h"
#include "CpptrajStdio.h"
#include "Constants.h" // PI

/// Strings corresponding to modes, used in output.
const char Analysis_Timecorr::ModeString[2][6] = {
  "Auto", "Cross"
};

// CONSTRUCTOR
Analysis_Timecorr::Analysis_Timecorr() :
  tstep_(1),
  tcorr_(10000),
  order_(2),
  mode_(AUTO),
  dplr_(false),
  norm_(false),
  drct_(false),
  table_(0),
  data1_(0),
  data2_(0),
  cf_(0),
  p2cf_(0),
  rcf_(0),
  vinfo1_(0),
  vinfo2_(0)
{}

// DESTRUCTOR
Analysis_Timecorr::~Analysis_Timecorr() {
  if (table_!=0) delete[] table_;
  if (data1_!=0) delete[] data1_;
  if (data2_!=0) delete[] data2_;
  if (cf_!=0) delete[] cf_;
  if (p2cf_!=0) delete[] p2cf_;
  if (rcf_!=0) delete[] rcf_;
}

// Analysis_Timecorr::CalcCorr()
// TODO: Move to DataSet_Vector
void Analysis_Timecorr::CalcCorr(int ndata, int nsteps, int frame) {
  if (drct_) {
    // Calc correlation function using direct approach
    DataSet_Vector::corfdir(ndata, data1_, data2_, nsteps, table_);
  } else {
    // Pad with zero's at the end
    for (int k = 2 * frame; k < ndata; ++k) {
      data1_[k] = 0;
      if (data2_ != 0) data2_[k] = 0;
    }
    // Calc correlation function using FFT
    pubfft_.CorF_FFT(ndata, data1_, data2_);
  }
}


// Analysis_Timecorr::Setup()
/** analyze timecorr vec1 <vecname1> [vec2 <vecname2>]
  *                  [relax] [freq <hz>] [NHdist <distnh>] [order <order>]
  *                  tstep <tstep> tcorr <tcorr> out <filename>
  *                  [ corrired modes <modesname> [beg <ibeg> end <iend>] ]
  */
int Analysis_Timecorr::Setup(DataSetList* DSLin) {
  // Ensure first 2 args (should be 'analyze' 'timecorr') are marked.
  analyzeArgs_.MarkArg(0);
  analyzeArgs_.MarkArg(1);
  // Get Vectors
  std::string vec1name = analyzeArgs_.GetStringKey("vec1");
  if (vec1name.empty()) {
    mprinterr("Error: analyze timecorr: no vec1 given, ignoring command\n");
    return 1;
  }
  vinfo1_ = (DataSet_Vector*)DSLin->FindSetOfType( vec1name, DataSet::VECTOR );
  if (vinfo1_==0) {
    mprinterr("Error: analyze timecorr: vec1: no vector with name %s found.\n", 
              vec1name.c_str());
    return 1;
  }
  std::string vec2name = analyzeArgs_.GetStringKey("vec2");
  if (!vec2name.empty()) {
    vinfo2_ = (DataSet_Vector*)DSLin->FindSetOfType( vec2name, DataSet::VECTOR );
    if (vinfo2_==0) {
      mprinterr("Error: analyze timecorr: vec2: no vector with name %s found.\n", 
                vec2name.c_str());
      return 1;
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
  order_ = analyzeArgs_.getKeyInt("order",2);
  if (order_ < 0 || order_ > 2) {
    mprintf("Warning: vector order out of bounds (<0 or >2), resetting to 2.\n");
    order_ = 2;
  }

  // Get tstep, tcorr, filename
  tstep_ = analyzeArgs_.getKeyDouble("tstep", 1.0);
  tcorr_ = analyzeArgs_.getKeyDouble("tcorr", 10000.0);
  filename_ = analyzeArgs_.GetStringKey("out");
  if (filename_.empty()) {
    mprinterr("Error: analyze timecorr: No outfile given ('out <filename>').\n");
    return 1;
  }

  // Get dplr, norm, drct
  dplr_ = analyzeArgs_.hasKey("dplr");
  norm_ = analyzeArgs_.hasKey("norm");
  drct_ = analyzeArgs_.hasKey("drct");

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

  return 0;
}

// Analysis_Timecorr::Analyze()
int Analysis_Timecorr::Analyze() {
  // If 2 vectors, ensure they have the same # of frames
  if (vinfo2_!=0) {
    if (vinfo1_->Size() != vinfo2_->Size()) {
      mprinterr("Error: # Frames in vec %s (%i) != # Frames in vec %s (%i)\n",
                vinfo1_->Legend().c_str(), vinfo1_->Size(), 
                vinfo2_->Legend().c_str(), vinfo2_->Size());
      return 1;
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
  // ndata
  int ndata = 0;
  if (drct_) {
    ndata = 2 * frame;
    table_ = new double[ 2 * nsteps ];
  } else {
    // Initialize FFT
    pubfft_.SetupFFT( frame );
    ndata = pubfft_.size() * 2;
  }
  int mtot = 2 * order_ + 1;
  mprintf("CDBG: frame=%i time=%i nsteps=%i ndata=%i mtot=%i\n",frame,time,nsteps,ndata,mtot);

  // ---------------------------------------------------------------------------
  if (dplr_) {
    vinfo1_->CalculateAverages();
    if (vinfo2_ != 0)
      vinfo2_->CalculateAverages();
  }
  // Real + Img. for each -order <= m <= order, spherical Harmonics for each frame
  vinfo1_->CalcSphericalHarmonics(order_);
  if (vinfo2_ != 0)
    vinfo2_->CalcSphericalHarmonics(order_);
  // ---------------------------------------------------------------------------

  // Allocate common memory and initialize arrays
  data1_ = new double[ ndata ];
  data2_ = 0;
  if (mode_ == CROSS)
    data2_ = new double[ ndata ];

  // Initialize memory
  p2cf_ = new double[ nsteps ];
  for (int i = 0; i < nsteps; ++i)
    p2cf_[i] = 0;
  if (dplr_) {
    cf_ = new double[ nsteps ];
    rcf_ = new double[ nsteps ];
    for (int i = 0; i < nsteps; ++i) {
      cf_[i] = 0;
      rcf_[i] = 0;
    }
  }

  // P2
  for (int midx = -order_; midx <= order_; ++midx) {
    vinfo1_->FillData( data1_, midx );
    if (vinfo2_ != 0)
      vinfo2_->FillData( data2_, midx);
    CalcCorr( ndata, nsteps, frame );
    for (int k = 0; k < nsteps; ++k)
      p2cf_[k] += data1_[2 * k];
  }
  if (dplr_) {
    // C
    for (int midx = -order_; midx <= order_; ++midx) {
      vinfo1_->FillData( data1_, midx );
      if (vinfo2_ != 0)
        vinfo2_->FillData( data2_, midx);
      for (int i = 0; i < frame; ++i) {
        int i2 = i * 2;
        double r3i = vinfo1_->R3( i ); 
        data1_[i2  ] *= r3i;
        data1_[i2+1] *= r3i;
        if ( vinfo2_ != 0 ) {
          r3i = vinfo2_->R3( i );
          data2_[i2  ] *= r3i;
          data2_[i2+1] *= r3i;
        }
      }
      CalcCorr( ndata, nsteps, frame );
      for (int k = 0; k < nsteps; ++k) 
        cf_[k] += data1_[2 * k];
    }
    // 1 / R^6
    for (int i = 0; i < frame; ++i) {
      int i2 = i * 2;
      data1_[i2  ] = vinfo1_->R3( i );
      data1_[i2+1] = 0.0;
      if ( vinfo2_ != 0 ) {
        data2_[i2  ] = vinfo2_->R3( i );
        data2_[i2+1] = 0.0;
      }
    }
    CalcCorr( ndata, nsteps, frame );
    for (int k = 0; k < nsteps; ++k)
      rcf_[k] = data1_[2 * k];
  }
    
  // ----- PRINT NORMAL -----
  CpptrajFile outfile;
  if (outfile.OpenWrite(filename_)) return 1;
  outfile.Printf("%s-correlation functions, normal type\n",ModeString[mode_]);
  if (dplr_) {
    outfile.Printf("***** Vector length *****\n");
    outfile.Printf("%10s %10s %10s %10s\n", "<r>", "<rrig>", "<1/r^3>", "<1/r^6>");
    vinfo1_->PrintAvgcrd( outfile );
    if (vinfo2_ != 0)
      vinfo2_->PrintAvgcrd( outfile );
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

  return 0;
}

