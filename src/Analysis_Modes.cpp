#include "Analysis_Modes.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Analysis_Modes::Analysis_Modes() :
  type_(FLUCT),
  beg_(0),
  end_(0),
  bose_(false),
  factor_(0),
  modinfo_(0),
  source_(ModesInfo::MS_UNKNOWN),
  results_(0)
{}

/// Strings describing modes calculation types.
const char Analysis_Modes::analysisTypeString[3][22] = {
  "rms fluctuations",
  "displacements",
  "correlation functions"
};

// DESTRUCTOR
Analysis_Modes::~Analysis_Modes() {
  // Only delete ModesInfo if it was read in from a file.
  if (modinfo_ != 0 && source_==ModesInfo::MS_FILE)
    delete modinfo_;
  if (results_!=0)
    delete[] results_;
}

// Analysis_Modes::Setup()
int Analysis_Modes::Setup(DataSetList* DSLin) {
  // Analysis type
  if (analyzeArgs_.hasKey("fluct"))
    type_ = FLUCT;
  else if (analyzeArgs_.hasKey("displ"))
    type_ = DISPLACE;
  else if (analyzeArgs_.hasKey("corr"))
    type_ = CORR;
  else {
    mprinterr("Error: analyze modes: no analysis type specified.\n");
    return 1;
  }

  // Get beg, end, factor, bose
  beg_ = analyzeArgs_.getKeyInt("beg",7);
  end_ = analyzeArgs_.getKeyInt("end", 50);
  bose_ = analyzeArgs_.hasKey("bose");
  factor_ = analyzeArgs_.getKeyDouble("factor",1.0);

  // Get modes
  std::string modesfile = analyzeArgs_.GetStringKey("file");
  std::string modesname = analyzeArgs_.GetStringKey("stack");
  if (!modesfile.empty()) {
    modinfo_ = new ModesInfo();
    if (modinfo_->ReadEvecFile( modesfile, beg_, end_ )) return 1;
    if (modinfo_->Nvect() != (end_ - beg_ + 1)) {
      mprintf("Warning: analyze modes: # read evecs is %i, # requested is %i\n",
              modinfo_->Nvect(), (end_ - beg_ + 1));
    }
  } else if (!modesname.empty()) {
    DSLin->VectorBegin();
    while ( (modinfo_ = (ModesInfo*)DSLin->NextModes()) != 0 ) {
      if (modinfo_->Name() == modesname)
        break;
    }
    if (modinfo_ == 0) {
      mprinterr("Error: analyze modes: modes with name %s not found.\n", modesname.c_str());
      return 1;
    }
  } else {
    mprinterr("Error: analyze modes: Source 'stack' or 'file' not specified.\n");
    return 1;
  }
  source_ = modinfo_->Source();

  // Check modes type
  if (modinfo_->Mtype()!=ModesInfo::MT_COVAR && 
      modinfo_->Mtype()!=ModesInfo::MT_MWCOVAR)
  {
    mprinterr("Error: analyze modes: evecs must be of type COVAR or MWCOVAR.\n");
    return 1;
  }

  // Get output filename
  filename_ = analyzeArgs_.GetStringKey("out");

  // Get mask pair info for ANALYZEMODES_CORR option and build the atom pair stack
  if ( type_ == CORR ) {
    while (analyzeArgs_.hasKey("maskp")) {
      // Next two arguments should be one-atom masks
      char* a1mask = analyzeArgs_.getNextMask();
      char* a2mask = analyzeArgs_.getNextMask();
      if (a1mask==NULL || a2mask==NULL) {
        mprinterr("Error: analyze modes: For 'corr' two 1-atom masks are expected.\n");
        return 1;
      }
      // Check that each mask is just 1 atom
      AtomMask m1( a1mask );
      AtomMask m2( a2mask );
      analyzeParm_->SetupIntegerMask( m1 ); 
      analyzeParm_->SetupIntegerMask( m2 );
      if ( m1.Nselected()==1 && m2.Nselected()==1 )
        // Store atom pair
        // NOTE: use push_front to mimic stack behavior
        atompairStack_.push_front( std::pair<int,int>( m1[0], m2[0] ) );
      else {
        mprinterr("Error: analyze modes: For 'corr', masks should specify only one atom.\n");
        mprinterr("\tM1[%s]=%i atoms, M2[%s]=%i atoms.\n", m1.MaskString(), m1.Nselected(),
                  m2.MaskString(), m2.Nselected());
        return 1;
      }
    }
    if ( atompairStack_.empty() ) {
      mprinterr("Error: analyze modes: no atom pairs found (use 'maskp <mask1> <mask2>')\n");
      return 1;
    }
  }

  // Memory allocation for results vector is done in ModesInfo.

  // Status
  mprintf("    ANALYZE MODES: Calculating %s using modes %i to %i from %s\n",
          analysisTypeString, beg_, end_, modinfo_->c_str());
  mprintf("\tResults are written to");
  if (filename_.empty())
    mprintf(" STDOUT\n");
  else
    mprintf(" %s\n", filename_.c_str());
  if (bose_)
    mprintf("\tBose statistics used.\n");
  else
    mprintf("\tBoltzmann statistics used.\n");
  if (type_ == DISPLACE)
    mprintf("\tFactor for displacement: %lf\n", factor_);
  if (type_ == CORR) {
    mprintf("\tUsing the following atom pairs:");
    for (ModesInfo::modestack_it apair = atompairStack_.begin();
                                 apair != atompairStack_.end(); ++apair)
      mprintf(" (%i,%i)", (*apair).first+1, (*apair).second+1 );
    mprintf("\n");
  }

  return 0;
}

// Analysis_Modes::Analyze()
int Analysis_Modes::Analyze() {
  CpptrajFile outfile;
  // ----- FLUCT PRINT -----
  if (type_ == FLUCT) {
    // Calc rms atomic fluctuations
    results_ = modinfo_->CalcRMSfluct(beg_, end_, bose_);
    if (results_ == 0) return 1;
    outfile.OpenWrite( filename_ );
    outfile.Printf("Analysis of modes: RMS FLUCTUATIONS\n");
    outfile.Printf("%10s %10s %10s %10s %10s\n", "Atom no.", "rmsX", "rmsY", "rmsZ", "rms");
    int anum = 1;
    for (int i4 = 0; i4 < modinfo_->Navgelem()*4/3; i4+=4) 
      outfile.Printf("%10i %10.3f %10.3f %10.3f %10.3f\n", anum++, results_[i4], 
                     results_[i4+1], results_[i4+2], results_[i4+3]); 
  // ----- DISPLACE PRINT -----
  } else if (type_ == DISPLACE) {
    // Calc displacement of coordinates along normal mode directions
    results_ = modinfo_->CalcDisplacement(beg_, end_, bose_, factor_);
    if (results_ == 0) return 1;
    outfile.OpenWrite( filename_ );
    outfile.Printf("Analysis of modes: DISPLACEMENT\n");
    outfile.Printf("%10s %10s %10s %10s\n", "Atom no.", "displX", "displY", "displZ");
    int anum = 1;
    for (int i3 = 0; i3 < modinfo_->Navgelem(); i3 += 3)
      outfile.Printf("%10i %10.3f %10.3f %10.3f\n", anum++, results_[i3], 
                     results_[i3+1], results_[i3+2]);
  // ----- CORR PRINT -----
  } else if (type_ == CORR) {
    // Calc dipole-dipole correlation functions
    results_ = modinfo_->CalcDipoleCorr( beg_, end_, bose_, atompairStack_);
    if (results_==0) return 1;
    outfile.OpenWrite( filename_ );
    outfile.Printf("Analysis of modes: CORRELATION FUNCTIONS\n");
    outfile.Printf("%10s %10s %10s %10s %10s %10s\n", "Atom1", "Atom2", "Mode", 
                   "Freq", "1-S^2", "P2(cum)");
    int ncnt = 0;
    for (ModesInfo::modestack_it apair = atompairStack_.begin();
                                 apair != atompairStack_.end(); ++apair)
    {
      outfile.Printf("%10i %10i\n", (*apair).first+1, (*apair).second+1);
      double val = 1.0;
      for (int i = 0; i < modinfo_->Nvect(); ++i) {
        if (i+1 >= beg_ && i+1 <= end_) { // TODO: Pre-calc the shift
          double frq = modinfo_->Freq(i);
          if (frq >= 0.5) {
            val += results_[ncnt];
            outfile.Printf("%10s %10s %10i %10.5f %10.5f %10.5f\n",
                           "", "", i, frq, results_[ncnt], val);
            ++ncnt;
          }
        }
      }
    } // END loop over CORR atom pairs
  } else // SANITY CHECK
    return 1;
  outfile.CloseFile();

  return 0;
} 
