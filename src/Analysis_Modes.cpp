#include <cmath> // tanh, sqrt
#include <cstring> // memset
#include "Analysis_Modes.h"
#include "CpptrajStdio.h"
#include "Constants.h" // TWOPI

// CONSTRUCTOR
Analysis_Modes::Analysis_Modes() :
  type_(FLUCT),
  beg_(0),
  end_(0),
  bose_(false),
  factor_(0),
  modinfo_(0),
  results_(0)
{}

//#define TWOPI 6.2832
/// hc/2kT in cm, with T=300K; use for quantum Bose statistics)
const double Analysis_Modes::CONSQ = 2.39805E-3;
const double Analysis_Modes::TKBC2 = 0.46105E-34;
/// Avogadros number
const double Analysis_Modes::AVO   = 6.023E23;
const double Analysis_Modes::CNST  = TKBC2 * AVO;
// cm to angstroms
const double Analysis_Modes::CMTOA = 1.000E8;
const double Analysis_Modes::CONT  = CMTOA / TWOPI;
// NOTE: Original TWOPI was 6.2832, results in small roundoff diffs from ptraj

/// Strings describing modes calculation types.
const char Analysis_Modes::analysisTypeString[3][22] = {
  "rms fluctuations",
  "displacements",
  "correlation functions"
};

// DESTRUCTOR
Analysis_Modes::~Analysis_Modes() {
  if (results_!=0)
    delete[] results_;
}

void Analysis_Modes::CheckDeprecated(std::string& modesname, const char* key) {
  std::string arg = analyzeArgs_.GetStringKey( key );
  if (!arg.empty()) {
    mprintf("Warning: analyze modes: Argument %s is deprecated, use 'name <modes>' instead.\n",
            key);
    if (modesname.empty()) modesname = arg;
  }
}

// Analysis_Modes::Setup()
/** analyze modes fluct|displ|corr
  *                            name <modesname> 
  *                            [beg <beg>] [end <end>] 
  *                            [bose] [factor <factor>]
  *                            [out <outfile>] [maskp <mask1> <mask2> [...]]
  *  - fluct: rms fluctations from normal modes
  *  - displ: displacement of cartesian coordinates along normal mode directions
  *
  * results vector usage:
  *  - fluct:
  *      [rmsx(atom1), rmsy(atom1), rmsz(atom1), rms(atom1), ..., rmsx(atomN), ..., rms(atomN)]
  *  - displ:
  *      [displx(atom1), disply(atom1), displz(atom1), ..., displx(atomN), ..., displz(atomN)]
  *  - corr:
  *      [corr(pair1, vec1), ..., corr(pair1, vecN), ..., corr(pairM, vec1), ..., corr(pairM, vecN)
  */
int Analysis_Modes::Setup(DataSetList* DSLin) {
  // Ensure first 2 args (should be 'analyze' 'modes') are marked.
  analyzeArgs_.MarkArg(0);
  analyzeArgs_.MarkArg(1);
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
  bose_ = analyzeArgs_.hasKey("bose");
  factor_ = analyzeArgs_.getKeyDouble("factor",1.0);

  // Get modes name
  std::string modesfile = analyzeArgs_.GetStringKey("name");
  if (modesfile.empty()) {
    // Check for deprecated args
    CheckDeprecated(modesfile, "file");
    CheckDeprecated(modesfile, "stack");
    if (modesfile.empty()) {
      mprinterr("Error: analyze modes: No 'modes <name>' arg given.\n");
      return 1;
    }
  }
  beg_ = analyzeArgs_.getKeyInt("beg",7) - 1; // Args start at 1
  end_ = analyzeArgs_.getKeyInt("end", 50);
  // Check if modes name exists on the stack
  modinfo_ = (DataSet_Modes*)DSLin->FindSetOfType( modesfile, DataSet::MODES );
  if (modinfo_ == 0) {
    // If not on stack, check for file.
    if ( fileExists(modesfile.c_str()) ) {
      modinfo_ = (DataSet_Modes*)DSLin->AddSet( DataSet::MODES, modesfile, "Modes" );
      if (modinfo_->ReadEvecFile( modesfile, beg_+1, end_ )) return 1;
    }
  }
  if (modinfo_ == 0) {
    mprinterr("Error: analyze modes: Modes %s DataSet/file not found.\n");
    return 1;
  }

  // Check modes type
  if (modinfo_->Type() != DataSet_Matrix::COVAR && 
      modinfo_->Type() != DataSet_Matrix::MWCOVAR)
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
      ArgList::ConstArg a1mask = analyzeArgs_.getNextMask();
      ArgList::ConstArg a2mask = analyzeArgs_.getNextMask();
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

  // Status
  mprintf("    ANALYZE MODES: Calculating %s using modes %i to %i from %s\n",
          analysisTypeString[type_], beg_+1, end_, modinfo_->Legend().c_str());
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
    for (modestack_it apair = atompairStack_.begin();
                      apair != atompairStack_.end(); ++apair)
      mprintf(" (%i,%i)", (*apair).first+1, (*apair).second+1 );
    mprintf("\n");
  }

  return 0;
}

// Analysis_Modes::Analyze()
int Analysis_Modes::Analyze() {
  CpptrajFile outfile;
  // Check # of modes
  if (beg_ < 0 || beg_ >= modinfo_->Nmodes()) {
    mprinterr("Error: analyze modes: 'beg %i' is out of bounds.\n", beg_+1);
    return 1;
  }
  if (end_ >= modinfo_->Nmodes()) {
    mprintf("Warning: analyze modes: 'end %i' is > # of modes, setting to %i\n",
            end_, modinfo_->Nmodes());
    end_ = modinfo_->Nmodes();
  }
  if (end_ <= beg_) {
    mprinterr("Warning: analyze modes: beg must be <= end, (%i -- %i)\n", beg_+1, end_);
    return 1;
  }

  // ----- FLUCT PRINT -----
  if (type_ == FLUCT) {
    // Calc rms atomic fluctuations
    int natoms = modinfo_->NavgCrd() / 3; // Possible because COVAR/MWCOVAR
    results_ = new double[ natoms * 4 ];
    memset(results_, 0, natoms * 4 * sizeof(double));
    double* Ri = results_;
    for (int vi = 0; vi < natoms; ++vi) {
      double sumx = 0.0;
      double sumy = 0.0;
      double sumz = 0.0;
      // Loop over all eigenvector elements for this atom
      const double* Vec = modinfo_->Eigenvector(beg_) + vi*3;
      for (int mode = beg_; mode < end_; ++mode) {
        double frq = modinfo_->Eigenvalue(mode);
        if (frq >= 0.5) {
          // Don't use eigenvectors associated with zero or negative eigenvalues
          double distx = Vec[0] * Vec[0];
          double disty = Vec[1] * Vec[1];
          double distz = Vec[2] * Vec[2];
          double fi = 1.0 / (frq*frq);
          if (bose_) {
            double argq = CONSQ * frq;
            fi *= (argq / tanh(argq));
          }
          sumx += distx * fi;
          sumy += disty * fi;
          sumz += distz * fi;
        }
        Vec += modinfo_->VectorSize();
      }
      sumx *= CNST;
      sumy *= CNST;
      sumz *= CNST;
      Ri[0] = sqrt(sumx) * CONT;
      Ri[1] = sqrt(sumy) * CONT;
      Ri[2] = sqrt(sumz) * CONT;
      Ri[3] = sqrt(sumx + sumy + sumz) * CONT;
      Ri += 4;
    }
    // Output
    outfile.OpenWrite( filename_ );
    outfile.Printf("Analysis of modes: RMS FLUCTUATIONS\n");
    outfile.Printf("%10s %10s %10s %10s %10s\n", "Atom no.", "rmsX", "rmsY", "rmsZ", "rms");
    int anum = 1;
    for (int i4 = 0; i4 < modinfo_->NavgCrd()*4/3; i4+=4) 
      outfile.Printf("%10i %10.3f %10.3f %10.3f %10.3f\n", anum++, results_[i4], 
                     results_[i4+1], results_[i4+2], results_[i4+3]); 
  // ----- DISPLACE PRINT -----
  } else if (type_ == DISPLACE) {
    // Calc displacement of coordinates along normal mode directions
    results_ = new double[ modinfo_->NavgCrd() ];
    memset(results_, 0, modinfo_->NavgCrd()*sizeof(double));
    double sqrtcnst = sqrt(CNST) * CONT * factor_;
    // Loop over all modes
    for (int mode = beg_; mode < end_; ++mode) {
      double frq = modinfo_->Eigenvalue(mode);
      if (frq >= 0.5) {
        // Don't use eigenvectors associated with zero or negative eigenvalues
        double fi = 1.0 / frq;
        if (bose_) {
          double argq = CONSQ * frq;
          fi *= (fi * argq / tanh(argq));
          fi = sqrt(fi);
        }
        fi *= sqrtcnst; // * CONT * factor
        // Loop over all vector elements
        const double* Vec = modinfo_->Eigenvector(mode);
        for (int vi = 0; vi < modinfo_->NavgCrd(); vi += 3) {
          results_[vi  ] += Vec[vi  ] * fi;
          results_[vi+1] += Vec[vi+1] * fi;
          results_[vi+2] += Vec[vi+2] * fi;
        }
      }
    }
    // Output
    outfile.OpenWrite( filename_ );
    outfile.Printf("Analysis of modes: DISPLACEMENT\n");
    outfile.Printf("%10s %10s %10s %10s\n", "Atom no.", "displX", "displY", "displZ");
    int anum = 1;
    for (int i3 = 0; i3 < modinfo_->NavgCrd(); i3 += 3)
      outfile.Printf("%10i %10.3f %10.3f %10.3f\n", anum++, results_[i3], 
                     results_[i3+1], results_[i3+2]);
  // ----- CORR PRINT -----
  } else if (type_ == CORR) {
    // Calc dipole-dipole correlation functions
    //results_ = modinfo_->CalcDipoleCorr( beg_, end_, bose_, atompairStack_); // TODO
    if (results_==0) return 1;
    outfile.OpenWrite( filename_ );
    outfile.Printf("Analysis of modes: CORRELATION FUNCTIONS\n");
    outfile.Printf("%10s %10s %10s %10s %10s %10s\n", "Atom1", "Atom2", "Mode", 
                   "Freq", "1-S^2", "P2(cum)");
    int ncnt = 0;
    for (modestack_it apair = atompairStack_.begin();
                      apair != atompairStack_.end(); ++apair)
    {
      outfile.Printf("%10i %10i\n", (*apair).first+1, (*apair).second+1);
      double val = 1.0;
      for (int i = 0; i < modinfo_->Nmodes(); ++i) {
        if (i+1 >= beg_ && i+1 <= end_) { // TODO: Pre-calc the shift
          double frq = modinfo_->Eigenvalue(i);
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
