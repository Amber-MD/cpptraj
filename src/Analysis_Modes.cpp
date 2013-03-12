#include <cmath> // tanh, sqrt
#include <cstring> // memset
#include "Analysis_Modes.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // fileExists
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

void Analysis_Modes::Help() {
  mprintf("\t{fluct|displ|corr} name <modesname>\n"); 
  mprintf("\t[beg <beg>] [end <end>] [bose] [factor <factor>]\n");
  mprintf("\t[out <outfile>] [maskp <mask1> <mask2> [...]]\n");
  mprintf("\tPerform one of the following analysis on calculated Eigenmodes.\n");
  mprintf("\tfluct: rms fluctations from normal modes\n");
  mprintf("\tdispl: displacement of cartesian coordinates along normal mode directions\n");
  mprintf(" Results vector usage:\n");
  mprintf("   fluct:\n");
  mprintf("\t[rmsx(at1), rmsy(at1), rmsz(at1), rms(at1), ..., rmsx(atN), ..., rms(atN)]\n");
  mprintf("   displ:\n");
  mprintf("\t[displx(at1), disply(at1), displz(at1), ..., displx(atN), ..., displz(atN)]\n");
  mprintf("   corr:\n");
  mprintf("\t[corr(pair1, vec1), ..., corr(pair1, vecN), ..., corr(pairM, vec1), ..., corr(pairM, vecN)\n");
}

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
const char* Analysis_Modes::analysisTypeString[] = {
  "rms fluctuations",
  "displacements",
  "correlation functions"
};

// DESTRUCTOR
Analysis_Modes::~Analysis_Modes() {
  if (results_!=0)
    delete[] results_;
}

// Analysis_Modes::CheckDeprecated()
void Analysis_Modes::CheckDeprecated(ArgList& analyzeArgs, std::string& modesname, 
                                     const char* key) {
  std::string arg = analyzeArgs.GetStringKey( key );
  if (!arg.empty()) {
    mprintf("Warning: analyze modes: Argument '%s' is deprecated, use 'name <modes>' instead.\n",
            key);
    if (modesname.empty()) modesname = arg;
  }
}

// Analysis_Modes::Setup()
Analysis::RetType Analysis_Modes::Setup(ArgList& analyzeArgs, DataSetList* DSLin,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  // Analysis type
  if (analyzeArgs.hasKey("fluct"))
    type_ = FLUCT;
  else if (analyzeArgs.hasKey("displ"))
    type_ = DISPLACE;
  else if (analyzeArgs.hasKey("corr"))
    type_ = CORR;
  else {
    mprinterr("Error: analyze modes: no analysis type specified.\n");
    return Analysis::ERR;
  }

  // Get beg, end, factor, bose
  bose_ = analyzeArgs.hasKey("bose");
  factor_ = analyzeArgs.getKeyDouble("factor",1.0);

  // Get modes name
  std::string modesfile = analyzeArgs.GetStringKey("name");
  if (modesfile.empty()) {
    // Check for deprecated args
    CheckDeprecated(analyzeArgs, modesfile, "file");
    CheckDeprecated(analyzeArgs, modesfile, "stack");
    if (modesfile.empty()) {
      mprinterr("Error: analyze modes: No 'modes <name>' arg given.\n");
      return Analysis::ERR;
    }
  }
  beg_ = analyzeArgs.getKeyInt("beg",7) - 1; // Args start at 1
  end_ = analyzeArgs.getKeyInt("end", 50);
  // Check if modes name exists on the stack
  modinfo_ = (DataSet_Modes*)DSLin->FindSetOfType( modesfile, DataSet::MODES );
  if (modinfo_ == 0) {
    // If not on stack, check for file.
    if ( fileExists(modesfile.c_str()) ) {
      modinfo_ = (DataSet_Modes*)DSLin->AddSet( DataSet::MODES, modesfile, "Modes" );
      if (modinfo_->ReadEvecFile( modesfile, beg_+1, end_ )) return Analysis::ERR;
    }
  }
  if (modinfo_ == 0) {
    mprinterr("Error: analyze modes: Modes %s DataSet/file not found.\n");
    return Analysis::ERR;
  }

  // Check modes type
  if (modinfo_->Type() != DataSet_Matrix::COVAR && 
      modinfo_->Type() != DataSet_Matrix::MWCOVAR)
  {
    mprinterr("Error: analyze modes: evecs must be of type COVAR or MWCOVAR.\n");
    return Analysis::ERR;
  }

  // Get output filename
  filename_ = analyzeArgs.GetStringKey("out");

  // Get mask pair info for ANALYZEMODES_CORR option and build the atom pair stack
  if ( type_ == CORR ) {
    Topology* analyzeParm = PFLin->GetParm( analyzeArgs );
    if (analyzeParm == 0) {
      mprinterr("Error: analyze modes: corr requires topology (parm <file>, parmindex <#>).\n");
      return Analysis::ERR;
    }
    while (analyzeArgs.hasKey("maskp")) {
      // Next two arguments should be one-atom masks
      std::string a1mask = analyzeArgs.GetMaskNext();
      std::string a2mask = analyzeArgs.GetMaskNext();
      if (a1mask.empty() || a2mask.empty()) {
        mprinterr("Error: analyze modes: For 'corr' two 1-atom masks are expected.\n");
        return Analysis::ERR;
      }
      // Check that each mask is just 1 atom
      AtomMask m1( a1mask );
      AtomMask m2( a2mask );
      analyzeParm->SetupIntegerMask( m1 ); 
      analyzeParm->SetupIntegerMask( m2 );
      if ( m1.Nselected()==1 && m2.Nselected()==1 )
        // Store atom pair
        atompairStack_.push_back( std::pair<int,int>( m1[0], m2[0] ) );
      else {
        mprinterr("Error: analyze modes: For 'corr', masks should specify only one atom.\n");
        mprinterr("\tM1[%s]=%i atoms, M2[%s]=%i atoms.\n", m1.MaskString(), m1.Nselected(),
                  m2.MaskString(), m2.Nselected());
        return Analysis::ERR;
      }
    }
    if ( atompairStack_.empty() ) {
      mprinterr("Error: analyze modes: no atom pairs found (use 'maskp <mask1> <mask2>')\n");
      return Analysis::ERR;
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

  return Analysis::OK;
}

// Analysis_Modes::Analyze()
Analysis::RetType Analysis_Modes::Analyze() {
  CpptrajFile outfile;
  // Check # of modes
  if (beg_ < 0 || beg_ >= modinfo_->Nmodes()) {
    mprinterr("Error: analyze modes: 'beg %i' is out of bounds.\n", beg_+1);
    return Analysis::ERR;
  }
  if (end_ > modinfo_->Nmodes()) {
    mprintf("Warning: analyze modes: 'end %i' is > # of modes, setting to %i\n",
            end_, modinfo_->Nmodes());
    end_ = modinfo_->Nmodes();
  }
  if (end_ <= beg_) {
    mprinterr("Warning: analyze modes: beg must be <= end, (%i -- %i)\n", beg_+1, end_);
    return Analysis::ERR;
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
    CalcDipoleCorr();
    if (results_==0) return Analysis::ERR;
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
      for (int mode = beg_; mode < end_; ++mode) {
        double frq = modinfo_->Eigenvalue(mode);
        if (frq >= 0.5) {
          val += results_[ncnt];
          outfile.Printf("%10s %10s %10i %10.5f %10.5f %10.5f\n",
                         "", "", mode, frq, results_[ncnt], val);
          ++ncnt;
        }
      }
    } // END loop over CORR atom pairs
  } else // SANITY CHECK
    return Analysis::ERR;
  outfile.CloseFile();

  return Analysis::OK;
}

// Analysis_Modes::CalcDipoleCorr()
void Analysis_Modes::CalcDipoleCorr() {
  double qcorr = 1.0; // For when bose is false
  int rsize = atompairStack_.size() * (end_ - beg_ + 1);
  results_ = new double[ rsize ];
  memset(results_, 0, rsize*sizeof(double));
  // Loop over atom pairs 
  double* Res = results_;
  const double* Avg = modinfo_->AvgCrd();
  for (modestack_it apair = atompairStack_.begin(); apair != atompairStack_.end(); ++apair)
  {
    int idx1 = (*apair).first  * 3;
    int idx2 = (*apair).second * 3;
    // Check if coordinates are out of bounds
    if (idx1 >= modinfo_->NavgCrd() || idx2 >= modinfo_->NavgCrd()) {
      mprintf("Warning: Atom pair %i -- %i is out of bounds (# avg atoms in modes = %i)\n",
              (*apair).first + 1, (*apair).second + 1, modinfo_->NavgCrd() / 3);
      continue;
    }
    // Calc unit vector along at2->at1 bond
    Vec3 vec = Vec3(Avg + idx1) - Vec3(Avg + idx2);
    double dnorm = sqrt( vec.Magnitude2() );
    vec /= dnorm;
    // Precalc certain values
    dnorm = 3.0 / (dnorm * dnorm);
    double vx2 = vec[0] * vec[0];
    double vxy = vec[0] * vec[1];
    double vxz = vec[0] * vec[2];
    double vy2 = vec[1] * vec[1];
    double vyz = vec[1] * vec[2];
    double vz2 = vec[2] * vec[2]; 
    // Loop over desired modes
    for (int mode = beg_; mode < end_; ++mode) {
      double eval = modinfo_->Eigenvalue(mode);
      if (eval >= 0.5) {
        // Don't use eigenvectors associated with zero or negative eigenvalues
        // NOTE: Where is 11791.79 from? Should it be a const?
        double frq = eval * eval / 11791.79;
        if (bose_) {
          double argq = CONSQ * eval;
          qcorr = argq / tanh(argq);
        }
        /* Calc the correlation matrix for delta
         *    as in eq. 7.16 of lamm and szabo, J Chem Phys 1986, 85, 7334.
         *  Note that the rhs of this eq. should be multiplied by kT
         */
        double qcorrf = (qcorr / frq) * 0.6;
        const double* Evec = modinfo_->Eigenvector(mode);
        double dx = (Evec[idx1  ] - Evec[idx2  ]);
        double dy = (Evec[idx1+1] - Evec[idx2+1]);
        double dz = (Evec[idx1+2] - Evec[idx2+2]);
        double delx2 = qcorrf * dx * dx;
        double delxy = qcorrf * dx * dy;
        double delxz = qcorrf * dx * dz;
        double dely2 = qcorrf * dy * dy;
        double delyz = qcorrf * dy * dz;
        double delz2 = qcorrf * dz * dz;
        // Correlation in length, eq. 10.2 of lamm and szabo
        // NOTE: Commented out in PTRAJ
        /*****
        rtr0 = 0.0;
        for(j = 0; j < 3; j++)
          for(k = 0; k < 3; k++)
            rtr0 += e[j] * e[k] * del[j][k];
        *****/
        // Librational correlation function, using eq. 7.12 of lamm and szabo
        // (w/o beta on the lhs).
        *(Res++) = (-delx2 + (vx2 * delx2) + (vxy * delxy) + (vxz * delxz)
                    -dely2 + (vxy * delxy) + (vy2 * dely2) + (vyz * delyz)
                    -delz2 + (vxz * delxz) + (vyz * delyz) + (vz2 * delz2))
                   * dnorm; // Here dnorm is actually (3 / dnorm^2)
      } // END if positive definite eigenvalue
    } // END loop over modes
  } // END loop over atom pairs
}
