#include <cmath> // tanh, sqrt
#include "Analysis_Modes.h"
#include "CpptrajStdio.h"
#include "Constants.h" // TWOPI

// CONSTRUCTOR
Analysis_Modes::Analysis_Modes() :
  debug_(0),
  type_(FLUCT),
  beg_(0),
  end_(0),
  bose_(false),
  calcAll_(false),
  factor_(0),
  modinfo_(0),
  modinfo2_(0),
  outfile_(0),
  tOutParm_(0),
  tMode_(0),
  pcmin_(0.0),
  pcmax_(0.0)
{}

void Analysis_Modes::Help() const {
  mprintf("\t{fluct|displ|corr|eigenval|trajout|rmsip} name <modesname> [name2 <modesname>]\n"
          "\t[beg <beg>] [end <end>] [bose] [factor <factor>] [calcall]\n"
          "\t[out <outfile>] [setname <name>]\n"
          "    Options for 'trajout': (Generate pseudo-trajectory)\n"
          "\t[trajout <name> %s\n", DataSetList::TopArgs);
  mprintf("\t  [trajoutfmt <format>] [trajoutmask <mask>]\n"
          "\t  [pcmin <pcmin>] [pcmax <pcmax>] [tmode <mode>]]\n"
          "    Options for 'corr': (Calculate dipole correlation)\n"
          "\t{ maskp <mask1> <mask2> [...] | mask1 <mask> mask2 <mask> }\n"
          "\t%s\n", DataSetList::TopArgs);
  mprintf("  Perform one of the following analysis on calculated Eigenmodes.\n"
          "    fluct    : RMS fluctations from normal modes.\n"
          "    displ    : Displacement of cartesian coordinates along normal mode directions.\n"
          "    corr     : Calculate dipole-dipole correlation functions.\n"
          "    eigenval : Calculate eigenvalue fraction for all modes.\n"
          "    trajout  : Calculate pseudo-trajectory along given mode 'tmode'.\n"
          "    rmsip    : Root mean square inner product between 'name' and 'name2'.\n");
}

/// hc/2kT in cm, with T=300K; use for quantum Bose statistics)
const double Analysis_Modes::CONSQ = 2.39805E-3;
/// kT/c^2 in cgs units (grams), with T=300K
const double Analysis_Modes::TKBC2 = 0.46105E-34;
/// Avogadros number
const double Analysis_Modes::AVO   = 6.023E23;
const double Analysis_Modes::CNST  = TKBC2 * AVO;
// cm to angstroms
const double Analysis_Modes::CMTOA = 1.000E8;
const double Analysis_Modes::CONT  = CMTOA / Constants::TWOPI;
// NOTE: Original TWOPI was 6.2832, results in small roundoff diffs from ptraj

/// Strings describing modes calculation types.
const char* Analysis_Modes::analysisTypeString[] = {
  "rms fluctuations",
  "displacements",
  "correlation functions",
  "coordinate projection",
  "eigenvalue fraction",
  "root mean square inner product"
};

// DESTRUCTOR
Analysis_Modes::~Analysis_Modes() {
  if (tOutParm_ != 0)
    delete tOutParm_;
}

// Analysis_Modes::CheckDeprecated()
void Analysis_Modes::CheckDeprecated(ArgList& analyzeArgs, std::string& modesname, 
                                     const char* key) {
  std::string arg = analyzeArgs.GetStringKey( key );
  if (!arg.empty()) {
    mprintf("Warning: Argument '%s' is deprecated, use 'name <modes>' instead.\n",
            key);
    if (modesname.empty()) modesname = arg;
  }
}

// Analysis_Modes::Setup()
Analysis::RetType Analysis_Modes::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  debug_ = debugIn;
  // Analysis type
  if (analyzeArgs.hasKey("fluct"))
    type_ = FLUCT;
  else if (analyzeArgs.hasKey("displ"))
    type_ = DISPLACE;
  else if (analyzeArgs.hasKey("corr"))
    type_ = CORR;
  else if (analyzeArgs.Contains("trajout"))
    type_ = TRAJ;
  else if (analyzeArgs.hasKey("eigenval"))
    type_ = EIGENVAL;
  else if (analyzeArgs.hasKey("rmsip"))
    type_ = RMSIP;
  else {
    mprinterr("Error: No analysis type specified.\n");
    return Analysis::ERR;
  }

  // Get modes name
  std::string modesfile = analyzeArgs.GetStringKey("name");
  if (modesfile.empty()) {
    // Check for deprecated args
    CheckDeprecated(analyzeArgs, modesfile, "file");
    CheckDeprecated(analyzeArgs, modesfile, "stack");
    if (modesfile.empty()) {
      mprinterr("Error: No 'name <modes data set name>' argument given.\n");
      return Analysis::ERR;
    }
  }
  // Get second modes name for RMSIP
  std::string modesfile2 = analyzeArgs.GetStringKey("name2");
  if (type_ == RMSIP) {
    if (modesfile2.empty()) {
      mprinterr("Error: 'rmsip' requires second modes data 'name2 <modes>'\n");
      return Analysis::ERR;
    }
  } else
    modesfile2.clear();
  // Get topology for TRAJ/CORR
  Topology* analyzeParm = setup.DSL().GetTopology( analyzeArgs ); 

  if (type_ == TRAJ ) {
    // Get trajectory format args for projected traj
    beg_ = analyzeArgs.getKeyInt("beg",1) - 1; // Args start at 1
    std::string tOutName = analyzeArgs.GetStringKey("trajout");
    if (tOutName.empty()) {
      mprinterr("Error: Require output trajectory filename, 'trajout <name>'\n");
      return Analysis::ERR;
    }
    TrajectoryFile::TrajFormatType tOutFmt = TrajectoryFile::UNKNOWN_TRAJ;
    if ( analyzeArgs.Contains("trajoutfmt") )
      tOutFmt = TrajectoryFile::WriteFormatFromString( analyzeArgs.GetStringKey("trajoutfmt"),
                                                       TrajectoryFile::AMBERTRAJ );
    if (analyzeParm == 0) {
      mprinterr("Error: Could not get topology for output trajectory.\n");
      return Analysis::ERR;
    }
    AtomMask tOutMask( analyzeArgs.GetStringKey("trajoutmask") );
    if ( analyzeParm->SetupIntegerMask( tOutMask ) || tOutMask.None() ) {
      mprinterr("Error: Could not setup output trajectory mask.\n");
      return Analysis::ERR;
    }
    tOutMask.MaskInfo();
    // Strip topology to match mask.
    if (tOutParm_ != 0) delete tOutParm_;
    tOutParm_ = analyzeParm->modifyStateByMask( tOutMask );
    if (tOutParm_ == 0) {
      mprinterr("Error: Could not create topology to match mask.\n");
      return Analysis::ERR;
    }
    // Setup output traj
    if (trajout_.InitTrajWrite( tOutName, ArgList(), tOutFmt ) != 0) {
      mprinterr("Error: Could not init output trajectory.\n");
      return Analysis::ERR;
    }
    // Get min and max for PC
    pcmin_ = analyzeArgs.getKeyDouble("pcmin", -10.0);
    pcmax_ = analyzeArgs.getKeyDouble("pcmax",  10.0);
    if (pcmax_ < pcmin_ || pcmax_ - pcmin_ < Constants::SMALL) {
      mprinterr("Error: pcmin must be less than pcmax\n");
      return Analysis::ERR;
    }
    tMode_ = analyzeArgs.getKeyInt("tmode", 1);
  } else {
    // Args for everything else
    beg_ = analyzeArgs.getKeyInt("beg",7) - 1; // Args start at 1
    bose_ = analyzeArgs.hasKey("bose");
    calcAll_ = analyzeArgs.hasKey("calcall");
  }
  end_ = analyzeArgs.getKeyInt("end", 50);
  factor_ = analyzeArgs.getKeyDouble("factor",1.0);
  std::string setname = analyzeArgs.GetStringKey("setname");

  // Check if modes name exists on the stack
  modinfo_ = (DataSet_Modes*)setup.DSL().FindSetOfType( modesfile, DataSet::MODES );
  if (modinfo_ == 0) {
    mprinterr("Error: '%s' not found: %s\n", modesfile.c_str(), DataSet_Modes::DeprecateFileMsg);
    return Analysis::ERR;
  }
  if (!modesfile2.empty()) {
    modinfo2_ = (DataSet_Modes*)setup.DSL().FindSetOfType( modesfile2, DataSet::MODES );
    if (modinfo2_ == 0) {
      mprinterr("Error: Set %s not found.\n", modesfile2.c_str());
      return Analysis::ERR;
    }
  }

  // Check modes type for specified analysis
  if (type_ == FLUCT || type_ == DISPLACE || type_ == CORR || type_ == TRAJ) {
    if (modinfo_->Meta().ScalarType() != MetaData::COVAR && 
        modinfo_->Meta().ScalarType() != MetaData::MWCOVAR)
    {
      mprinterr("Error: Modes must be of type COVAR or MWCOVAR for %s.\n",
                analysisTypeString[type_]);
      return Analysis::ERR;
    }
  }

  // Get output filename for types that use DataSets
  std::string outfilename = analyzeArgs.GetStringKey("out"); // TODO all datafile?
  DataFile* dataout = 0;
  if (type_ == FLUCT || type_ == DISPLACE || type_ == EIGENVAL || type_ == RMSIP)
    dataout = setup.DFL().AddDataFile( outfilename, analyzeArgs );
  else if (type_ == CORR) {
    // CORR-specific setup
    outfile_ = setup.DFL().AddCpptrajFile( outfilename, "Modes analysis",
                                           DataFileList::TEXT, true );
    if (outfile_ == 0) return Analysis::ERR;
    // Get list of atom pairs
    if (analyzeParm == 0) {
      mprinterr("Error: 'corr' requires topology.\n");
      return Analysis::ERR;
    }
    std::string maskexp = analyzeArgs.GetStringKey("mask1");
    if (maskexp.empty()) {
      while (analyzeArgs.hasKey("maskp")) {
        // Next two arguments should be one-atom masks
        std::string a1mask = analyzeArgs.GetMaskNext();
        std::string a2mask = analyzeArgs.GetMaskNext();
        if (a1mask.empty() || a2mask.empty()) {
          mprinterr("Error: For 'corr' two 1-atom masks are expected.\n");
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
          mprinterr("Error: For 'corr', masks should specify only one atom.\n"
                    "\tM1[%s]=%i atoms, M2[%s]=%i atoms.\n", m1.MaskString(), m1.Nselected(),
                    m2.MaskString(), m2.Nselected());
          return Analysis::ERR;
        }
      }
    } else {
      AtomMask mask1( maskexp );
      maskexp = analyzeArgs.GetStringKey("mask2");
      if (maskexp.empty()) {
        mprinterr("Error: 'mask2' must be specified if 'mask1' is.\n");
        return Analysis::ERR;
      }
      AtomMask mask2( maskexp );
      if ( analyzeParm->SetupIntegerMask( mask1 ) ) return Analysis::ERR;
      if ( analyzeParm->SetupIntegerMask( mask2 ) ) return Analysis::ERR;
      mask1.MaskInfo();
      mask2.MaskInfo();
      if (mask1.None() || mask2.None()) {
        mprinterr("Error: One or both masks are empty.\n");
        return Analysis::ERR;
      }
      if (mask1.Nselected() != mask2.Nselected()) {
        mprinterr("Error: # atoms in mask 1 not equal to # atoms in mask 2.\n");
        return Analysis::ERR;
      }
      for (int idx = 0; idx != mask1.Nselected(); idx++)
        atompairStack_.push_back( std::pair<int,int>( mask1[idx], mask2[idx] ) );
    }
    if ( atompairStack_.empty() ) {
      mprinterr("Error: No atom pairs found (use 'maskp' or 'mask1'/'mask2' keywords.)\n");
      return Analysis::ERR;
    }
  }

  // Set up data sets
  Dimension Xdim;
  if (type_ == FLUCT) {
    if (setname.empty()) setname = setup.DSL().GenerateDefaultName("FLUCT");
    MetaData md(setname, "rmsX");
    OutSets_.resize( 4, 0 );
    OutSets_[RMSX] = setup.DSL().AddSet( DataSet::DOUBLE, md );
    md.SetAspect("rmsY");
    OutSets_[RMSY] = setup.DSL().AddSet( DataSet::DOUBLE, md );
    md.SetAspect("rmsZ");
    OutSets_[RMSZ] = setup.DSL().AddSet( DataSet::DOUBLE, md );
    md.SetAspect("rms");
    OutSets_[RMS]  = setup.DSL().AddSet( DataSet::DOUBLE, md );
    Xdim = Dimension(1, 1, "Atom_no.");
  } else if (type_ == DISPLACE) {
    if (setname.empty()) setname = setup.DSL().GenerateDefaultName("DISPL");
    MetaData md(setname, "displX");
    OutSets_.resize( 3, 0 );
    OutSets_[RMSX] = setup.DSL().AddSet( DataSet::DOUBLE, md );
    md.SetAspect("displY");
    OutSets_[RMSY] = setup.DSL().AddSet( DataSet::DOUBLE, md );
    md.SetAspect("displZ");
    OutSets_[RMSZ] = setup.DSL().AddSet( DataSet::DOUBLE, md );
    Xdim = Dimension(1, 1, "Atom_no.");
  } else if (type_ == EIGENVAL) {
    if (setname.empty()) setname = setup.DSL().GenerateDefaultName("XEVAL");
    MetaData md(setname, "Frac");
    OutSets_.resize( 3, 0 );
    OutSets_[0] = setup.DSL().AddSet( DataSet::DOUBLE, md );
    md.SetAspect("Cumulative");
    OutSets_[1] = setup.DSL().AddSet( DataSet::DOUBLE, md );
    md.SetAspect("Eigenval");
    OutSets_[2] = setup.DSL().AddSet( DataSet::DOUBLE, md );
    Xdim = Dimension( 1, 1, "Mode" );
  } else if (type_ == RMSIP) {
    if (setname.empty()) setname = setup.DSL().GenerateDefaultName("RMSIP");
    OutSets_.push_back( setup.DSL().AddSet( DataSet::DOUBLE, setname ) );
    if (dataout != 0) dataout->ProcessArgs("noxcol");
    OutSets_[0]->SetupFormat() = TextFormat(TextFormat::GDOUBLE);
    OutSets_[0]->SetLegend( modinfo_->Meta().Legend() + "_X_" + modinfo2_->Meta().Legend() );
  }
  for (std::vector<DataSet*>::const_iterator set = OutSets_.begin(); set != OutSets_.end(); ++set)
  {
    if (*set == 0) return Analysis::ERR;
    if (dataout != 0) dataout->AddDataSet( *set );
    (*set)->SetDim(Dimension::X, Xdim);
  }

  // Status
  mprintf("    ANALYZE MODES: Calculating %s using modes from %s", 
          analysisTypeString[type_], modinfo_->legend());
  if ( type_ != TRAJ ) {
    if (type_ != EIGENVAL)
      mprintf(", modes %i to %i", beg_+1, end_);
    if (outfile_ != 0)
      mprintf("\n\tResults are written to %s\n", outfile_->Filename().full());
    else if (dataout != 0)
      mprintf("\n\tResults are written to '%s'\n", dataout->DataFilename().full());
    if (type_ != EIGENVAL && type_ != RMSIP) {
      if (bose_)
        mprintf("\tBose statistics used.\n");
      else
        mprintf("\tBoltzmann statistics used.\n");
      if (calcAll_)
        mprintf("\tEigenvectors associated with zero or negative eigenvalues will be used.\n");
      else
        mprintf("\tEigenvectors associated with zero or negative eigenvalues will be skipped.\n");
    }
    if (type_ == DISPLACE)
      mprintf("\tFactor for displacement: %f\n", factor_);
    if (type_ == CORR) {
      mprintf("\tUsing the following atom pairs:");
      for (modestack_it apair = atompairStack_.begin();
                        apair != atompairStack_.end(); ++apair)
        mprintf(" (%i,%i)", apair->first+1, apair->second+1 );
      mprintf("\n");
    }
    if (type_ == RMSIP)
      mprintf("\tRMSIP calculated to modes in %s\n", modinfo2_->legend());
  } else {
    mprintf("\n\tCreating trajectory for mode %i\n"
              "\tWriting to trajectory %s\n"
              "\tPC range: %f to %f\n"
              "\tScaling factor: %f\n", tMode_, 
            trajout_.Traj().Filename().full(), pcmin_, pcmax_, factor_);
  }

  return Analysis::OK;
}

// Analysis_Modes::Analyze()
Analysis::RetType Analysis_Modes::Analyze() {
  // Check # of modes
  if (type_ != TRAJ && type_ != EIGENVAL) {
    if (beg_ < 0 || beg_ >= modinfo_->Nmodes()) {
      mprinterr("Error: 'beg %i' is out of bounds.\n", beg_+1);
      return Analysis::ERR;
    }
    if (end_ > modinfo_->Nmodes()) {
      mprintf("Warning: 'end %i' is > # of modes, setting to %i\n",
              end_, modinfo_->Nmodes());
      end_ = modinfo_->Nmodes();
    }
    if (end_ <= beg_) {
      mprinterr("Warning: beg must be <= end, (%i -- %i)\n", beg_+1, end_);
      return Analysis::ERR;
    }
  }
  mprintf("\tModes '%s'", modinfo_->legend());
  if (modinfo_->EvalsAreFreq())
    mprintf(", eigenvalues are in cm^-1");
  if (modinfo_->EvecsAreMassWtd())
    mprintf(", eigenvectors are mass-weighted");
  mprintf("\n");
  if (!modinfo_->EvalsAreFreq() && type_ == CORR)
    mprintf("Warning: 'corr' analysis expects eigenvalues in cm^-1.\n");

  int err = 0;
  switch ( type_ ) {
    case FLUCT:    CalcFluct( *modinfo_ ); break;
    case DISPLACE: CalcDisplacement( *modinfo_ ); break;
    case CORR:     CalcDipoleCorr( *modinfo_ ); break;
    case TRAJ:     err = ProjectCoords( *modinfo_ ); break;
    case EIGENVAL: CalcEvalFrac( *modinfo_ ); break;
    case RMSIP:    err = CalcRMSIP( *modinfo_, *modinfo2_ ); break;
  }
  if (err != 0) return Analysis::ERR;

  return Analysis::OK;
}

// Analysis_Modes::CalcFluct()
void Analysis_Modes::CalcFluct(DataSet_Modes const& modes) {
  int natoms = modes.NavgCrd() / 3; // Possible because COVAR/MWCOVAR
  for (int vi = 0; vi < natoms; ++vi) {
    double sumx = 0.0;
    double sumy = 0.0;
    double sumz = 0.0;
    // Loop over all eigenvector elements for this atom
    const double* Vec = modes.Eigenvector(beg_) + vi*3;
    for (int mode = beg_; mode < end_; ++mode) {
      if (modes.EvalsAreFreq()) {
        double frq = modes.Eigenvalue(mode);
        // Don't use eigenvectors associated with zero or negative eigenvalues
        if (frq >= 0.5 || calcAll_) {
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
      } else {
        double ev = modes.Eigenvalue(mode);
        if ( ev > 0.0 || calcAll_) {
          double distx = Vec[0] * Vec[0];
          double disty = Vec[1] * Vec[1];
          double distz = Vec[2] * Vec[2];
          sumx += distx * ev;
          sumy += disty * ev;
          sumz += distz * ev;
        }
      }
      Vec += modes.VectorSize();
    }
    double c1 = 1.0;
    double c2 = 1.0;
    if (modes.EvalsAreFreq()) {
      c1 = CNST;
      c2 = CONT;
    }
    sumx *= c1;
    sumy *= c1;
    sumz *= c1;
    double rmsX = sqrt(sumx) * c2;
    OutSets_[RMSX]->Add( vi, &rmsX );
    double rmsY = sqrt(sumy) * c2;
    OutSets_[RMSY]->Add( vi, &rmsY );
    double rmsZ = sqrt(sumz) * c2;
    OutSets_[RMSZ]->Add( vi, &rmsZ );
    double rms = sqrt(sumx + sumy + sumz) * c2;
    OutSets_[RMS]->Add( vi, &rms );
  }
}

// Analysis_Modes::CalcDisplacement()
void Analysis_Modes::CalcDisplacement( DataSet_Modes const& modes ) {
  std::vector<double> results( modes.NavgCrd(), 0.0 );
  double dfactor = factor_;
  if (modes.EvalsAreFreq())
    dfactor *= (sqrt(CNST) * CONT);
  // Loop over all modes
  for (int mode = beg_; mode < end_; ++mode) {
    double fi;
    if (modes.EvalsAreFreq()) {
      double frq = modes.Eigenvalue(mode);
      // Don't use eigenvectors associated with zero or negative eigenvalues
      if (frq < 0.5 && !calcAll_) continue;
      fi = 1.0 / frq;
      if (bose_) {
        double argq = CONSQ * frq;
        fi *= (fi * argq / tanh(argq));
        fi = sqrt(fi);
      }
    } else {
      fi = modes.Eigenvalue(mode);
      if (fi < Constants::SMALL && !calcAll_) continue;
      fi = sqrt( fi );
    }
    fi *= dfactor; // * CONT * factor
    // Loop over all vector elements
    const double* Vec = modes.Eigenvector(mode);
    for (int vi = 0; vi < modes.NavgCrd(); vi += 3) {
      results[vi  ] += Vec[vi  ] * fi;
      results[vi+1] += Vec[vi+1] * fi;
      results[vi+2] += Vec[vi+2] * fi;
    }
  }
  int anum = 0;
  for (int i3 = 0; i3 < modes.NavgCrd(); i3 += 3, anum++) {
    OutSets_[RMSX]->Add( anum, &results[i3] );
    OutSets_[RMSY]->Add( anum, &results[i3+1] );
    OutSets_[RMSZ]->Add( anum, &results[i3+2] );
  }
}

// Analysis_Modes::CalcDipoleCorr()
void Analysis_Modes::CalcDipoleCorr(DataSet_Modes const& modes) {
  double qcorr = 1.0; // For when bose is false
  int rsize = atompairStack_.size() * (end_ - beg_ + 1);
  std::vector<double> results( rsize, 0.0 );
  // Loop over atom pairs 
  std::vector<double>::iterator Res = results.begin();
  DataSet_Modes::Darray const& Avg = modes.AvgCrd();
  for (modestack_it apair = atompairStack_.begin(); apair != atompairStack_.end(); ++apair)
  {
    int idx1 = apair->first  * 3;
    int idx2 = apair->second * 3;
    // Check if coordinates are out of bounds
    if (idx1 >= modes.NavgCrd() || idx2 >= modes.NavgCrd()) {
      mprintf("Warning: Atom pair %i -- %i is out of bounds (# avg atoms in modes = %i)\n",
              apair->first + 1, apair->second + 1, modes.NavgCrd() / 3);
      continue;
    }
    // Calc unit vector along at2->at1 bond
    Vec3 vec = Vec3(Avg[idx1], Avg[idx1+1], Avg[idx1+2]) -
               Vec3(Avg[idx2], Avg[idx2+1], Avg[idx2+2]);
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
      double eval = modes.Eigenvalue(mode);
      if (eval >= 0.5 || calcAll_) {
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
        const double* Evec = modes.Eigenvector(mode);
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
  outfile_->Printf("#Analysis of modes: CORRELATION FUNCTIONS\n");
  outfile_->Printf("%-10s %10s %10s %10s %10s %10s\n", "#Atom1", "Atom2", "Mode", 
                 "Freq", "1-S^2", "P2(cum)");
  int ncnt = 0;
  for (modestack_it apair = atompairStack_.begin();
                    apair != atompairStack_.end(); ++apair)
  {
    outfile_->Printf("%10i %10i\n", apair->first+1, apair->second+1);
    double val = 1.0;
    for (int mode = beg_; mode < end_; ++mode) {
      double frq = modes.Eigenvalue(mode);
      if (frq >= 0.5 || calcAll_) {
        val += results[ncnt];
        outfile_->Printf("%10s %10s %10i %10.5f %10.5f %10.5f\n",
                       "", "", mode, frq, results[ncnt], val);
        ++ncnt;
      }
    }
  } // END loop over CORR atom pairs
}

/** Calculate projection of coords along given mode. */
void Analysis_Modes::CalculateProjection(int set, Frame const& Crd, int mode) const {
  double proj = 0.0;
  DataSet_Modes::Darray const& Avg = modinfo_->AvgCrd();
  const double* Vec = modinfo_->Eigenvector(mode);
  for (int idx = 0; idx < Crd.size(); ++idx)
    proj += (Crd[idx] - Avg[idx]) * 1.0 * Vec[idx];
  mprintf("\tFrame %i mode %i projection = %f\n", set, mode, proj);
}

/** Project average coords along eigenvectors */
int Analysis_Modes::ProjectCoords(DataSet_Modes const& modes) {
  double scale = factor_;
  int max_it = (int)(pcmax_ - pcmin_);
  // Check that size of eigenvectors match # coords
  int ncoord = tOutParm_->Natom() * 3;
  if (ncoord != modes.NavgCrd()) {
    mprinterr("Error: # selected coords (%i) != eigenvector size (%i)\n",
               ncoord, modes.NavgCrd());
    return 1;
  }
  // Check that mode is valid.
  if (tMode_ < 1 || tMode_ > modes.Nmodes() ) {
    mprinterr("Error: mode %i is out of bounds.\n", tMode_);
    return Analysis::ERR;
  }
  // Setup frame to hold output coords, initalized to avg coords.
  Frame outframe;
  outframe.SetupFrameXM( modes.AvgCrd(), modes.Mass() );
  // Point to correct eigenvector
  const double* Vec = modes.Eigenvector(tMode_-1);
  // Initialize coords to pcmin
  for (int idx = 0; idx < ncoord; idx++)
    outframe[idx] += pcmin_ * Vec[idx];
  if (debug_>0) CalculateProjection(0, outframe, tMode_-1);
  if (trajout_.SetupTrajWrite( tOutParm_, CoordinateInfo(), max_it )) {
    mprinterr("Error: Could not open output modes traj '%s'\n", trajout_.Traj().Filename().full());
    return 1;
  }
  // Write first frame with coords at pcmin.
  int set = 0;
  trajout_.WriteSingle(set++, outframe);
  // Main loop
  for (int it = 0; it < max_it; it++) {
    double* crd = outframe.xAddress();
    // Move coordinates along eigenvector.
    for (int idx = 0; idx < ncoord; ++idx)
      crd[idx] += scale * Vec[idx];
    if (debug_>0) CalculateProjection(set, outframe, tMode_-1);
    // DEBUG: calc PC projection for first 3 modes
    //for (int m = 0; m < 3; m++)
    //  CalculateProjection(set, outframe, m);
    // Write frame
    trajout_.WriteSingle(set++, outframe);
  }
  trajout_.EndTraj();
  return 0;
}

// Analysis_Modes::CalcEvalFrac()
void Analysis_Modes::CalcEvalFrac(DataSet_Modes const& modes) {
  double sum = 0.0;
  for (unsigned int mode = 0; mode != modes.Size(); mode++)
    sum += modes.Eigenvalue( mode );
  mprintf("\t%zu eigenvalues, sum is %f\n", modes.Size(), sum);
  double cumulative = 0.0;
  for (unsigned int mode = 0; mode != modes.Size(); mode++) {
    double frac = modes.Eigenvalue( mode ) / sum;
    cumulative += frac;
    OutSets_[0]->Add( mode, &frac );
    OutSets_[1]->Add( mode, &cumulative );
    frac = modes.Eigenvalue( mode );
    OutSets_[2]->Add( mode, &frac );
  }
}

// Analysis_Modes::CalcRMSIP()
int Analysis_Modes::CalcRMSIP(DataSet_Modes const& modes1, DataSet_Modes const& modes2) {
  if (modes1.VectorSize() != modes2.VectorSize()) {
    mprinterr("Error: '%s' vector size (%i) != '%s' vector size (%i)\n",
              modes1.legend(), modes1.VectorSize(),
              modes2.legend(), modes2.VectorSize());
    return 1;
  }
  if ( beg_ >= modes2.Nmodes() || end_ > modes2.Nmodes() ) {
    mprinterr("Error: beg/end out of range for %s (%i modes)\n",
              modes2.legend(), modes2.Nmodes());
    return 1;
  }
  double sumsq = 0.0;
  if (modes1.Meta().ScalarType() != modes2.Meta().ScalarType())
    mprintf("Warning: Modes '%s' type (%s) does not match modes '%s' type (%s)\n"
            "Warning; RMSIP value may not make sense.\n",
            modes1.legend(), modes1.Meta().TypeString(),
            modes2.legend(), modes2.Meta().TypeString());
  if (modes1.EvecsAreMassWtd() && modes1.Mass().empty()) {
     mprintf("Warning: Modes '%s' have been mass-weighted but no mass information present.\n",
             modes1.legend());
  }
  if (modes2.EvecsAreMassWtd() && modes2.Mass().empty()) {
     mprintf("Warning: Modes '%s' have been mass-weighted but no mass information present.\n",
             modes2.legend());
  }

  for (int m1 = beg_; m1 < end_; m1++) {
    const double* ev1 = modes1.Eigenvector(m1);
    for (int m2 = beg_; m2 < end_; m2++) {
      const double* ev2 = modes2.Eigenvector(m2);
      double dot = 0.0;
      for (int iv = 0; iv < modes1.VectorSize(); iv++)
        dot += ev1[iv] * ev2[iv];
      sumsq += (dot * dot);
    }
  }
  sumsq /= (double)(end_ - beg_);
  double rmsip = sqrt( sumsq );
  OutSets_[0]->Add( 0, &rmsip );
  return 0;
}
