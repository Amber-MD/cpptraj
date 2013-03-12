#include <cstring> // memset
#include "Analysis_IRED.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // fileExists
#include "Constants.h" // PI
#include "ComplexArray.h"

// CONSTRUCTOR
Analysis_IRED::Analysis_IRED() :
  freq_(-1),
  tstep_(1),
  tcorr_(10000),
  distnh_(1.02),
  order_(2),
  debug_(0),
  relax_(false),
  norm_(false),
  drct_(false),
  data1_(0),
  cf_(0),
  cf_cjt_(0),
  cfinf_(0),
  taum_(0),
  modinfo_(0)
{}

void Analysis_IRED::Help() {
  mprintf("\t[relax freq <hz> [NHdist <distnh>]] [order <order>]\n");
  mprintf("\ttstep <tstep> tcorr <tcorr> out <filename> [norm] [drct]\n");
  mprintf("\tmodes <modesname> [beg <ibeg> end <iend>]\n");
  mprintf("\tPerform isotropic reorientational Eigenmode dynamics analysis.\n");
}

// DESTRUCTOR
Analysis_IRED::~Analysis_IRED() {
  if (cf_!=0) delete[] cf_;
  if (cf_cjt_!=0) delete[] cf_cjt_;
  if (cfinf_!=0) delete[] cfinf_;
  if (taum_!=0) delete[] taum_;
}

// Analysis_IRED::Setup()
Analysis::RetType Analysis_IRED::Setup(ArgList& analyzeArgs, DataSetList* DSLin,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  debug_ = debugIn;
  int ibeg=0, iend=0;
  // Count and store the number of previously defined IRED vectors.
  for ( DataSetList::const_iterator DS = DSLin->begin(); DS != DSLin->end(); ++DS) {
    if ( (*DS)->Type() == DataSet::VECTOR ) {
      DataSet_Vector* Vtmp = (DataSet_Vector*)(*DS);
      if (Vtmp->IsIred())
        IredVectors_.push_back( Vtmp );
    }
  }
  if (IredVectors_.empty()) {
    mprinterr("Error: analyze ired: corrired: No IRED vectors defined.\n");
    return Analysis::ERR;
  }
  // Get order for Legendre polynomial
  order_ = analyzeArgs.getKeyInt("order",2);
  if (order_ < 0 || order_ > 2) {
    mprintf("Warning: vector order out of bounds (<0 or >2), resetting to 2.\n");
    order_ = 2;
  }
  // Get modes name
  std::string modesfile = analyzeArgs.GetStringKey("modes");
  if (modesfile.empty()) {
    mprinterr("Error: analyze ired: No 'modes <name>' arg given, needed for 'corrired'.\n");
    return Analysis::ERR;
  }
  // Check if modes name exists on the stack
  bool modesFromFile = false;
  modinfo_ = (DataSet_Modes*)DSLin->FindSetOfType( modesfile, DataSet::MODES );
  if (modinfo_ == 0) {
    // If not on stack, check for file.
    if ( fileExists(modesfile.c_str()) ) {
      ibeg = analyzeArgs.getKeyInt("beg",1);
      iend = analyzeArgs.getKeyInt("end", 50);
      modinfo_ = (DataSet_Modes*)DSLin->AddSet( DataSet::MODES, modesfile, "Modes" );
      if (modinfo_->ReadEvecFile( modesfile, ibeg, iend )) return Analysis::ERR;
      modesFromFile = true;
    }
  }
  if (modinfo_ == 0) {
    mprinterr("Error: analyze ired: Modes %s DataSet/file not found.\n",modesfile.c_str());
    return Analysis::ERR;
  }
  // TODO: Check that number of evecs match number of IRED vecs
  orderparamfile_ = analyzeArgs.GetStringKey("orderparamfile");

  // Get tstep, tcorr, filenames
  tstep_ = analyzeArgs.getKeyDouble("tstep", 1.0);
  tcorr_ = analyzeArgs.getKeyDouble("tcorr", 10000.0);
  noeFilename_ = analyzeArgs.GetStringKey("noefile");
  filename_ = analyzeArgs.GetStringKey("out");
  if (filename_.empty()) {
    mprinterr("Error: analyze ired: No outfile given ('out <filename>').\n");
    return Analysis::ERR;
  }

  // Get norm, drct, relax
  norm_ = analyzeArgs.hasKey("norm");
  drct_ = analyzeArgs.hasKey("drct");
  relax_ = analyzeArgs.hasKey("relax");

  // Relax parameters
  if (relax_) {
    // Get freq, NH distance
    freq_ = analyzeArgs.getKeyDouble("freq", -1.0);
    if (freq_ == -1.0) {
      mprinterr("Error: analyze ired: No frequency for calculation of relaxation\n");
      mprinterr("       parameters given ('freq <frequency>').\n");
      return Analysis::ERR;
    }
    // 1.02 * 10**(-10) in Angstroms
    distnh_ = analyzeArgs.getKeyDouble("NHdist", 1.02);
  }

  // Print Status
  mprintf("    ANALYZE IRED: IRED calculation, %u vectors.\n", IredVectors_.size());
  if (!orderparamfile_.empty())
    mprintf("\tOrder parameters will be written to %s\n",orderparamfile_.c_str());
  mprintf("\tCorrelation time %lf, time step %lf\n", tcorr_, tstep_);
  mprintf("\tCorrelation functions are");
  if (norm_)
    mprintf(" normalized.\n");
  else
    mprintf(" not normalized.\n");
  mprintf("\tCorrelation functions are calculated using the");
  if (drct_)
    mprintf(" direct approach.\n");
  else
    mprintf(" FFT approach.\n");
  if (modesFromFile) 
    mprintf("\tIRED modes %i to %i are read from %s,\n", ibeg, iend, modesfile.c_str());
  else
    mprintf("\tIRED modes taken from DataSet %s\n", modinfo_->Legend().c_str());
  if (relax_) {
    mprintf("\t\tTauM, relaxation rates, and NOEs are calculated using the iRED\n");
    mprintf("\t\t  approach using an NH distance of %lf and a frequency of %lf\n",
            distnh_, freq_);
  }
  if (!noeFilename_.empty())
    mprintf("\t\tNOEs and relaxation rates will be written to %s\n",
            noeFilename_.c_str());
  mprintf("\t\tResults are written to %s\n", filename_.c_str());

  return Analysis::OK;
}

double Analysis_IRED::calc_spectral_density(int vi, double omega) {
  // Loop over all eigenvector elements vi for all modes
  double Jval = 0.0;
  const double* Vec = modinfo_->Eigenvectors() + vi;
  for (int mode = 0; mode < modinfo_->Nmodes(); ++mode) {
    Jval += (modinfo_->Eigenvalue(mode) * (*Vec * *Vec)) * 2.0 * taum_[mode] /
            (1.0 + omega*omega * taum_[mode]*taum_[mode]);
    Vec += modinfo_->VectorSize();
  }
  return Jval;
}           

// Analysis_IRED::Analyze()
Analysis::RetType Analysis_IRED::Analyze() {
  CorrF_FFT pubfft_;
  CorrF_Direct corfdir_;
  if (!orderparamfile_.empty()) {
    // Calculation of S2 order parameters according to 
    //   Prompers & Brüschweiler, JACS  124, 4522, 2002; 
    // Originally added by A.N. Koller & H. Gohlke.
    CpptrajFile orderout;
    if (orderout.OpenWrite(orderparamfile_)) {
      mprinterr("Error: Could not set up order parameter file.\n");
      return Analysis::ERR;
    }
    orderout.Printf("\n\t************************************\n");
    orderout.Printf("\t- Calculated iRed order parameters -\n");
    orderout.Printf("\t************************************\n\n");
    orderout.Printf("vector    S2\n----------------------\n");
    // Loop over all vector elements
    for (int vi = 0; vi < modinfo_->VectorSize(); ++vi) {
      // Sum according to Eq. A22 in Prompers & Brüschweiler, JACS 124, 4522, 2002
      double sum = 0.0;
      // Loop over all eigenvectors except the first five ones
      const double* evectorElem = modinfo_->Eigenvector(5) + vi;
      for (int mode = 5; mode < modinfo_->Nmodes(); ++mode) {
        sum += modinfo_->Eigenvalue(mode) * (*evectorElem) * (*evectorElem);
        evectorElem += modinfo_->VectorSize();
      }
      orderout.Printf(" %4i  %10.5f\n", vi, 1.0 - sum);
    }
    orderout.CloseFile();
  }

  if (modinfo_->Nmodes() != (int)IredVectors_.size()) {
    mprinterr("Error: analyze ired: # Modes in %s (%i) does not match # of Ired Vecs (%u)\n",
              modinfo_->Legend().c_str(), modinfo_->Nmodes(), IredVectors_.size());
    return Analysis::ERR;
  }

  // All IRED vectors must have the same size
  int Nframes_ = -1;
  for (std::vector<DataSet_Vector*>::iterator Vtmp = IredVectors_. begin();
                                              Vtmp != IredVectors_.end(); ++Vtmp)
  { 
    if (Nframes_ == -1)
      Nframes_ = (*Vtmp)->Size();
    else if (Nframes_ != (*Vtmp)->Size()) {
      mprinterr("Error: All IRED vectors must have the same size.\n");
      mprinterr("Error: Vector %s size = %i, first vector size = %i\n",
                (*Vtmp)->Legend().c_str(), (*Vtmp)->Size(), Nframes_);
      return Analysis::ERR;
    }
  }

  // Determine sizes
  int time = (int)(tcorr_ / tstep_) + 1;
  // nsteps
  int nsteps = 0;
  if (time > Nframes_)
    nsteps = Nframes_;
  else
    nsteps = time;
  // Allocate memory to hold complex numbers for direct or FFT
  if (drct_) {
    data1_.Allocate( Nframes_ );
    corfdir_.Allocate( nsteps );
  } else {
    // Initialize FFT
    pubfft_.Allocate( Nframes_ );
    data1_ = pubfft_.Array();
  }

  // -------------------- IRED CALCULATION ---------------------------
  // Store Modes Info
  int nvect = modinfo_->Nmodes();
  int nvectelem = modinfo_->VectorSize();
  const double* eigval = modinfo_->Eigenvalues();
  const double* vout = modinfo_->Eigenvectors();
  // Initialize memory
  cf_ = new double[ nvect * nsteps ];
  cf_cjt_ = new double[ nvect * nsteps ];
  for (int i = 0; i < nvect*nsteps; ++i) {
    cf_[i] = 0;
    cf_cjt_[i] = 0;
  }
  cfinf_ = new double[ nvect ];
  for (int i = 0; i < nvect; ++i)
    cfinf_[i] = 0;
  if (relax_) {
    taum_ = new double[ nvect ];
    for (int i = 0; i < nvect; ++i)
      taum_[i] = 0;
  }

  // Allocate memory to project spherical harmonics on eigenmodes
  int mtot = 2 * order_ + 1;
  int p2blocksize = 2 * mtot;                // Real + Img. for each -order <= m <= order
  int nsphereharm = Nframes_ * p2blocksize;  // Spherical Harmonics for each frame
  int ntotal = nvect * nsphereharm;          // Each vector has set of spherical harmonics
  double *cftmp1 = new double[ ntotal ];
  memset(cftmp1, 0, ntotal * sizeof(double));
  // Project spherical harmonics for each IRED vector on eigenmodes
  int n_ivec = 0;
  for (std::vector<DataSet_Vector*>::iterator ivec = IredVectors_.begin();
                                              ivec != IredVectors_.end(); ++ivec)
  {
    double* CF = cftmp1;
    (*ivec)->CalcSphericalHarmonics( order_ );
    // Loop over all eigenvectors
    for (int veci = 0; veci < nvect; ++veci) {
      double Qvec = vout[veci * nvectelem + n_ivec];
      // Loop over all m = -L, ...., L
      for (int midx = -order_; midx <= order_; ++midx) {
        // Loop over spherical harmonic coords for this m (Complex, [Real][Img])
        for ( int sidx = 2 * (midx + order_); sidx < nsphereharm; sidx += p2blocksize) {
          *(CF++) += (Qvec * (*ivec)->SphereHarm(sidx  ));
          *(CF++) += (Qvec * (*ivec)->SphereHarm(sidx+1));
        }
      }
    }
    ++n_ivec;
  }
  // cftmp1 now contains avgs over every IRED vector for each eigenvector:
  //   [m-2R0][m-2I0][m-2R1][m-2I1] ... [m-2RN][m-2IN][m-1R0][m-1I0] ... 
  
  // Loop over all eigenvectors
  double* CF = cftmp1;
  for (int veci = 0; veci < nvect; ++veci) {
    // Loop over all m = -L, ...., L
    for (int midx = -order_; midx <= order_; ++midx) {
      // Loop over all snapshots
      double cfinfavgreal = 0;
      double cfinfavgimg = 0;
      for (int k = 0; k < Nframes_*2; k += 2) {
        data1_[k  ]   = *CF;
        cfinfavgreal += *(CF++);
        data1_[k+1]   = *(CF);
        cfinfavgimg  += *(CF++);
        //mprintf("CDBG:\tVec=%i Frame=%i data1[%i]=%lf data1_[%i]=%lf\n",i,k,
        //        2*k, data1_[2*k], 2*k+1, data1_[2*k+1]);
      }
      cfinfavgreal /= Nframes_;
      cfinfavgimg  /= Nframes_;
      //mprintf("CDBG:\tVect[%i] cfinfavgreal=%lf cfinfavgimg=%lf\n",i,cfinfavgreal,cfinfavgimg);
      // Calc plateau value of correlation function (= C(m,t->T) in Bruschweiler paper (A20))
      cfinf_[veci] += (cfinfavgreal * cfinfavgreal) + (cfinfavgimg * cfinfavgimg);
      if (drct_) {
        // Calc correlation function (= C(m,l,t) in Bruschweiler paper) using direct approach
        corfdir_.AutoCorr( data1_ );
      } else {
        // Pad with zero's at the end
        data1_.PadWithZero( Nframes_ );
        // Calc correlation function (= C(m,l,t) in Bruschweiler paper) using FFT
        pubfft_.AutoCorr(data1_);
      }
      // Sum into cf (= C(m,t) in Bruschweiler paper)
      for (int k = 0; k < nsteps; ++k) {
        //mprintf("CDBG: cf[%i] += data1[%i] (%lf)\n",idx1+k, 2*k, data1_[2 * k]);
        cf_[nsteps * veci + k] += data1_[2 * k];
      }
    }
  }
  delete[] cftmp1;
  // Calculate correlation function for each vector:
  // Cj(t) according to eq. A23 in Prompers & Brüschweiler, JACS  124, 4522, 2002; 
  // added by A.N. Koller & H. Gohlke
  for (int i = 0; i < nvect; ++i) {
    for (int k = 0; k < nsteps; ++k) {
      double sum = 0;
      for (int j = 0; j < nvect; ++j) {
        //mprintf("CDBG: eigval[%i]=%lf vout[%i]=%lf cf[%i]=%lf frame-k=%i\n",
        //        j, eigval[j],
        //        j * nvectelem + i, vout[j * nvectelem + i],
        //        nsteps * j + k, cf_[nsteps * j + k], frame - k);
        sum += (eigval[j] * (vout[j * nvectelem + i] * vout[j * nvectelem + i])) *
               (cf_[nsteps * j + k] / (Nframes_ - k));
      }
      //mprintf("CDBG:\tVec=%i Step=%i sum=%lf\n",i,k,sum);
      cf_cjt_[nsteps * i + k] = sum;
    }
  }

  if (relax_) {
    // TauM calculation: eq. A19 from Prompers & Brüschweiler, JACS  124, 4522, 2002 is used.
    // Values are calculated from normalized correlation functions (Cm(t)s).
    // Added by Alrun N. Koller & H. Gohlke
    // Conversion to picoseconds
    double deltat = tstep_ * 1.0E-12;
    if (nvectelem != nvect) {
      mprinterr("Error: analyze ired: Different number of eigenmodes (%i) and\n", nvect);
      mprinterr("       eigenmode elements (%i)\n", nvectelem);
      return Analysis::ERR;
    }
    // consider only first third of the correlation curve to avoid fitting errors
    // NOTE: Done this way to be consistent with PTRAJ behavior. This really should
    //       be cast back to an integer.
    double new_nsteps = nsteps / 3.0;
    //mprintf("CDBG: new_nsteps= %lf\n",new_nsteps);
    for (int i = 0; i < nvect; ++i) {
      double integral = 0;
      double cfinfval = cfinf_[i] / cf_[nsteps * i] * Nframes_;
      //mprintf("CDBG: cfinfval= %.10lE\n",cfinfval);
      for ( int j = 0 ; j < new_nsteps; ++j ) {
          double cfk = cf_[nsteps*i + j] * Nframes_ / (cf_[nsteps * i] * (Nframes_ - j));
          double cfk1 = cf_[nsteps*i + j + 1 ] * Nframes_ / (cf_[nsteps * i] * (Nframes_ - (j+1)));
          double fa = cfk - cfinfval;
          double fb = cfk1 - cfinfval;
          integral += deltat * ( fa + fb ) * 0.5;
          //mprintf("CDBG:\tintegral[%i]= %.10lE\n",j,integral);
      }
      double taum_val = integral / ( 1.0 - cfinfval );
      taum_[i] = taum_val;
      //mprintf("CDBG: taum[%i]= %.10lE\n", i, taum_[i]);
    }

    // Relaxation calculation. Added by Alrun N. Koller & H. Gohlke
    CpptrajFile noefile;
    int err = 0;
    if (noeFilename_.empty())
      err = noefile.SetupWrite(0, debug_);
    else
      err = noefile.SetupWrite(noeFilename_, debug_);
    if (err != 0) {
      mprinterr("Error: analyze ired: Could not open NOE file for write.\n");
      return Analysis::ERR;
    }
    noefile.OpenFile();
    noefile.Printf("\n\t****************************************");
    noefile.Printf("\n\t- Calculated relaxation rates and NOEs -");
    noefile.Printf("\n\t****************************************\n\n");
    noefile.Printf("vector   %10s   %10s   %10s\n","R1","R2","NOE");
    // conversion from Angstrom to meter
    double rnh = distnh_ * 1.0E-10;
    // ---------- CONSTANTS ----------
    // in H m^-1; permeability
    const double mu_zero = 4.0 * PI * 1.0E-7;
    // in m^2 kg s1-1 ; Js ; Planck's constant
    const double ha = 6.626176 * 1.0E-34;
    // in rad s^-1 T^-1 ; gyromagnetic ratio; T = absolut temperature in K
    const double gamma_h = 2.6751987 * 1.0E8;
    // in rad s^-1 T^-1 ; gyromagnetic ratio
    const double gamma_n = -2.7126 * 1.0E7;
    const double csa = -170.0 * 1.0E-6;
    // conversion from MHz to rad s^-1
    //double spec_freq = freq_ * TWOPI * 1.0E6;
    // in T (Tesla)
    double b_zero = TWOPI * freq_ * 1.0E6 / gamma_h;
    // Next two both in rad s^-1
    const double lamfreqh = -1 * gamma_h * b_zero;
    const double lamfreqn = -1 * gamma_n * b_zero;
    double c2 = lamfreqn*lamfreqn * csa*csa;
    double d2 = (mu_zero * ha * gamma_n * gamma_h)/( 8.0 * (PI*PI) * (rnh*rnh*rnh));
    d2 = d2*d2; // fix from Junchao
    // -------------------------------
    // loop over all vector elements --> have only one element/vector here; nvectelem = nelem
    for (int i = 0; i < nvectelem; ++i) {
      // Eq. A1 in Prompers & Brüschweiler, JACS  124, 4522, 2002
      double R1 = d2 / 20.0 * (
            calc_spectral_density( i, lamfreqh - lamfreqn ) 
            + 3.0 * calc_spectral_density( i, lamfreqn ) 
            + 6.0 * calc_spectral_density( i, lamfreqh + lamfreqn)
          ) + c2 / 15.0 * calc_spectral_density( i, lamfreqn );
      // Eq. A2 in Prompers & Brüschweiler, JACS  124, 4522, 2002
      double R2 = d2 / 40.0 * ( 4.0 *
            calc_spectral_density( i, 0.0 ) 
            + calc_spectral_density( i, lamfreqh - lamfreqn )
            + 3.0 * calc_spectral_density( i, lamfreqn )
            + 6.0 * calc_spectral_density( i, lamfreqh ) 
            + 6.0 * calc_spectral_density( i, lamfreqh + lamfreqn)
          ) + c2 / 90.0 * ( 4.0 *
            calc_spectral_density( i, 0.0 ) 
            + 3.0 * calc_spectral_density( i, lamfreqn) );
      // Eq. A3 in Prompers & Brüschweiler, JACS  124, 4522, 2002
      double Tj = d2 / 20.0 * ( 6.0 *
            calc_spectral_density( i, lamfreqh + lamfreqn ) 
            - calc_spectral_density( i, lamfreqh + lamfreqn ) );

      double Noe = 1.0 + ( gamma_h * 1.0/gamma_n ) * ( 1.0 / R1 ) * Tj;
      noefile.Printf("%6i   %10.5f   %10.5f   %10.5f\n", i, R1, R2, Noe);
    }
    noefile.Printf("\n\n");
    noefile.CloseFile();
  } // END if (relax_)

  // ----- PRINT IRED -----
  // Setup and open files with .cjt/.cmt extensions.
  std::string cmtname = filename_ + ".cmt";
  CpptrajFile cmtfile;
  if (cmtfile.OpenWrite( cmtname )) return Analysis::ERR;
  std::string cjtname = filename_ + ".cjt";
  CpptrajFile cjtfile;
  if (cjtfile.OpenWrite( cjtname )) return Analysis::ERR;  
  // Print headers
  cmtfile.Printf("Auto-correlation functions Cm(t) for each eigenmode m, IRED type according to eq. A18 Prompers & Brüschweiler, JACS  124, 4522, 2002\n");
  cjtfile.Printf("Auto-correlation functions Cj(t) for each ired vector j, IRED type according to eq. A23 Prompers & Brüschweiler, JACS  124, 4522, 2002\n");
  cmtfile.Printf("%12s","XXX");
  cjtfile.Printf("%12s","XXX");
  int colwidth = 11;
  int tgti = 10;
  for (int i = 1; i <= nvect; ++i) {
    if (i == tgti) {
      --colwidth;
      if (colwidth < 7) colwidth = 7;
      tgti *= 10;
    }
    cmtfile.Printf("%*s%i", colwidth, "Mode",   i);
    cjtfile.Printf("%*s%i", colwidth, "Vector", i);
  }
  cmtfile.Printf("\n");
  cjtfile.Printf("\n");
  // Print cfinf
  cmtfile.Printf("%12s", "C(m,t->T)");
  for (int i = 0; i < nvect; ++i) {
    if (norm_)
      cmtfile.Printf("%12.8f", cfinf_[i] / cf_[nsteps * i] * Nframes_);
    else
      cmtfile.Printf("%12.8f", cfinf_[i]);
  }
  cmtfile.Printf("\n");
  // Print Relaxation
  if (relax_) {
    // Print Taum in ps
    cmtfile.Printf("%12s", "Tau_m [ps]");
    for (int i = 0; i < nvect; ++i)
      cmtfile.Printf("%12.6f", taum_[i]* 1.0E12);
    cmtfile.Printf("\n");
  }
  // Print cf
  for (int i = 0; i < nsteps; ++i) {
    cmtfile.Printf("%12.8f", (double)i * tstep_);
    cjtfile.Printf("%12.8f", (double)i * tstep_);
    for (int j = 0; j < nvect; ++j) {
      if (norm_) {
        cmtfile.Printf("%12.8f", cf_[nsteps*j + i] * Nframes_ / (cf_[nsteps * j] * (Nframes_ - i)));
        cjtfile.Printf("%12.8f", cf_cjt_[nsteps*j + i] / cf_cjt_[nsteps*j]);
      } else {
        // 4/5*PI due to spherical harmonics addition theorem
        cmtfile.Printf("%12.8f", FOURFIFTHSPI * cf_[nsteps*j + i] / (Nframes_ - i));
        cjtfile.Printf("%12.8f", FOURFIFTHSPI * cf_cjt_[nsteps*j + i]);
      }
    }
    cmtfile.Printf("\n");
    cjtfile.Printf("\n");
  }
  cmtfile.CloseFile();
  cjtfile.CloseFile();

  return Analysis::OK;
}
