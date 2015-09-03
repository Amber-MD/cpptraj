#include "Analysis_IRED.h"
#include "CpptrajStdio.h"
#include "Constants.h" // PI
#include "Corr.h"
#include "StringRoutines.h" // DEBUG
#include "DataSet_Mesh.h" // DEBUG

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
  cf_(0),
  cf_cjt_(0),
  cfinf_(0),
  taum_(0),
  orderout_(0),
  noefile_(0),
  cmtfile_(0),
  cjtfile_(0),
  data_s2_(0),
  modinfo_(0)
{}

void Analysis_IRED::Help() {
  mprintf("\t[relax freq <MHz> [NHdist <distnh>]] [order <order>]\n"
          "\ttstep <tstep> tcorr <tcorr> out <filename> [norm] [drct]\n"
          "\tmodes <modesname> [name <output sets name>]\n"
          "  Perform isotropic reorientational Eigenmode dynamics analysis.\n");
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
  // Count and store the number of previously defined IRED vectors.
  for ( DataSetList::const_iterator DS = DSLin->begin(); DS != DSLin->end(); ++DS) {
    if ( (*DS)->Type() == DataSet::VECTOR ) {
      DataSet_Vector* Vtmp = (DataSet_Vector*)(*DS);
      if (Vtmp->IsIred())
        IredVectors_.push_back( Vtmp );
    }
  }
  if (IredVectors_.empty()) {
    mprinterr("Error: No IRED vectors defined.\n");
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
    mprinterr("Error: No modes data specified: use 'modes <name>'.\n");
    return Analysis::ERR;
  }
  // Check if modes name exists on the stack
  modinfo_ = (DataSet_Modes*)DSLin->FindSetOfType( modesfile, DataSet::MODES );
  if (modinfo_ == 0) {
    mprinterr("Error: %s\n", DataSet_Modes::DeprecateFileMsg);
    return Analysis::ERR;
  }
  // Get tstep, tcorr, filenames
  tstep_ = analyzeArgs.getKeyDouble("tstep", 1.0);
  tcorr_ = analyzeArgs.getKeyDouble("tcorr", 10000.0);
  std::string orderParamName = analyzeArgs.GetStringKey("orderparamfile");
  if (!orderParamName.empty()) {
    orderout_ = DFLin->AddCpptrajFile(orderParamName, "IRED order parameters");
    if (orderout_ == 0) return Analysis::ERR;
  }
  noefile_ = DFLin->AddCpptrajFile(analyzeArgs.GetStringKey("noefile"), "IRED NOEs",
                                   DataFileList::TEXT, true);
  if (noefile_ == 0) return Analysis::ERR;
  std::string filename = analyzeArgs.GetStringKey("out");
  if (filename.empty()) {
    mprinterr("Error: No outfile given ('out <filename>').\n");
    return Analysis::ERR;
  }
  cmtfile_ = DFLin->AddCpptrajFile(filename + ".cmt", "Auto-correlation functions Cm(t)");
  cjtfile_ = DFLin->AddCpptrajFile(filename + ".cjt", "Auto-correlation functions Cj(t)");
  if (cmtfile_ == 0 || cjtfile_ == 0) return Analysis::ERR;
  // Output data sets
  std::string dsname = analyzeArgs.GetStringKey("name");
  if (dsname.empty()) dsname = DSLin->GenerateDefaultName("IRED");
  data_s2_ = DSLin->AddSetAspect(DataSet::FLOAT, dsname, "S2");
  if (data_s2_ == 0) return Analysis::ERR;

  // Get norm, drct, relax
  norm_ = analyzeArgs.hasKey("norm");
  drct_ = analyzeArgs.hasKey("drct");
  relax_ = analyzeArgs.hasKey("relax");

  // Relax parameters
  if (relax_) {
    // Get freq, NH distance
    freq_ = analyzeArgs.getKeyDouble("freq", -1.0);
    if (freq_ == -1.0) {
      mprinterr("Error: No frequency for calculation of relaxation\n"
                "Error:   parameters given ('freq <frequency>').\n");
      return Analysis::ERR;
    }
    // 1.02 * 10**(-10) in Angstroms
    distnh_ = analyzeArgs.getKeyDouble("NHdist", 1.02);
  }

  // Print Status
  mprintf("    IRED: %u IRED vectors.\n", IredVectors_.size());
  if (orderout_ != 0)
    mprintf("\tOrder parameters will be written to %s\n", orderout_->Filename().full());
  mprintf("\tCorrelation time %f, time step %lf\n", tcorr_, tstep_);
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
  mprintf("\tIRED modes will be taken from DataSet %s\n", modinfo_->legend());
  if (relax_) {
    mprintf("\t\tTauM, relaxation rates, and NOEs are calculated using the iRED\n"
            "\t\t  approach using an NH distance of %lf Ang. and a frequency of %lf MHz\n",
            distnh_, freq_);
  
    mprintf("\t\tNOEs and relaxation rates will be written to %s\n",
            noefile_->Filename().full());
  }
  mprintf("\t\tResults are written to %s and %s\n", cmtfile_->Filename().full(),
          cjtfile_->Filename().full());
  mprintf("#Citation: Prompers, J. J.; Brüschweiler, R.; \"General framework for\n"
          "#          studying the dynamics of folded and nonfolded proteins by\n"
          "#          NMR relaxation spectroscopy and MD simulation\"\n"
          "#          J. Am. Chem. Soc. (2002) V.124 pp.4522-4534\n");

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

// Calculate spectral density for given vector and omega.
double Analysis_IRED::Jw(int ivec, double omega, std::vector<double> TauM) const
{
  double Jval = 0.0;
  for (int mode = 0; mode != modinfo_->Nmodes(); ++mode)
  {
    // Get element ivec of current eigenvector
    double evectorElement = *(modinfo_->Eigenvector(mode) + ivec);
    Jval += ((modinfo_->Eigenvalue(mode) * evectorElement * evectorElement)) *
            (2.0 * TauM[mode]) / (1.0 + omega*omega * TauM[mode]*TauM[mode]);
  }
  return Jval;
}

// Analysis_IRED::Analyze()
Analysis::RetType Analysis_IRED::Analyze() {
  CorrF_FFT pubfft_;
  CorrF_Direct corfdir_;
  ComplexArray data1_;
  mprintf("\t'%s' has %zu modes.\n", modinfo_->legend(), modinfo_->Size());
  if ( modinfo_->Size() != IredVectors_.size() )
    mprintf("Warning: Number of IRED vectors (%zu) does not equal number of modes (%zu).\n",
            IredVectors_.size(), modinfo_->Size());
  if (orderout_ != 0) {
    // Calculation of S2 order parameters according to 
    //   Prompers & Brüschweiler, JACS  124, 4522, 2002; 
    // Originally added by A.N. Koller & H. Gohlke.
    orderout_->Printf("\n\t************************************\n"
                        "\t- Calculated iRed order parameters -\n"
                        "\t************************************\n\n"
                      "vector    S2\n----------------------\n");
    // Loop over all vector elements
    for (int vi = 0; vi < modinfo_->VectorSize(); ++vi) {
      // Sum according to Eq. A22 in Prompers & Brüschweiler, JACS 124, 4522, 2002
      double sum = 0.0;
      // Loop over all eigenvectors except the first five ones, i.e.
      // sum over all internal modes only.
      const double* evectorElem = modinfo_->Eigenvector(5) + vi;
      for (int mode = 5; mode < modinfo_->Nmodes(); ++mode) {
        sum += modinfo_->Eigenvalue(mode) * (*evectorElem) * (*evectorElem);
        evectorElem += modinfo_->VectorSize();
      }
      orderout_->Printf(" %4i  %10.5f\n", vi, 1.0 - sum);
      float fval = (float)(1.0 - sum);
      data_s2_->Add(vi, &fval);
    }
  }

  if (modinfo_->Nmodes() != (int)IredVectors_.size()) {
    mprinterr("Error: # Modes in %s (%i) does not match # of Ired Vecs (%u)\n",
              modinfo_->legend(), modinfo_->Nmodes(), IredVectors_.size());
    return Analysis::ERR;
  }

  // All IRED vectors must have the same size
  int Nframes_ = -1;
  for (std::vector<DataSet_Vector*>::const_iterator Vtmp = IredVectors_.begin();
                                                    Vtmp != IredVectors_.end(); ++Vtmp)
  { 
    if (Nframes_ == -1)
      Nframes_ = (*Vtmp)->Size();
    else if (Nframes_ != (int)(*Vtmp)->Size()) {
      mprinterr("Error: All IRED vectors must have the same size.\n"
                "Error:   Vector %s size = %i, first vector size = %i\n",
                (*Vtmp)->legend(), (*Vtmp)->Size(), Nframes_);
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
  // Ensure SH coords calculated for each ired vector.
  for (std::vector<DataSet_Vector*>::const_iterator iredvec = IredVectors_.begin();
                                                    iredvec != IredVectors_.end();
                                                  ++iredvec)
    (*iredvec)->CalcSphericalHarmonics( order_ );

  // Calculate Cm(t) for each mode
  std::vector<double> Plateau( modinfo_->Nmodes(), 0.0 );
  typedef std::vector< std::vector<double> > CmtArrayType;
  CmtArrayType CmtArray( modinfo_->Nmodes() );
  for (int mode = 0; mode != modinfo_->Nmodes(); mode++)
  {
    std::vector<double>& cm_t = CmtArray[mode];
    cm_t.assign( nsteps, 0.0 );
    // Loop over L = -order ... +order
    for (int Lval = -order_; Lval <= order_; Lval++)
    {
      // Values for determining plateau value of Cm(t)
      double plateau_r = 0.0;
      double plateau_i  = 0.0;
      // Loop over all frames.
      for (int frame = 0; frame != Nframes_; frame++)
      {
        const double* eigenvec = modinfo_->Eigenvector( mode );
        //DataSet_Modes::AvgIt Avg = modinfo_->AvgBegin(); //FIXME: Is this needed?
        int cidx = frame*2; // Index into complex arrays
        data1_[cidx  ] = 0.0; // Real component
        data1_[cidx+1] = 0.0; // Imaginary component
        // Loop over each eigenvector/ired vector element for this frame and L
        //for (int i = 0; i != modinfo_->VectorSize(); i++, Avg++)
        for (int i = 0; i != modinfo_->VectorSize(); i++)
        {
          ComplexArray const& SH_L = IredVectors_[i]->SphericalHarmonics( Lval );
          double alpha_r = SH_L[cidx    ] * eigenvec[i];
          plateau_r += alpha_r;
          data1_[cidx  ] += alpha_r;
          double alpha_i = SH_L[cidx + 1] * eigenvec[i];
          plateau_i += alpha_i;
          data1_[cidx+1] += alpha_i;
        }
      }
      plateau_r /= (double)Nframes_;
      plateau_i /= (double)Nframes_;
      // Calc contribution of this L to plateau value of correlation function 
      // (= C(m,t->T) in Bruschweiler paper (A20))
      Plateau[mode] += (plateau_r*plateau_r) + (plateau_i*plateau_i);
      // data1_ should now contain projected SH coords for this mode and L.
      // Calculate autocorrelation of projected coords.
      if (drct_)
        corfdir_.AutoCorr( data1_ );
      else {
        // Pad with zeros at the end
        data1_.PadWithZero( Nframes_ );
        pubfft_.AutoCorr( data1_ );
      }
      // Sum this L into Cm(t)
      for (int k = 0; k < nsteps; ++k)
        cm_t[k] += data1_[2 * k];
    }
    // Normalize
    // 4*PI / ((2*order)+1) due to spherical harmonics addition theorem
    double Kn = DataSet_Vector::SphericalHarmonicsNorm( order_ ); 
    // Normalize to 1.0
//    double Kn = (double)Nframes_ / cm_t[0];
    for (int k = 0; k < nsteps; ++k)
      cm_t[k] *= (Kn / (double)(Nframes_ - k));
//      cm_t[k] /= (double)(Nframes_ - k);
    // DEBUG: Write cm_t
    CpptrajFile cmt_out;
    cmt_out.OpenWrite( "dbg_cmt." + integerToString(mode) + ".dat" );
    cmt_out.Printf("%-12s Cm(t)_%i\n", "#Time", mode);
    //mprintf("DEBUG: Mode %i Kn= %g  cm_t[0]= %g\n", mode, Kn, cm_t[0]);
    for (int k = 0; k < nsteps; ++k)
      cmt_out.Printf("%12.4f %g\n", (double)k * tstep_, cm_t[k]);
      //cmt_out.Printf("%12.4f %g\n", (double)k * tstep_, cm_t[k] * Kn);
    cmt_out.CloseFile();
  }

  // Calculate tau_m for each mode.
  CpptrajFile tauout;
  tauout.OpenWrite("dbg_taum.dat");
  tauout.Printf("# Taum\n");
  std::vector<double> TauM( modinfo_->Nmodes(), 0.0 );
  for (int mode = 0; mode != modinfo_->Nmodes(); mode++)
  {
    std::vector<double> const& cm_t = CmtArray[mode];
    // Only integrate a third of the way in. 
    //int maxsteps = nsteps / 3;
    int maxsteps = nsteps;
    // Get previously calculated plateau value for this.
    double Cplateau = Plateau[mode];
    mprintf("Cplateau[%i]= %g\n", mode, Cplateau);
    // Integrate Cm(t) - Cplateau
    double sum = 0.0;
    for (int i = 1; i < maxsteps; i++)
    {
      //double b_minus_a = ((double)i * tstep_) - ((double)(i-1) * tstep_);
      double curr_val = cm_t[i  ] - Cplateau;
      double prev_val = cm_t[i-1] - Cplateau;
      sum += (tstep_ * (prev_val + curr_val) * 0.5);
    }
    TauM[mode] = sum / (cm_t[0] - Cplateau);
    tauout.Printf("%i %g\n", mode, TauM[mode]);
    // Convert tau from ps to s
    TauM[mode] *= 1.0E-12;
  }

/*
  // Create X mesh
  DataSet_Mesh mesh(nsteps, 0.0, tstep_ * nsteps);
  for (int mode = 0; mode != modinfo_->Nmodes(); mode++)
  {
    std::vector<double> const& cm_t = CmtArray[mode];
    for (int i = 0; i < nsteps; i++)
      mesh.SetY(i, cm_t[i]);
    double intercept, slope, corr;
    mesh.SingleExpRegression(slope, intercept, corr, false);
    TauM[mode] = -1.0 / slope;
    tauout.Printf("%i %g\n", mode, TauM[mode]);
    // Convert tau from ps to s
    TauM[mode] *= 1.0E-12;
  }
  tauout.CloseFile();
*/
  // Calculate Cj(t) for each vector j as weighted sum over Cm(t) arrays.
  // Cj(t) = SUM(m)[ dSjm^2 * Cm(t) ]
  // dSjm^2 = EVALm * (EVECm[i])^2
  // Cm(t) must be normalized to 1.0.
  for (unsigned int ivec = 0; ivec != IredVectors_.size(); ivec++)
  {
    std::vector<double> cj_t( nsteps, 0.0 );
    // Calculate dS^2 for this vector and mode.
    for (int mode = 0; mode != modinfo_->Nmodes(); ++mode)
    {
      std::vector<double> const& cm_t = CmtArray[mode];
      double Norm1 = 1.0 / cm_t[0];
      // Get element ivec of current eigenvector
      double evectorElement = *(modinfo_->Eigenvector(mode) + ivec);
      double dS2 = modinfo_->Eigenvalue(mode) * evectorElement * evectorElement;
      // Calculate weighted contribution of this mode to Cj(t) for this vector.
      for (int k = 0; k < nsteps; k++)
      {
        cj_t[k] += (dS2 * cm_t[k] * Norm1);
      }
    }
    // DEBUG: Write cj_t
    CpptrajFile cjt_out;
    cjt_out.OpenWrite( "dbg_cjt." + integerToString(ivec) + ".dat" );
    cjt_out.Printf("%-12s Cj(t)_%i\n", "#Time", ivec);
    for (int k = 0; k < nsteps; ++k)
      cjt_out.Printf("%12.4f %g\n", (double)k * tstep_, cj_t[k]);
  }

  // Define constants used in NMR relaxation calc.
  // Gyromagnetic ratio for H in rad s^-1 T^-1
  const double gamma_h = 2.6751987 * 1.0E8;
  // Gyromagnetic ratio for N in rad s^-1 T^-1 
  const double gamma_n = -2.7126 * 1.0E7;
  // Calculate magnetic field B in Teslas. Convert frequency to Hz (s^-1)
  const double Bzero = (Constants::TWOPI * freq_ * 1.0E6) / gamma_h;
  // Calculate Larmor freq. for H in rad s^-1
  const double omega_h = -gamma_h * Bzero;
  // Calculate Larmor freq. for N in rad s^-1
  const double omega_n = -gamma_n * Bzero;
  // Permeability of vacuum in T^2 m^3 J^-1
  const double mu_zero = Constants::FOURPI * 1.0E-7;
  // Planck's constant in J * s
  const double ha = 6.626176 * 1.0E-34;
  // Chemical shielding anisotropy constant in Hz (s^-1). 
  const double csa = -170.0 * 1.0E-6;
  // Convert distance in Angstroms to meters
  const double rnh = distnh_ * 1.0E-10; 
  // Calculate useful prefactors
  double FAC = (mu_zero * ha * gamma_n * gamma_h) /
               (8.0 * Constants::PI * Constants::PI * (rnh*rnh*rnh));
  FAC = (FAC * FAC) / 20.0; // Square it.
  double on2c2 = omega_n*omega_n * csa*csa;
  mprintf("DEBUG: omega_h= %g\nDEBUG: omega_n= %g\nDEBUG: c2= %g\nDEBUG: d2= %g\n",
          omega_h, omega_n, on2c2, FAC);
  mprintf("DEBUG: Jw(0, omega_h-omega_n)= %g\n", Jw(0, omega_h - omega_n, TauM));
  // Calculate Spectral Densities and NMR relaxation parameters
  CpptrajFile noeout;
  noeout.OpenWrite("dbg_noe.dat");
  noeout.Printf("# T1 T2 NOE\n");
  for (unsigned int ivec = 0; ivec != IredVectors_.size(); ivec++)
  {
    double JomegaN = Jw(ivec, omega_n, TauM);
    double JN3 = 3.0 * JomegaN;
    double JHminusN = Jw(ivec, omega_h - omega_n, TauM);
    double JHplusN6 = 6.0 * Jw(ivec, omega_h + omega_n, TauM);
    // Calculate T1
    double R1 = (FAC * (JN3 + JHminusN + JHplusN6))
                + (1/15.0) * on2c2 * JomegaN;
    // Calculate T2
    double Jzero4 = 4.0 * Jw(ivec, 0.0, TauM);
    double JH6 = 6.0 * Jw(ivec, omega_h, TauM);
    double R2 = ((FAC*0.5) * ( Jzero4 + JN3 + JHminusN + JH6 + JHplusN6))
                + ((1/90.0) * on2c2 * (Jzero4 + JN3));
    // Calculate NOE
    double cross_relax = FAC * (JHplusN6 - JHminusN);
    double NOE = 1.0 + ((gamma_h / gamma_n) * (1.0 / R1) * cross_relax);
    noeout.Printf("%u %g %g %g\n", ivec, 1.0/R1, 1.0/R2, NOE);
  }
  noeout.CloseFile();

/*
    // Relaxation calculation. Added by Alrun N. Koller & H. Gohlke
    CpptrajFile noefile;
    if (noefile.OpenWrite("dbg_noe.dat") != 0) {
      mprinterr("Error: Could not open NOE file for write.\n");
      return Analysis::ERR;
    }
    noefile.Printf("#vec     %10s   %10s   %10s\n","T1","T2","NOE");
    // conversion from Angstrom to meter
    double rnh = distnh_ * 1.0E-10;
    // ---------- CONSTANTS ----------
    // in H m^-1; permeability
    const double mu_zero = Constants::FOURPI * 1.0E-7;
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
    double b_zero = Constants::TWOPI * freq_ * 1.0E6 / gamma_h;
    // Next two both in rad s^-1
    const double lamfreqh = -1 * gamma_h * b_zero;
    const double lamfreqn = -1 * gamma_n * b_zero;
    double c2 = lamfreqn*lamfreqn * csa*csa;
    double d2 = (mu_zero * ha * gamma_n * gamma_h)/( 8.0 * (Constants::PI*Constants::PI) *
                                                     (rnh*rnh*rnh) );
    d2 = d2*d2; // fix from Junchao
    mprintf("DEBUG: omega_h= %g\nDEBUG: omega_n= %g\nDEBUG: c2= %g\nDEBUG: d2= %g\n",
            lamfreqh, lamfreqn, c2, d2);
    // -------------------------------
    taum_ = new double[ IredVectors_.size() ];
    std::copy( TauM.begin(), TauM.end(), taum_ ); 
    mprintf("DEBUG: Jw(0, omega_h-omega_n)= %g\n", calc_spectral_density(0, lamfreqh - lamfreqn));
    // loop over all vector elements --> have only one element/vector here; nvectelem = nelem
    for (int i = 0; i < (int)IredVectors_.size(); ++i) {
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
            - calc_spectral_density( i, lamfreqh - lamfreqn ) );

      double Noe = 1.0 + ( gamma_h * 1.0/gamma_n ) * ( 1.0 / R1 ) * Tj;
      noefile.Printf("%i %g %g %g\n", i, 1.0/R1, 1.0/R2, Noe);
    }
    noefile.CloseFile();
*/
  return Analysis::OK;
}
