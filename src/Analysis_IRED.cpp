#include "Analysis_IRED.h"
#include "CpptrajStdio.h"
#include "Constants.h" // PI
#include "Corr.h"
#include "DataSet_double.h" // Access to Resize and [] op
#include "DataSet_MatrixDbl.h" // Access to AddElement
#ifdef TIMER
# include "Timer.h"
#endif

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
  cmtfile_(0),
  cjtfile_(0),
  data_s2_(0),
  data_plateau_(0),
  data_tauM_(0),
  data_noe_(0),
  data_t1_(0),
  data_t2_(0),
  data_ds2_mat_(0),
  masterDSL_(0),
  modinfo_(0)
{}

void Analysis_IRED::Help() const {
  mprintf("\t[relax freq <MHz> [NHdist <distnh>]] [order <order>]\n"
          "\ttstep <tstep> tcorr <tcorr> out <filename> [norm] [drct]\n"
          "\tmodes <modesname> [name <output sets name>] [ds2matrix <file>]\n"
          "  Perform isotropic reorientational Eigenmode dynamics analysis.\n");
}

// Analysis_IRED::Setup()
Analysis::RetType Analysis_IRED::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  debug_ = debugIn;
  // Count and store the number of previously defined IRED vectors.
  for ( DataSetList::const_iterator DS = setup.DSL().begin(); DS != setup.DSL().end(); ++DS) {
    if ( (*DS)->Type() == DataSet::VECTOR && (*DS)->Meta().ScalarType() == MetaData::IREDVEC)
      IredVectors_.push_back( (DataSet_Vector*)*DS );
  }
  if (IredVectors_.empty()) {
    mprinterr("Error: No iRED vectors defined.\n");
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
  modinfo_ = (DataSet_Modes*)setup.DSL().FindSetOfType( modesfile, DataSet::MODES );
  if (modinfo_ == 0) {
    mprinterr("Error: %s\n", DataSet_Modes::DeprecateFileMsg);
    return Analysis::ERR;
  }
  // Get tstep, tcorr, S2 and Cm(t)/Cj(t) filenames
  tstep_ = analyzeArgs.getKeyDouble("tstep", 1.0);
  tcorr_ = analyzeArgs.getKeyDouble("tcorr", 10000.0);
  DataFile* orderout = setup.DFL().AddDataFile(analyzeArgs.GetStringKey("orderparamfile"), analyzeArgs);
  DataFile* outfile = 0;
  std::string filename = analyzeArgs.GetStringKey("out");
  if (!filename.empty()) {
    outfile  = setup.DFL().AddDataFile(filename);
    cmtfile_ = setup.DFL().AddDataFile(filename + ".cmt");
    cjtfile_ = setup.DFL().AddDataFile(filename + ".cjt");
  }
  DataFile* ds2matfile = setup.DFL().AddDataFile(analyzeArgs.GetStringKey("ds2matrix"));
  // Output data sets
  dsname_ = analyzeArgs.GetStringKey("name");
  if (dsname_.empty()) dsname_ = setup.DSL().GenerateDefaultName("IRED");
  data_s2_ = setup.DSL().AddSet(DataSet::FLOAT, MetaData(dsname_, "S2"));
  if (data_s2_ == 0) return Analysis::ERR;
  data_s2_->SetupFormat().SetFormatWidthPrecision(10,5);
  if (orderout != 0) orderout->AddDataSet( data_s2_ );
  data_plateau_ = setup.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "Plateau"));
  if (data_plateau_ == 0) return Analysis::ERR;
  data_plateau_->SetupFormat().SetFormatWidthPrecision(12,8);
  data_tauM_ = setup.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "TauM"));
  if (data_tauM_ == 0) return Analysis::ERR;
  data_tauM_->SetupFormat().SetFormatWidthPrecision(12,6);
  if (outfile != 0) {
    outfile->AddDataSet( data_plateau_ );
    outfile->AddDataSet( data_tauM_ );
  }
  if (ds2matfile != 0) {
    data_ds2_mat_ = setup.DSL().AddSet(DataSet::MATRIX_DBL, MetaData(dsname_, "dS2"));
    if (data_ds2_mat_ == 0) return Analysis::ERR;
    data_ds2_mat_->SetupFormat().SetFormatWidthPrecision(10,5);
    ds2matfile->ProcessArgs("square2d");
    ds2matfile->AddDataSet( data_ds2_mat_ );
  }
  // Get norm, drct, relax
  norm_ = analyzeArgs.hasKey("norm");
  drct_ = analyzeArgs.hasKey("drct");
  relax_ = analyzeArgs.hasKey("relax");
  // Relax parameters
  DataFile* noefile = 0;
  if (relax_) {
    noefile = setup.DFL().AddDataFile(analyzeArgs.GetStringKey("noefile"), analyzeArgs);
    data_t1_  = setup.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "T1"));
    data_t2_  = setup.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "T2"));
    data_noe_ = setup.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "NOE"));
    if (data_t1_ == 0 || data_t2_ == 0 || data_noe_ == 0) return Analysis::ERR;
    data_t1_->SetupFormat().SetFormatWidthPrecision(10,5);
    data_t2_->SetupFormat().SetFormatWidthPrecision(10,5);
    data_noe_->SetupFormat().SetFormatWidthPrecision(10,5);
    if (noefile != 0) {
      noefile->AddDataSet( data_t1_ );
      noefile->AddDataSet( data_t2_ );
      noefile->AddDataSet( data_noe_ );
    }
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
  mprintf("    IRED: %u iRED vectors.\n", IredVectors_.size());
  mprintf("\tData set name: %s\n", dsname_.c_str());
  if (orderout != 0)
    mprintf("\tOrder parameters will be written to '%s'\n", orderout->DataFilename().full());
  mprintf("\tCorrelation time %g, time step %g\n", tcorr_, tstep_);
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
  if (cmtfile_ != 0)
    mprintf("\tCm(t) functions will be written to '%s'\n", cmtfile_->DataFilename().full());
  if (outfile != 0)
    mprintf("\tCm(t->T) and TauM values will be written to '%s'\n", outfile->DataFilename().full());
  if (cjtfile_ != 0)
    mprintf("\tCj(t) functions will be written to '%s'\n", cjtfile_->DataFilename().full());
  if (data_ds2_mat_ != 0)
     mprintf("\tFull delta*S^2 matrix (# iRED vec rows by # eigenmodes cols) will be calcd\n"
             "\t  and written to '%s'\n", ds2matfile->DataFilename().full());
  mprintf("\tiRED modes will be taken from DataSet '%s'\n", modinfo_->legend());
  if (relax_) {
    mprintf("\tRelaxation rates and NOEs will be calculated using the iRED\n"
            "\t  approach using an NH distance of %.2f Ang. and a frequency of %.2f MHz\n",
            distnh_, freq_);
    if (noefile != 0)
      mprintf("\tNOEs and relaxation rates will be written to '%s'\n",
              noefile->DataFilename().full());
  }
  mprintf("#Citation: Prompers, J. J.; Brüschweiler, R.; \"General framework for\n"
          "#          studying the dynamics of folded and nonfolded proteins by\n"
          "#          NMR relaxation spectroscopy and MD simulation\"\n"
          "#          J. Am. Chem. Soc. (2002) V.124 pp.4522-4534\n");
  masterDSL_ = setup.DSL_Ptr();
  return Analysis::OK;
}

// Calculate spectral density for given vector and omega.
double Analysis_IRED::Jw(int ivec, double omega, std::vector<double> TauM) const
{
  double Jval = 0.0;
  for (int mode = 0; mode != modinfo_->Nmodes(); ++mode)
  {
    // Get element ivec of current eigenvector
    double evectorElement = modinfo_->Eigenvector(mode)[ivec];
    Jval += ((modinfo_->Eigenvalue(mode) * evectorElement * evectorElement)) *
            (2.0 * TauM[mode]) / (1.0 + omega*omega * TauM[mode]*TauM[mode]);
  }
  return Jval;
}

// Analysis_IRED::Analyze()
Analysis::RetType Analysis_IRED::Analyze() {
# ifdef TIMER
  Timer time_total, time_SH, time_cmt, time_tau, time_cjt, time_relax;
  time_total.Start();
# endif
  mprintf("\t'%s' has %zu eigenmodes.\n", modinfo_->legend(), modinfo_->Size());
  if ( modinfo_->Size() != IredVectors_.size() ) {
    mprinterr("Error: Number of iRED vectors (%zu) does not equal number of eigenmodes (%zu).\n",
              IredVectors_.size(), modinfo_->Size());
    return Analysis::ERR;
  }
  if ( modinfo_->VectorSize() != (int)IredVectors_.size() ) {
    mprinterr("Error: Number of iRED vectors (%zu) does not equal eigenvector length (%i).\n",
              IredVectors_.size(), modinfo_->VectorSize());
    return Analysis::ERR;
  }
  // ----- Calculation of S2 order parameters ----
  mprintf("Info: Calculation of S2 parameters does not include first five modes.\n");
  int startMode = 5;
  if (data_ds2_mat_ != 0) {
    if (((DataSet_MatrixDbl*)data_ds2_mat_)->Allocate2D(modinfo_->Nmodes(), modinfo_->VectorSize()))
      return Analysis::ERR;
    startMode = 0;
  }
  // Loop over all vector elements
  for (int vi = 0; vi < modinfo_->VectorSize(); ++vi) {
    // Sum according to Eq. A22 in Prompers & Brüschweiler, JACS 124, 4522, 2002
    double sum = 0.0;
    // Loop over all eigenvectors to calculate dS^2. Sum over all except the first
    // five ones (i.e. sum over all internal modes only).
    for (int mode = startMode; mode < modinfo_->Nmodes(); ++mode) {
      double Qvec = modinfo_->Eigenvector(mode)[vi];
      double ds2 = modinfo_->Eigenvalue(mode) * Qvec * Qvec;
      if (data_ds2_mat_ != 0)
        ((DataSet_MatrixDbl*)data_ds2_mat_)->AddElement( ds2 );
      if (mode > 4) sum += ds2;
    }
    float fval = (float)(1.0 - sum);
    data_s2_->Add(vi, &fval);
  }

  // All IRED vectors must have the same size
  int Nframes = -1;
  for (std::vector<DataSet_Vector*>::const_iterator Vtmp = IredVectors_.begin();
                                                    Vtmp != IredVectors_.end(); ++Vtmp)
  { 
    if (Nframes == -1)
      Nframes = (*Vtmp)->Size();
    else if (Nframes != (int)(*Vtmp)->Size()) {
      mprinterr("Error: All iRED vectors must have the same size.\n"
                "Error:   Vector %s size = %i, first vector size = %i\n",
                (*Vtmp)->legend(), (*Vtmp)->Size(), Nframes);
      return Analysis::ERR;
    }
  }

  // Determine max length of correlation functions
  int time = (int)(tcorr_ / tstep_) + 1;
  // nsteps
  int nsteps = 0;
  if (time > Nframes)
    nsteps = Nframes;
  else
    nsteps = time;
  mprintf("\tCorrelation functions calculated from t=0 to %g\n", (double)nsteps * tstep_);
  // Allocate memory to hold complex numbers for direct or FFT
  CorrF_FFT pubfft_;
  CorrF_Direct corfdir_;
  ComplexArray data1_;
  if (drct_) {
    data1_.Allocate( Nframes );
    corfdir_.Allocate( nsteps );
  } else {
    // Initialize FFT
    pubfft_.Allocate( Nframes );
    data1_ = pubfft_.Array();
  }
  // ----- Cm(t) CALCULATION ---------------------
# ifdef TIMER
  time_cmt.Start();
# endif
  // Allocate temporary storage for projection of SH for every l on eigenvecs.
  // Each SH value has a real + imaginary component.
  // [m-2R0][m-2I0][m-2R1][m-2I1] ... [m-2RN][m-2IN][m-1R0][m-1I0] ... [m+2RN][m+2IN]
  int ltot = 2 * order_ + 1; // Total # l values
  std::vector<double> cf_tmp( modinfo_->Nmodes() * ltot * Nframes * 2, 0.0 );
  // Project SH for each IRED vector on eigenvectors
  for (unsigned int vidx = 0; vidx != IredVectors_.size(); vidx++)
  {
    std::vector<double>::iterator CF = cf_tmp.begin();
    // Ensure SH calced for order
#   ifdef TIMER
    time_SH.Start();
#   endif
    IredVectors_[vidx]->CalcSphericalHarmonics( order_ );
#   ifdef TIMER
    time_SH.Stop();
#   endif
    // Loop over all eigenvectors
    for (int mode = 0; mode != modinfo_->Nmodes(); mode++)
    {
      double Qvec = modinfo_->Eigenvector(mode)[vidx];
      // Loop over L = -order ... +order
      for (int Lval = -order_; Lval <= order_; ++Lval)
      {
        // Loop over SH coords for this l value (real, img)
        for (ComplexArray::iterator sh = IredVectors_[vidx]->SphericalHarmonics(Lval).begin();
                                    sh != IredVectors_[vidx]->SphericalHarmonics(Lval).end(); ++sh)
          *(CF++) += (Qvec * (*sh));
      }
    }
  }
  // Now sum each Cml(t) into total Cm(t) for each mode.
  DataSet_double& Plateau = static_cast<DataSet_double&>( *data_plateau_ );
  Plateau.Resize( modinfo_->Nmodes() ); // Sets all elements to 0.0
  CmtArray_.resize( modinfo_->Nmodes(), 0 );
  Dimension Tdim(0.0, tstep_);
  std::vector<double>::const_iterator CF = cf_tmp.begin();
  for (int mode = 0; mode != modinfo_->Nmodes(); mode++)
  {
    // Add DataSet for Cm(t)
    CmtArray_[mode] = masterDSL_->AddSet(DataSet::DOUBLE, MetaData(dsname_, "Cm(t)", mode));
    if (cmtfile_ != 0) cmtfile_->AddDataSet( CmtArray_[mode] );
    DataSet_double& cm_t = static_cast<DataSet_double&>( *CmtArray_[mode] );
    cm_t.SetupFormat().SetFormatWidthPrecision(12,8);
    cm_t.SetDim(Dimension::X, Tdim);
    cm_t.Resize( nsteps ); // Sets all elements to 0.0
    // Loop over L = -order ... +order
    for (int Lval = -order_; Lval <= order_; Lval++)
    {
      // Loop over all snapshots
      double plateau_r = 0;
      double plateau_i = 0;
      for (int k = 0; k < Nframes*2; k += 2) {
        data1_[k  ] = *CF;
        plateau_r  += *(CF++);
        data1_[k+1] = *CF;
        plateau_i  += *(CF++);
      }
      plateau_r /= (double)Nframes;
      plateau_i /= (double)Nframes;
      // Calc plateau value of correlation function, Cm(t->T) in Bruschweiler paper (A20)
      Plateau[mode] += (plateau_r * plateau_r) + (plateau_i * plateau_i);
      // Calc correlation function for this mode and l, Cml(t)
      if (drct_)
        corfdir_.AutoCorr( data1_ );
      else {
        // Pad with zeros after Nframes 
        data1_.PadWithZero( Nframes );
        pubfft_.AutoCorr(data1_);
      }
      // Sum into Cm(t)
      for (int k = 0; k < nsteps; ++k)
        cm_t[k] += data1_[2 * k];
    }
  }
  cf_tmp.clear();
# ifdef TIMER
  time_cmt.Stop();
  time_tau.Start();
#endif
  // ----- TauM CALCULATION ----------------------
  // Calculate tau_m for each mode.
  DataSet_double& TauM = static_cast<DataSet_double&>( *data_tauM_ );
  TauM.Resize( modinfo_->Nmodes() ); // Sets all elements to 0.0
  for (int mode = 0; mode != modinfo_->Nmodes(); mode++)
  {
    DataSet_double const& cm_t = static_cast<DataSet_double const&>( *CmtArray_[mode] );
    // Only integrate a third of the way in. 
    //int maxsteps = nsteps / 3;
    int maxsteps = nsteps;
    // Get previously calculated plateau value for this.
    double Cplateau = Plateau[mode];
    // Integrate Cm(t) - Cplateau. Cm(t) has "standard" SH normalization.
    double sum = 0.0;
    double Norm1 = DataSet_Vector::SphericalHarmonicsNorm( order_ );
    double cm0 = cm_t[0] * (Norm1 / Nframes);
    double prev_val = cm0 - Cplateau; 
    for (int i = 1; i < maxsteps; i++)
    {
      double curr_val = (cm_t[i] * (Norm1 / (double)(Nframes - i))) - Cplateau;
      //mprintf("\tcm_t-T[%i]= %g  cm_t-T[%i]=%g\n", i-1, prev_val, i, curr_val);
      sum += (tstep_ * (prev_val + curr_val) * 0.5);
      prev_val = curr_val;
    }
    if (debug_ > 0)
      mprintf("Mode %i : Cm(0)= %g  Cplateau= %g  Sum= %g\n", mode, cm0, Cplateau, sum);
    TauM[mode] = sum / (cm0 - Cplateau);
  }
# ifdef TIMER
  time_tau.Stop();
  time_cjt.Start();
# endif
  // ----- Cj(t) CALCULATION ---------------------
  // Calculate Cj(t) for each vector j as weighted sum over Cm(t) arrays.
  // Cj(t) = SUM(m)[ dSjm^2 * Cm(t) ]
  // dSjm^2 = EVALm * (EVECm[i])^2
  CjtArray_.resize( IredVectors_.size(), 0 );
  for (unsigned int ivec = 0; ivec != IredVectors_.size(); ivec++)
  {
    // Add DataSet for Cj(t)
    CjtArray_[ivec] = masterDSL_->AddSet(DataSet::DOUBLE, MetaData(dsname_, "Cj(t)", ivec));
    if (cjtfile_ != 0) cjtfile_->AddDataSet( CjtArray_[ivec] );
    DataSet_double& cj_t = static_cast<DataSet_double&>( *CjtArray_[ivec] );
    cj_t.SetupFormat().SetFormatWidthPrecision(12,8);
    cj_t.Resize( nsteps ); // Set all elements to 0.0
    cj_t.SetDim(Dimension::X, Tdim);
    // Calculate dS^2 for this vector and mode.
    for (int mode = 0; mode != modinfo_->Nmodes(); ++mode)
    {
      DataSet_double const& cm_t = static_cast<DataSet_double const&>( *CmtArray_[mode] );
      double Norm1 = Nframes / cm_t[0];
      // Get element ivec of current eigenvector
      double evectorElement = *(modinfo_->Eigenvector(mode) + ivec);
      double dS2 = modinfo_->Eigenvalue(mode) * evectorElement * evectorElement;
      // Calculate weighted contribution of this mode to Cj(t) for this vector.
      // Cm(t) must be normalized to 1.0.
      for (int k = 0; k < nsteps; k++)
      {
        cj_t[k] += (dS2 * cm_t[k] * (Norm1 / (Nframes - k)));
      }
    }
  }
# ifdef TIMER
  time_cjt.Stop();
# endif
  // Normalize Cm(t) for output (also plateau values if necessary)
  if (norm_) {
    for (int mode = 0; mode != modinfo_->Nmodes(); mode++)
    {
      DataSet_double& cm_t = static_cast<DataSet_double&>( *CmtArray_[mode] );
      double Norm1 = Nframes / cm_t[0];
      Plateau[mode] *= Norm1;
      for (int k = 0; k < nsteps; ++k)
        cm_t[k] *= (Norm1 / (Nframes - k));
    }
  } else {
    for (int mode = 0; mode != modinfo_->Nmodes(); mode++)
    {
      DataSet_double& cm_t = static_cast<DataSet_double&>( *CmtArray_[mode] );
      // 4*PI / ((2*order)+1) due to spherical harmonics addition theorem
      double Norm1 = DataSet_Vector::SphericalHarmonicsNorm( order_ );
      for (int k = 0; k < nsteps; ++k)
        cm_t[k] *= (Norm1 / (Nframes - k));
    }
  }
  // ----- T1/T2/NOE CALCULATION -----------------
  if (relax_) {
#   ifdef TIMER
    time_relax.Start();
#   endif
    // Convert Tau from ps to s
    std::vector<double> TauM_s;
    for (unsigned int i = 0; i != TauM.Size(); ++i)
      TauM_s.push_back( TauM[i] * 1.0E-12 );
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
    if (debug_ > 0) {
      mprintf("DEBUG: omega_h= %g\nDEBUG: omega_n= %g\nDEBUG: c2= %g\nDEBUG: d2= %g\n",
              omega_h, omega_n, on2c2, FAC);
      mprintf("DEBUG: Jw(0, omega_h-omega_n)= %g\n", Jw(0, omega_h - omega_n, TauM_s));
    }
    // Calculate Spectral Densities and NMR relaxation parameters
    for (unsigned int ivec = 0; ivec != IredVectors_.size(); ivec++)
    {
      double JomegaN = Jw(ivec, omega_n, TauM_s);
      double JN3 = 3.0 * JomegaN;
      double JHminusN = Jw(ivec, omega_h - omega_n, TauM_s);
      double JHplusN6 = 6.0 * Jw(ivec, omega_h + omega_n, TauM_s);
      // Calculate T1
      double R1 = (FAC * (JN3 + JHminusN + JHplusN6))
                  + (1/15.0) * on2c2 * JomegaN;
      double T1 = 1.0 / R1;
      // Calculate T2
      double Jzero4 = 4.0 * Jw(ivec, 0.0, TauM_s);
      double JH6 = 6.0 * Jw(ivec, omega_h, TauM_s);
      double R2 = ((FAC*0.5) * ( Jzero4 + JN3 + JHminusN + JH6 + JHplusN6))
                  + ((1/90.0) * on2c2 * (Jzero4 + JN3));
      double T2 = 1.0 / R2;
      // Calculate NOE
      double cross_relax = FAC * (JHplusN6 - JHminusN);
      double NOE = 1.0 + ((gamma_h / gamma_n) * (1.0 / R1) * cross_relax);
      data_t1_->Add( ivec, &T1 );
      data_t2_->Add( ivec, &T2 );
      data_noe_->Add( ivec, &NOE );
    }
#   ifdef TIMER
    time_relax.Stop();
#   endif
  }
# ifdef TIMER
  time_total.Stop();
  time_SH.WriteTiming(2, "Spherical Harmonics", time_total.Total());
  time_cmt.WriteTiming(2, "Cm(t)", time_total.Total());
  time_tau.WriteTiming(2, "TauM", time_total.Total());
  time_cjt.WriteTiming(2, "Cj(t)", time_total.Total());
  time_relax.WriteTiming(2, "Relax", time_total.Total());
  time_total.WriteTiming(1, "Total IRED");
# endif
  return Analysis::OK;
}
