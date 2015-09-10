#include "Analysis_IRED.h"
#include "CpptrajStdio.h"
#include "Constants.h" // PI
#include "Corr.h"
#include "DataSet_double.h"
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
  masterDSL_(0),
  modinfo_(0)
{}

void Analysis_IRED::Help() {
  mprintf("\t[relax freq <MHz> [NHdist <distnh>]] [order <order>]\n"
          "\ttstep <tstep> tcorr <tcorr> out <filename> [norm] [drct]\n"
          "\tmodes <modesname> [name <output sets name>]\n"
          "  Perform isotropic reorientational Eigenmode dynamics analysis.\n");
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
  // Get tstep, tcorr, S2 and Cm(t)/Cj(t) filenames
  tstep_ = analyzeArgs.getKeyDouble("tstep", 1.0);
  tcorr_ = analyzeArgs.getKeyDouble("tcorr", 10000.0);
  DataFile* orderout = DFLin->AddDataFile(analyzeArgs.GetStringKey("orderparamfile"), analyzeArgs);
  DataFile* outfile = 0;
  std::string filename = analyzeArgs.GetStringKey("out");
  if (!filename.empty()) {
    outfile  = DFLin->AddDataFile(filename);
    cmtfile_ = DFLin->AddDataFile(filename + ".cmt");
    cjtfile_ = DFLin->AddDataFile(filename + ".cjt");
  }
  // Output data sets
  dsname_ = analyzeArgs.GetStringKey("name");
  if (dsname_.empty()) dsname_ = DSLin->GenerateDefaultName("IRED");
  data_s2_ = DSLin->AddSetAspect(DataSet::FLOAT, dsname_, "S2");
  if (data_s2_ == 0) return Analysis::ERR;
  data_s2_->SetPrecision(10,5);
  if (orderout != 0) orderout->AddSet( data_s2_ );
  data_plateau_ = DSLin->AddSetAspect(DataSet::DOUBLE, dsname_, "Plateau");
  if (data_plateau_ == 0) return Analysis::ERR;
  data_plateau_->SetPrecision(12,8);
  data_tauM_ = DSLin->AddSetAspect(DataSet::DOUBLE, dsname_, "TauM");
  if (data_tauM_ == 0) return Analysis::ERR;
  data_tauM_->SetPrecision(12,6);
  if (outfile != 0) {
    outfile->AddSet( data_plateau_ );
    outfile->AddSet( data_tauM_ );
  }
  // Get norm, drct, relax
  norm_ = analyzeArgs.hasKey("norm");
  drct_ = analyzeArgs.hasKey("drct");
  relax_ = analyzeArgs.hasKey("relax");
  // Relax parameters
  DataFile* noefile = 0;
  if (relax_) {
    noefile = DFLin->AddDataFile(analyzeArgs.GetStringKey("noefile"), analyzeArgs);
    data_t1_  = DSLin->AddSetAspect(DataSet::DOUBLE, dsname_, "T1");
    data_t2_  = DSLin->AddSetAspect(DataSet::DOUBLE, dsname_, "T2");
    data_noe_ = DSLin->AddSetAspect(DataSet::DOUBLE, dsname_, "NOE");
    if (data_t1_ == 0 || data_t2_ == 0 || data_noe_ == 0) return Analysis::ERR;
    data_t1_->SetPrecision(10,5);
    data_t2_->SetPrecision(10,5);
    data_noe_->SetPrecision(10,5);
    if (noefile != 0) {
      noefile->AddSet( data_t1_ );
      noefile->AddSet( data_t2_ );
      noefile->AddSet( data_noe_ );
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
  mprintf("    IRED: %u IRED vectors.\n", IredVectors_.size());
  mprintf("\tData set name: %s\n", dsname_.c_str());
  if (orderout != 0)
    mprintf("\tOrder parameters will be written to %s\n", orderout->DataFilename().full());
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
  if (cmtfile_ != 0)
    mprintf("\tCm(t) functions will be written to %s\n", cmtfile_->DataFilename().full());
  if (outfile != 0)
    mprintf("\tCm(t->T) and TauM values will be written to %s\n", outfile->DataFilename().full());
  if (cjtfile_ != 0)
    mprintf("\tCj(t) functions will be written to %s\n", cjtfile_->DataFilename().full());
  mprintf("\tIRED modes will be taken from DataSet %s\n", modinfo_->legend());
  if (relax_) {
    mprintf("\t\tTauM, relaxation rates, and NOEs are calculated using the iRED\n"
            "\t\t  approach using an NH distance of %lf Ang. and a frequency of %lf MHz\n",
            distnh_, freq_);
    if (noefile != 0)
      mprintf("\t\tNOEs and relaxation rates will be written to %s\n",
              noefile->DataFilename().full());
  }
  mprintf("#Citation: Prompers, J. J.; Brüschweiler, R.; \"General framework for\n"
          "#          studying the dynamics of folded and nonfolded proteins by\n"
          "#          NMR relaxation spectroscopy and MD simulation\"\n"
          "#          J. Am. Chem. Soc. (2002) V.124 pp.4522-4534\n");
  masterDSL_ = DSLin;
  return Analysis::OK;
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
# ifdef TIMER
  Timer time_total, time_SH, time_cmt, time_tau, time_cjt, time_relax;
  time_total.Start();
# endif
  CorrF_FFT pubfft_;
  CorrF_Direct corfdir_;
  ComplexArray data1_;
  mprintf("\t'%s' has %zu modes.\n", modinfo_->legend(), modinfo_->Size());
  if ( modinfo_->Size() != IredVectors_.size() )
    mprintf("Warning: Number of IRED vectors (%zu) does not equal number of modes (%zu).\n",
            IredVectors_.size(), modinfo_->Size());
  // Calculation of S2 order parameters according to 
  //   Prompers & Brüschweiler, JACS  124, 4522, 2002; 
  // Loop over all vector elements
  mprintf("Info: Calculation of S2 parameters does not include first five modes.\n");
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
    float fval = (float)(1.0 - sum);
    data_s2_->Add(vi, &fval);
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
# ifdef TIMER
  time_SH.Start();
# endif
  // Ensure SH coords calculated for each ired vector.
  for (std::vector<DataSet_Vector*>::const_iterator iredvec = IredVectors_.begin();
                                                    iredvec != IredVectors_.end();
                                                  ++iredvec)
    (*iredvec)->CalcSphericalHarmonics( order_ );
# ifdef TIMER
  time_SH.Stop();
  time_cmt.Start();
# endif
  // Calculate Cm(t) for each mode
  DataSet_double& Plateau = static_cast<DataSet_double&>( *data_plateau_ );
  Plateau.Resize( modinfo_->Nmodes() ); // Sets all elements to 0.0
  CmtArray_.resize( modinfo_->Nmodes(), 0 );
  Dimension Tdim(0.0, tstep_, nsteps);
  for (int mode = 0; mode != modinfo_->Nmodes(); mode++)
  {
    // Add DataSet for Cm(t)
    CmtArray_[mode] = masterDSL_->AddSetIdxAspect(DataSet::DOUBLE, dsname_, mode, "Cm(t)");
    if (cmtfile_ != 0) cmtfile_->AddSet( CmtArray_[mode] );
    DataSet_double& cm_t = static_cast<DataSet_double&>( *CmtArray_[mode] );
    cm_t.SetPrecision(12,8);
    cm_t.SetDim(Dimension::X, Tdim);
    cm_t.Resize( nsteps ); // Sets all elements to 0.0
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
  }
# ifdef TIMER
  time_cmt.Stop();
  time_tau.Start();
#endif
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
    for (int i = 1; i < maxsteps; i++)
    {
      //double b_minus_a = ((double)i * tstep_) - ((double)(i-1) * tstep_);
      double curr_val = (cm_t[i  ] * (Norm1 / (double)(Nframes_ -  i     ))) - Cplateau;
      double prev_val = (cm_t[i-1] * (Norm1 / (double)(Nframes_ - (i - 1)))) - Cplateau;
      //mprintf("\tcm_t-T[%i]= %g  cm_t-T[%i]=%g\n", i-1, prev_val, i, curr_val);
      sum += (tstep_ * (prev_val + curr_val) * 0.5);
    }
    double cm0 = cm_t[0] * (Norm1 / Nframes_);
    mprintf("Mode %i : Cm(0)= %g  Cplateau= %g  Sum= %g\n", mode, cm0, Cplateau, sum);
    TauM[mode] = sum / (cm0 - Cplateau);
  }
# ifdef TIMER
  time_tau.Stop();
  time_cjt.Start();
# endif
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
  CjtArray_.resize( IredVectors_.size(), 0 );
  for (unsigned int ivec = 0; ivec != IredVectors_.size(); ivec++)
  {
    // Add DataSet for Cj(t)
    CjtArray_[ivec] = masterDSL_->AddSetIdxAspect(DataSet::DOUBLE, dsname_, ivec, "Cj(t)");
    if (cjtfile_ != 0) cjtfile_->AddSet( CjtArray_[ivec] );
    DataSet_double& cj_t = static_cast<DataSet_double&>( *CjtArray_[ivec] );
    cj_t.SetPrecision(12,8);
    cj_t.Resize( nsteps ); // Set all elements to 0.0
    cj_t.SetDim(Dimension::X, Tdim);
    // Calculate dS^2 for this vector and mode.
    for (int mode = 0; mode != modinfo_->Nmodes(); ++mode)
    {
      DataSet_double const& cm_t = static_cast<DataSet_double const&>( *CmtArray_[mode] );
      double Norm1 = Nframes_ / cm_t[0];
      // Get element ivec of current eigenvector
      double evectorElement = *(modinfo_->Eigenvector(mode) + ivec);
      double dS2 = modinfo_->Eigenvalue(mode) * evectorElement * evectorElement;
      // Calculate weighted contribution of this mode to Cj(t) for this vector.
      // Cm(t) must be normalized to 1.0.
      for (int k = 0; k < nsteps; k++)
      {
        cj_t[k] += (dS2 * cm_t[k] * (Norm1 / (Nframes_ - k)));
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
      double Norm1 = Nframes_ / cm_t[0];
      Plateau[mode] *= Norm1;
      for (int k = 0; k < nsteps; ++k)
        cm_t[k] *= (Norm1 / (Nframes_ - k));
    }
  } else {
    for (int mode = 0; mode != modinfo_->Nmodes(); mode++)
    {
      DataSet_double& cm_t = static_cast<DataSet_double&>( *CmtArray_[mode] );
      // 4*PI / ((2*order)+1) due to spherical harmonics addition theorem
      double Norm1 = DataSet_Vector::SphericalHarmonicsNorm( order_ );
      for (int k = 0; k < nsteps; ++k)
        cm_t[k] *= (Norm1 / (Nframes_ - k));
    }
  }

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
    mprintf("DEBUG: omega_h= %g\nDEBUG: omega_n= %g\nDEBUG: c2= %g\nDEBUG: d2= %g\n",
            omega_h, omega_n, on2c2, FAC);
    mprintf("DEBUG: Jw(0, omega_h-omega_n)= %g\n", Jw(0, omega_h - omega_n, TauM_s));
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
