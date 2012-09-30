#include <cstring> // memset
#include "Analysis_Timecorr.h"
#include "CpptrajStdio.h"
#include "Constants.h" // PI
// FIXME: Currently cross correlation FFT is inverted w.r.t. old ptraj was (corfdrct)
/// Strings corresponding to modes, used in output.
const char Analysis_Timecorr::ModeString[3][6] = {
  "????", "Auto", "Cross"
};

// CONSTRUCTOR
Analysis_Timecorr::Analysis_Timecorr() :
  mode_(M_UNKNOWN),
  type_(T_UNKNOWN),
  dplr_(false),
  norm_(false),
  drct_(false),
  relax_(false),
  tstep_(1),
  tcorr_(10000),
  distnh_(1.02),
  freq_(-1),
  ndata_(0),
  order_(2),
  ibeg_(0),
  iend_(0),
//  npair_(0),
  table_(0),
  data1_(0),
  data2_(0),
  cf_(0),
  cfinf_(0),
  p2cf_(0),
  rcf_(0),
  cf_cjt_(0),
  modinfo_(0),
  taum_(0),
  vinfo1_(0),
  vinfo2_(0)
{}

// DESTRUCTOR
Analysis_Timecorr::~Analysis_Timecorr() {
  if (cf_!=0) delete[] cf_;
  if (cfinf_!=0) delete[] cfinf_;
  if (p2cf_!=0) delete[] p2cf_;
  if (rcf_!=0) delete[] rcf_;
  if (data1_!=0) delete[] data1_;
  if (data2_!=0) delete[] data2_;
  if (table_!=0) delete[] table_;
  if (cf_cjt_!=0) delete[] cf_cjt_;
  if (modinfo_!=0) {
    if (modinfo_->Source() == ModesInfo::MS_FILE)
      delete modinfo_;
  }
  if (taum_!=0) delete[] taum_;
}

// Analysis_Timecorr::Setup()
/** analyze timecorr vec1 <vecname1> [vec2 <vecname2>]
  *                  [relax] [freq <hz>] [NHdist <distnh>] [order <order>]
  *                  tstep <tstep> tcorr <tcorr> out <filename>
  *                  [ corrired modes <modesname> [beg <ibeg> end <iend>] ]
  */
int Analysis_Timecorr::Setup(DataSetList* DSLin) {
  std::string modesfile;
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
  // Get order for Legendre polynomial
  order_ = analyzeArgs_.getKeyInt("order",2);
  if (order_ < 0 || order_ > 2) {
    mprintf("Warning: vector order out of bounds (<0 or >2), resetting to 2.\n");
    order_ = 2;
  }
  // Determine mode and type
  if (analyzeArgs_.hasKey("corrired")) {
    // TODO: Check that vinfo1 is IRED?
    // CORRIRED
    mode_ = AUTO;
    type_ = IRED;
    if (vinfo2_!=0) {
      mprintf("Warning: analyze timecorr: only calculating IRED corr for vec1, ignoring vec2\n");
      vinfo2_ = 0;
    }
    // Count and store the number of previously defined IRED vectors.
    DataSet_Vector* Vtmp;
    DSLin->VectorBegin();
    while ( (Vtmp = (DataSet_Vector*)DSLin->NextVector()) != 0 ) {
      if (Vtmp->IsIred()) 
        IredVectors_.push_back( Vtmp );
    }
    if (IredVectors_.empty()) {
      mprinterr("Error: analyze timecorr: corrired: No IRED vectors defined.\n");
      return 1;
    }
    // Get Pair number
/*    npair_ = analyzeArgs_.getKeyInt("npair",0);
    if (npair_ == 0) {
      mprinterr("Error: analyze timecorr: No 'npair <#>' arg given, needed for 'corrired'.\n");
      return 1;
    }
    // Actual pair number needs to start from 0
    --npair_;*/
    ibeg_ = analyzeArgs_.getKeyInt("beg",1);
    iend_ = analyzeArgs_.getKeyInt("end", 50);
    // Get modes name
    modesfile = analyzeArgs_.GetStringKey("modes");
    if (modesfile.empty()) {
      mprinterr("Error: analyze timecorr: No 'modes' <name> arg given, needed for 'corrired'.\n");
      return 1;
    }
    // Check if modes name exists on the stack
    modinfo_ = (ModesInfo*)DSLin->FindSetOfType( modesfile, DataSet::MODES );
    if (modinfo_ == 0) {
      // If not on stack, check for file.
      if ( fileExists(modesfile.c_str()) ) {
        modinfo_ = new ModesInfo;
        if (modinfo_->ReadEvecFile( modesfile, ibeg_, iend_ )) return 1;
      }
    }
    if (modinfo_ == 0) {
      mprinterr("Error: analyze timecorr: Modes %s DataSet/file not found.\n");
      return 1;
    }
  } else {
    // CORR, CORRPLANE
    type_ = NORMAL;
    if (vinfo2_ == 0)
      mode_ = AUTO;
    else 
      mode_ = CROSS;
  }
  // Sanity check
  if (mode_ == M_UNKNOWN || type_ == T_UNKNOWN) {
    mprinterr("Internal Error: analyze timecorr: Unknown mode or type.\n");
    return 1;
  }

  // Get tstep, tcorr, filename
  tstep_ = analyzeArgs_.getKeyDouble("tstep", 1.0);
  tcorr_ = analyzeArgs_.getKeyDouble("tcorr", 10000.0);
  noeFilename_ = analyzeArgs_.GetStringKey("noefile");
  filename_ = analyzeArgs_.GetStringKey("out");
  if (filename_.empty()) {
    mprinterr("Error: analyze timecorr: No outfile given ('out <filename>').\n");
    return 1;
  }

  // Get dplr, norm, drct, relax
  dplr_ = analyzeArgs_.hasKey("dplr");
  norm_ = analyzeArgs_.hasKey("norm");
  drct_ = analyzeArgs_.hasKey("drct");
  relax_ = analyzeArgs_.hasKey("relax");

  // Relax parameters
  if (relax_) {
    // Get freq, NH distance
    freq_ = analyzeArgs_.getKeyDouble("freq", -1.0);
    if (freq_ == -1.0) {
      mprinterr("Error: analyze timecorr: No frequency for calculation of relaxation\n");
      mprinterr("       parameters given ('freq <frequency>').\n");
      return 1;
    }
    // 1.02 * 10**(-10) in Angstroms
    distnh_ = analyzeArgs_.getKeyDouble("NHdist", 1.02);
  }

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
  if (type_ == IRED) {
    mprintf("\tIRED modes %i to %i are read from %s,\n", ibeg_, iend_,modesfile.c_str());
//    mprintf("\tand the pair %i is considered\n", npair_+1);
  }
  if (relax_) {
    mprintf("\t\tTauM, relaxation rates, and NOEs are calculated using the iRED\n");
    mprintf("\t\t  approach using an NH distance of %lf and a frequency of %lf\n",
            distnh_, freq_);
  }
  if (!noeFilename_.empty())
    mprintf("\t\tNOEs and relaxation rates will be written to %s\n",
            noeFilename_.c_str());
  mprintf("\t\tResults are written to %s\n", filename_.c_str());

  return 0;
}

// TODO: Make static functions below part of the class

/** Calculates correlation functions using the "direct" approach
  * (s. Comp. Sim. of Liquids, p.185)
  * - the result is not yet normalized by (no_of_discrete_data - t)**-1 (!)
  */
static void corfdir(int ndata, double *data1, double *data2, int nsteps, double *dtmp)
{
  int i, j;
  int ind1, ind2;
  double dsum, dsumi;
  int ndata2 = ndata / 2; // TODO: Pass in

  if (data2 == NULL) {
    for(i = 0; i < ndata2; i++){
      dsum = 0.0;
      for(j = i; j < ndata2; j++){
        ind1 = 2 * j;
        ind2 = 2 * (j-i);
        dsum += data1[ind1] * data1[ind2] + data1[ind1+1] * data1[ind2+1];
      }
      if(i < nsteps){
        ind1 = 2 * i;
        dtmp[ind1  ] = dsum;
        dtmp[ind1+1] = 0.0;
      }
      else{
        break;
      }
    }
  } else {
    for(i = 0; i < ndata2; i++){
      dsum = 0.0;
      dsumi = 0.0;
      for(j = i; j < ndata2; j++){
        ind1 = 2 * j;
        ind2 = 2 * (j-i);
        dsum  += data2[ind1] * data1[ind2  ] + data2[ind1+1] * data1[ind2+1];
        dsumi += data2[ind1] * data1[ind2+1] - data2[ind1+1] * data1[ind2  ];
      }
      if(i < nsteps){
        ind1 = 2 * i;
        dtmp[ind1  ] = dsum;
        dtmp[ind1+1] = dsumi;
      }
      else{
        break;
      }
    }
  }

  for(i = 0; i < nsteps; i++){
    ind1 = 2 * i;
    data1[ind1  ] = dtmp[ind1  ];
    data1[ind1+1] = dtmp[ind1+1];
  }
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
  //int ndata2 = 0;
  if (drct_) {
    ndata = 2 * frame;
    //ndata2 = frame;
    table_ = new double[ 2 * nsteps ];
  } else {
  //  ndata = ldexp( 1.0, (int)(log((double)(4 * frame - 1)) / log(2.0)) + 1);
    // Initialize FFT
    pubfft_.SetupFFT( frame );
    ndata = pubfft_.size() * 2;
    //ndata2 = pubfft_.size();
  }
  int mtot = 2 * order_ + 1;
  mprintf("CDBG: frame=%i time=%i nsteps=%i ndata=%i mtot=%i\n",frame,time,nsteps,ndata,mtot);

  // Allocate common memory and initialize arrays
  data1_ = new double[ ndata ];
  if (mode_ == CROSS)
    data2_ = new double[ ndata ];

  // -------------------- IRED CALCULATION ---------------------------
  if (type_ == IRED) {
    // Store ModesInfo
    int nvect = modinfo_->Nvect();
    int nvectelem = modinfo_->NvectElem();
    double* eigval = modinfo_->Freq();
    double* vout = modinfo_->Evec();
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
    int p2blocksize = 2 * mtot;             // Real + Img. for each -order <= m <= order
    int nsphereharm = frame * p2blocksize;  // Spherical Harmonics for each frame
    int ntotal = nvect * nsphereharm;       // Each vector has set of spherical harmonics
    double *cftmp1 = new double[ ntotal ];
    memset(cftmp1, 0, ntotal * sizeof(double));
    // Project spherical harmonics for each IRED vector on eigenmodes
    int n_ivec = 0;
    for (std::vector<DataSet_Vector*>::iterator ivec = IredVectors_.begin();
                                                ivec != IredVectors_.end(); ++ivec)
    {
      double *CF = cftmp1;
      double* sphereHarm = (*ivec)->SphericalHarmonics( order_ );
      // Loop over all eigenvectors
      for (int veci = 0; veci < nvect; ++veci) {
        double Qvec = modinfo_->Evec(veci, n_ivec);
        // Loop over all m = -L, ...., L
        for (int midx = -order_; midx <= order_; ++midx) {
          // Loop over spherical harmonic coords for this m (Complex, [Real][Img])
          for ( int sidx = 2 * (midx + order_); sidx < nsphereharm; sidx += p2blocksize) {
            *(CF++) += (Qvec * sphereHarm[sidx  ]);
            *(CF++) += (Qvec * sphereHarm[sidx+1]);
          }
        }
      }
      delete[] sphereHarm;
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
        for (int k = 0; k < frame*2; k += 2) {
          data1_[k  ]   = *CF;
          cfinfavgreal += *(CF++);
          data1_[k+1]   = *(CF);
          cfinfavgimg  += *(CF++);
          //mprintf("CDBG:\tVec=%i Frame=%i data1[%i]=%lf data1_[%i]=%lf\n",i,k,
          //        2*k, data1_[2*k], 2*k+1, data1_[2*k+1]);
        }
        cfinfavgreal /= frame;
        cfinfavgimg  /= frame;
        //mprintf("CDBG:\tVect[%i] cfinfavgreal=%lf cfinfavgimg=%lf\n",i,cfinfavgreal,cfinfavgimg);
        // Calc plateau value of correlation function (= C(m,t->T) in Bruschweiler paper (A20))
        cfinf_[veci] += (cfinfavgreal * cfinfavgreal) + (cfinfavgimg * cfinfavgimg);
        if (drct_) {
          // Calc correlation function (= C(m,l,t) in Bruschweiler paper) using direct approach
          corfdir(ndata, data1_, NULL, nsteps, table_);
        } else {
          // Pad with zero's at the end
          for (int k = 2*frame; k < ndata; ++k)
            data1_[k] = 0;
          // Calc correlation function (= C(m,l,t) in Bruschweiler paper) using FFT
          //corffft(ndata, data1_, NULL, table_);
          pubfft_.CorF_FFT(ndata, data1_, NULL);
        }
        // Sum into cf (= C(m,t) in Bruschweiler paper)
        for (int k = 0; k < nsteps; ++k) {
          //mprintf("CDBG: cf[%i] += data1[%i] (%lf)\n",idx1+k, 2*k, data1_[2 * k]);
          cf_[nsteps * veci + k] += data1_[2 * k];
        }
      }
    }
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
                 (cf_[nsteps * j + k] / (frame - k));
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
        mprinterr("Error: analyze timecorr: Different number of eigenmodes (%i) and\n", nvect);
        mprinterr("       eigenmode elements (%i)\n", nvectelem);
        return 1;
      }
      // consider only first third of the correlation curve to avoid fitting errors
      // NOTE: Done this way to be consistent with PTRAJ behavior. This really should
      //       be cast back to an integer.
      double new_nsteps = nsteps / 3.0;
      //mprintf("CDBG: new_nsteps= %lf\n",new_nsteps);
      for (int i = 0; i < nvect; ++i) {
        double integral = 0;
        double cfinfval = cfinf_[i] / cf_[nsteps * i] * frame;
        //mprintf("CDBG: cfinfval= %.10lE\n",cfinfval);
        for ( int j = 0 ; j < new_nsteps; ++j ) {
            double cfk = cf_[nsteps * i + j] * frame / (cf_[nsteps * i] * (frame - j));
            double cfk1 = cf_[nsteps * i + j + 1 ] * frame / (cf_[nsteps * i] * (frame - (j+1)));
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
        err = noefile.SetupWrite(NULL, debug_);
      else
        err = noefile.SetupWrite(noeFilename_, debug_);
      if (err != 0) {
        mprinterr("Error: analyze timecorr: Could not open NOE file for write.\n");
        return 1;
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
              modinfo_->calc_spectral_density( taum_, i, lamfreqh - lamfreqn ) 
              + 3.0 * modinfo_->calc_spectral_density( taum_, i, lamfreqn ) 
              + 6.0 * modinfo_->calc_spectral_density( taum_, i, lamfreqh + lamfreqn)
            ) + c2 / 15.0 * modinfo_->calc_spectral_density( taum_, i, lamfreqn );
        // Eq. A2 in Prompers & Brüschweiler, JACS  124, 4522, 2002
        double R2 = d2 / 40.0 * ( 4.0 *
              modinfo_->calc_spectral_density(  taum_, i, 0.0 ) 
              + modinfo_->calc_spectral_density( taum_, i, lamfreqh - lamfreqn )
              + 3.0 * modinfo_->calc_spectral_density( taum_, i, lamfreqn )
              + 6.0 * modinfo_->calc_spectral_density( taum_, i, lamfreqh ) 
              + 6.0 * modinfo_->calc_spectral_density( taum_, i, lamfreqh + lamfreqn)
            ) + c2 / 90.0 * ( 4.0 *
              modinfo_->calc_spectral_density( taum_, i, 0.0 ) 
              + 3.0 * modinfo_->calc_spectral_density( taum_, i, lamfreqn) );
        // Eq. A3 in Prompers & Brüschweiler, JACS  124, 4522, 2002
        double Tj = d2 / 20.0 * ( 6.0 *
              modinfo_->calc_spectral_density( taum_, i, lamfreqh + lamfreqn ) 
              - modinfo_->calc_spectral_density( taum_, i, lamfreqh + lamfreqn ) );

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
    if (cmtfile.OpenWrite( cmtname )) return 1;
    std::string cjtname = filename_ + ".cjt";
    CpptrajFile cjtfile;
    if (cjtfile.OpenWrite( cjtname )) return 1;  
    // Print headers
    cmtfile.Printf("%s-correlation functions Cm(t) for each eigenmode m, IRED type according to eq. A18 Prompers & Brüschweiler, JACS  124, 4522, 2002\n", ModeString[mode_]);
    cjtfile.Printf("%s-correlation functions Cj(t) for each ired vector j, IRED type according to eq. A23 Prompers & Brüschweiler, JACS  124, 4522, 2002\n", ModeString[mode_]);
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
        cmtfile.Printf("%12.8f", cfinf_[i] / cf_[nsteps * i] * frame);
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
          cmtfile.Printf("%12.8f", cf_[nsteps*j + i] * frame / (cf_[nsteps * j] * (frame - i)));
          cjtfile.Printf("%12.8f", cf_cjt_[nsteps*j + i] / cf_cjt_[nsteps*j]);
        } else {
          // 4/5*PI due to spherical harmonics addition theorem
          cmtfile.Printf("%12.8f", FOURFIFTHSPI * cf_[nsteps*j + i] / (frame - i));
          cjtfile.Printf("%12.8f", FOURFIFTHSPI * cf_cjt_[nsteps*j + i]);
        }
      }
      cmtfile.Printf("\n");
      cjtfile.Printf("\n");
    }
    cmtfile.CloseFile();
    cjtfile.CloseFile();

  // -------------------- NORMAL CALCULATION -------------------------
  } else if (type_ == NORMAL) {
    return 1; // DEBUG
/*    // Initialize memory
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
    // Loop over all three correlation functions
    double* dpt1 = 0;
    double* dpt2 = 0;
    for (int i = 0; i < 3; ++i) {
      if (i==0 && dplr_) {
        dpt1 = vinfo1_->Cftmp();
        if (vinfo2_!=0)
          dpt2 = vinfo2_->Cftmp();
      } else if (i==1) {
        dpt1 = vinfo1_->P2cftmp();
        if (vinfo2_!=0)
          dpt2 = vinfo2_->P2cftmp();
      } else if (i==2 && dplr_) {
        dpt1 = vinfo1_->Rcftmp();
        if (vinfo2_!=0)
          dpt2 = vinfo2_->Rcftmp();
      } else
        continue;

      // Loop over all m=-L, ..., L
      for (int j = 0; j < mtot; ++j) {
        // Loop over all snapshots
        for (int k = 0; k < frame; ++k) {
          if (i < 2) {
            int idx1 = 2 * (mtot * k + j);
            int idx2 = 2 * k;
            data1_[idx2  ] = dpt1[idx1  ];
            data1_[idx2+1] = dpt1[idx1+1];
            if(vinfo2_ != 0){
              data2_[idx2  ] = dpt2[idx1  ];
              data2_[idx2+1] = dpt2[idx1+1];
            }
          } else if (i == 2 && j == 0) {
            int idx1 = 2 * k;
            data1_[idx1  ] = dpt1[idx1  ];
            data1_[idx1+1] = dpt1[idx1+1];
            if(vinfo2_ != 0){
              data2_[idx1  ] = dpt2[idx1  ];
              data2_[idx1+1] = dpt2[idx1+1];
            }
          }
        }

        if (drct_) {
          // Calc correlation function using direct approach
          if (vinfo2_ == 0)
            corfdir(ndata, data1_, NULL, nsteps, table_);
          else
            corfdir(ndata, data1_, data2_, nsteps, table_);
        } else {
          // Pad with zero's at the end
          for (int k = 2 * frame; k < ndata; ++k) {
            data1_[k] = 0;
            if (vinfo2_ != 0)
              data2_[k] = 0;
          }
          // Calc correlation function using FFT
          if(vinfo2_ == 0)
            //corffft(ndata, data1_, NULL, table_);
            pubfft_.CorF_FFT(ndata, data1_, NULL);
          else
            //corffft(ndata, data1_, data2_, table_);
            pubfft_.CorF_FFT(ndata, data1_, data2_);
        }
        // Sum into cf
        for (int k = 0; k < nsteps; ++k){
          if(i == 0 && dplr_)
            cf_[k] += data1_[2 * k];
          else if(i == 1)
            p2cf_[k] += data1_[2 * k];
          else if(i == 2 && j == 0 && dplr_)
            rcf_[k] += data1_[2 * k];
          else
            break;
        }
      } // END j loop over m 
    } // END i loop over correlation fns

    // ----- PRINT NORMAL -----
    CpptrajFile outfile;
    if (outfile.OpenWrite(filename_)) return 1;
    outfile.Printf("%s-correlation functions, %s type\n",ModeString[mode_],"normal");
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
*/
  }

  return 0;
}

