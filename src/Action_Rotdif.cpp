// Action_Rotdif
#include <cmath>
#include <cfloat> // DBL_MAX
#include <cstdio> //sscanf
#include <cstring> // memset
#include "Action_Rotdif.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // NumberFilename
#include "Constants.h" // TWOPI
#include "Integrate.h"
#include "ProgressBar.h"
#include "ComplexArray.h"

#ifndef NO_MATHLIB
// Definition of Fortran subroutines called from this class
extern "C" {
  // LAPACK
  void dgesvd_(char*, char*, int&, int&, double*,
               int&, double*, double*, int&, double*, int&,
               double*, int&, int& );
  void dsyev_(char*, char*, int&, double*, int&, double*,double*,int&,int&);
  // DEBUG
  //void dgemm_(char*,char*,int&,int&,int&,double&,
  //            double*,int&,double*,int&,double&,double*,int&);
  // DEBUG
  //void dgemv_(char*,int&,int&,double&,double*,int&,double*,int&,double&,double*,int&);
}
// ---------- Class: Vec6 ------------------------------------------------------
// Q_to_D()
/** Convert column vector Q with format:
  *   {Qxx, Qyy, Qzz, Qxy, Qyz, Qxz}
  * to diffusion tensor D via the formula:
  *   D = 3*Diso*I - 2Q
  * where
  *   Diso = trace(Q) / 3
  *
  * J. Biomol. NMR 9, 287 (1997) L. K. Lee, M. Rance, W. J. Chazin, A. G. Palmer
  * it has been assumed that the Q(l=1) has the same relationship to
  * D(l=1) as in the l=2 case; we have not yet proved this for the
  * non-symmetric case, but it is true for the axially symmetric case
  * also, Deff(l=1) = 1/(2*tau(l=1)) = e^T * Q(l=1) * e yields the 
  * correct values for tau(l=1) (known from Woessner model type 
  * correlation functions) along principal axe
  */
void Action_Rotdif::Vec6::Q_to_D(Matrix_3x3& D) const {
  double tq = Q_[0] + Q_[1] + Q_[2];
  D[0] = tq - (2 * Q_[0]); // tq-2Qxx
  D[1] = -2 * Q_[3];       // -2Qxy
  D[2] = -2 * Q_[5];       // -2Qxz
  D[3] = D[1];            // -2Qyx
  D[4] = tq - (2 * Q_[1]); // tq-2Qyy
  D[5] = -2 * Q_[4];       // -2Qyz
  D[6] = D[2];            // -2Qzx
  D[7] = D[5];            // -2Qzy
  D[8] = tq - (2 * Q_[2]); // tq-2Qzz
}

// D_to_Q()
/** Given diffusion tensor D[9], calculate Q[6] where Q is:
  *   Q = {Qxx, Qyy, Qzz, Qxy, Qyz, Qxz}
  * from
  *   Q = (3*Dav*I - D) / 2
  * where Dav = trace(D). See Q_to_D for more discussion.
  */
void Action_Rotdif::Vec6::D_to_Q(Matrix_3x3 const& D) {
  double td = D[0] + D[4] + D[8];
  Q_[0] = (td - D[0]) / 2; // Qxx
  Q_[1] = (td - D[4]) / 2; // Qyy
  Q_[2] = (td - D[8]) / 2; // Qzz
  Q_[3] = -D[1] / 2; // Qxy
  Q_[4] = -D[5] / 2; // Qyz
  Q_[5] = -D[2] / 2; // Qxz
}
#endif
// -----------------------------------------------------------------------------
// CONSTRUCTOR
Action_Rotdif::Action_Rotdif() :
  debug_(0),
  rseed_( 1 ),
  nvecs_( 0 ),
  tfac_( 0.0 ),
  ti_( 0.0 ),
  tf_( 0.0 ),
  NmeshPoints_( -1 ),
  itmax_( 0 ),
  delmin_( 0.0 ),
  d0_( 0.0 ),
  olegendre_( 2 ),
  ncorr_( 0 ),
  delqfrac_( 0 ),
  amoeba_ftol_(0.0000001 ),
  amoeba_itmax_(10000 ),
  do_gridsearch_( false ),
  useMass_(false),
  usefft_( true ),
  work_( 0 ),
  lwork_( 0 ),
  Tau_( 0 )
{ } 
// TODO: MAKE ANALYSIS
void Action_Rotdif::Help() {
  mprintf("\t[rseed <rseed>] [nvecs <nvecs>]\n");
  mprintf("\tref <refname> | refindex <refindex> | reference\n");
  mprintf("\t[<refmask>] [ncorr <ncorr>] dt <tfac> [ti <ti>] tf <tf>\n");
  mprintf("\t[itmax <itmax>] [tol <delmin>] [d0 <d0>] [order <olegendre>]\n");
  mprintf("\t[delqfrac <delqfrac>] [rvecout <randvecOut>]\n");
  mprintf("\t[rmout <rmOut>] [deffout <deffOut>] [outfile <outfilename>]\n");
  mprintf("\t[corrout <corrOut>] [usefft]\n");
  mprintf("\t[rvecin <randvecIn>]\n");
  mprintf("\t[gridsearch] [nmesh <NmeshPoints>]\n");
  mprintf("\tCalculate rotational diffusion tensor.\n");
}

// DESTRUCTOR
Action_Rotdif::~Action_Rotdif() {
  //fprintf(stderr,"Rotdif Destructor.\n");
  if (work_!=0) delete[] work_;
  // Close output file
  outfile_.CloseFile();
}

// Action_Rotdif::Init()
Action::RetType Action_Rotdif::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  debug_ = debugIn;
  // Get Keywords
  usefft_ = actionArgs.hasKey("usefft");
  nvecs_ = actionArgs.getKeyInt("nvecs",1000);
  rseed_ = actionArgs.getKeyInt("rseed",80531);
  ncorr_ = actionArgs.getKeyInt("ncorr",0);
  tfac_ = actionArgs.getKeyDouble("dt",0);
  if (tfac_<=0) {
    mprinterr("Error: Rotdif: 'dt <timestep>' must be specified and > 0.\n");
    return Action::ERR;
  }
  ti_ = actionArgs.getKeyDouble("ti",0);
  tf_ = actionArgs.getKeyDouble("tf",0);
  if (tf_ <= ti_) {
    mprinterr("Error: Rotdif: Initial time ti (%f) must be < final time tf (%f).\n",
              ti_, tf_);
    return Action::ERR;
  }
  NmeshPoints_ = actionArgs.getKeyInt("nmesh", -1);
  itmax_ = actionArgs.getKeyInt("itmax",500);
  delmin_ = actionArgs.getKeyDouble("tol",0.000001);
  d0_ = actionArgs.getKeyDouble("d0",0.03);
  olegendre_ = actionArgs.getKeyInt("order",2);
  if (olegendre_!=1 && olegendre_!=2) {
    mprinterr("Error: Rotdif: Order of legendre polynomial (%i) must be 1 or 2.\n",
              olegendre_);
    return Action::ERR;
  }
  if (usefft_ && olegendre_ != 2) {
    mprinterr("Error: Rotdif: FFT only supports legendre order 2.\n");
    return Action::ERR;
  }
  delqfrac_ = actionArgs.getKeyDouble("delqfrac",0.5);
  randvecOut_ = actionArgs.GetStringKey("rvecout");
  randvecIn_ = actionArgs.GetStringKey("rvecin");
  rmOut_ = actionArgs.GetStringKey("rmout");
  deffOut_ = actionArgs.GetStringKey("deffout");
  std::string outfilename = actionArgs.GetStringKey("outfile");
  corrOut_ = actionArgs.GetStringKey("corrout");
  do_gridsearch_ = actionArgs.hasKey("gridsearch");
  // Reference Keywords
  ReferenceFrame REF = FL->GetFrameFromArgs( actionArgs );
  if (REF.error()) return Action::ERR;
  // Get Masks
  AtomMask RefMask( actionArgs.GetMaskNext() );
  TargetMask_.SetMaskString( RefMask.MaskExpression() );
  
  // Initialize random number generator
  RNgen_.rn_set( rseed_ );

  // Set up reference for RMSD
  // Setup reference mask
  if (REF.Parm()->SetupIntegerMask( RefMask )) return Action::ERR;
  if (RefMask.None()) {
    mprintf("Error: Rotdif::init: No atoms in reference mask.\n");
    return Action::ERR;
  }
  // Allocate frame for selected reference atoms
  SelectedRef_.SetupFrameFromMask(RefMask, REF.Parm()->Atoms());
  // Set reference frame coordinates
  SelectedRef_.SetCoordinates(*(REF.Coord()), RefMask);
  // Always fitting; Pre-center reference frame
  SelectedRef_.CenterOnOrigin(useMass_); 

  // Open output file. Defaults to stdout if no name specified
  if (outfile_.OpenWrite(outfilename)) {
    mprinterr("Error opening Rotdif output file %s.\n", outfilename.c_str());
    return Action::ERR;
  }

  mprintf("    ROTDIF: Random seed %i, # of random vectors to generate: %i\n",rseed_,nvecs_);
  mprintf("\tMax length to compute time correlation functions:");
  if (ncorr_ == 0)
    mprintf(" Total # of frames\n");
  else
    mprintf(" %i\n",ncorr_);
  mprintf("\tTimestep = %.4lf, T0 = %.4lf, TF = %.4lf\n",tfac_,ti_,tf_);
  if (usefft_)
    mprintf("\tVector correlation functions will be calculated using FFT\n");
  else
    mprintf("\tVector correlation functions will be calculated directly\n");
  if (NmeshPoints_ != -1)
    mprintf("\tNumber of mesh points for interpolation is %i\n", NmeshPoints_);
  mprintf("\tIterative solver: Max iterations = %i, tol = %lf, initial guess = %lf\n",
          itmax_, delmin_, d0_);
  mprintf("\tOrder of Legendre polynomial = %i\n",olegendre_);
  mprintf("\tSimplex scaling factor=%.4lf\n",delqfrac_);
  if (do_gridsearch_)
    mprintf("\tGrid search will be performed for Q with full anisotropy (time consuming)\n");
  if (!randvecIn_.empty())
    mprintf("\tRandom vectors will be read from %s\n",randvecIn_.c_str());
  if (!randvecOut_.empty())
    mprintf("\tRandom vectors will be written to %s\n",randvecOut_.c_str());
  if (!rmOut_.empty())
    mprintf("\tRotation matrices will be written out to %s\n",rmOut_.c_str());
  if (!deffOut_.empty())
    mprintf("\tDeff will be written out to %s\n",deffOut_.c_str());
  if (!corrOut_.empty())
    mprintf("\tTime correlation for l=1 and l=2 for vector 0 will be written to %s\n",
            corrOut_.c_str());
#ifdef NO_MATHLIB
  mprintf("------------------------------------------------------\n");
  mprintf("Warning: Cpptraj was compiled with -DNO_MATHLIB.\n");
  mprintf("         The final tensor fit cannot be performed.\n");
  mprintf("         Only Deffs will be calculated.\n");
  mprintf("------------------------------------------------------\n");
#else
  if (!outfilename.empty())
    mprintf("            Diffusion constants and tau will be written to %s\n",
            outfilename.c_str());
  else
    mprintf("            Diffusion constants and tau will be written to STDOUT.\n");
#endif
  return Action::OK;
}

// Action_Rotdif::setup()
/** Determine what atoms each mask pertains to for the current parm file.
  * Also determine whether imaging should be performed.
  */
Action::RetType Action_Rotdif::Setup(Topology* currentParm, Topology** parmAddress) {

  if ( currentParm->SetupIntegerMask( TargetMask_ ) ) return Action::ERR;
  if ( TargetMask_.None() ) {
    mprintf("    Error: Rotdif::setup: No atoms in mask.\n");
    return Action::ERR;
  }
  // Allocate space for selected atoms in the frame. This will also put the
  // correct masses in based on the mask.
  SelectedTgt_.SetupFrameFromMask(TargetMask_, currentParm->Atoms());
  // Check that num atoms in frame mask from this parm match ref parm mask
  if ( SelectedRef_.Natom() != TargetMask_.Nselected() ) {
    mprintf( "    Error: Number of atoms in RMS mask (%i) does not \n",TargetMask_.Nselected());
    mprintf( "           equal number of atoms in Ref mask (%i).\n",SelectedRef_.Natom());
    return Action::ERR;
  }
  
  // Print info for this parm
  mprintf("    ROTDIF: %i atoms selected for RMS fit.\n",TargetMask_.Nselected());
        
  return Action::OK;  
}

// Action_Rotdif::action()
/** Calculate and store the rotation matrix for frame to reference.
  */
Action::RetType Action_Rotdif::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  Matrix_3x3 U;
  Vec3 Trans; // Unused

  // Set selected frame atoms. Masses have already been set.
  SelectedTgt_.SetCoordinates(*currentFrame, TargetMask_);
  SelectedTgt_.RMSD_CenteredRef(SelectedRef_, U, Trans, useMass_);
  Rmatrices_.push_back( U );

  return Action::OK;
} 

// ---------- ROTATIONAL DIFFUSION CALC ROUTINES -------------------------------
// Action_Rotdif::RandomVectors()
/** If no input file is specified by randvecIn, generate nvecs vectors of length 
  * 1.0 centered at the coordinate origin that are randomly oriented. The x,
  * y, and z components of each vector are generated from a polar coordinate 
  * system in which phi is randomly chosen in the range 0 to 2*PI and theta
  * is randomly chosen in the range 0 to PI/2 (single hemisphere).
  *   R = 1.0
  *   x = R * sin(theta) * cos(phi)
  *   y = R * sin(theta) * sin(phi)
  *   z = R * cos(theta)
  * \return Array containing random vectors; array will be empty on error.
  */
// NOTE: Theta could also be generated in the same way as phi. Currently done
//       to be consistent with the original implementation in randvec.F90
std::vector<Vec3> Action_Rotdif::RandomVectors() {
  std::vector<Vec3> XYZ;
  XYZ.reserve( nvecs_ );
  // ----- Read nvecs vectors from a file
  if (!randvecIn_.empty()) {
    CpptrajFile vecIn;
    if (vecIn.OpenRead(randvecIn_)) {
      mprinterr("Error: Could not open random vectors input file %s",randvecIn_.c_str());
      return XYZ;
    }
    //vecIn.SetupBuffer(); // NOTE: Use CpptrajFile internal buffer (1024 chars)
    double xIn[3];
    for (int i = 0; i < nvecs_; i++) {
      const char* buffer = vecIn.NextLine();
      if ( buffer == 0 ) {
        mprinterr("Error: Could not read vector %i from file %s\n", i+1,randvecIn_.c_str());
        XYZ.clear();
        return XYZ;
      }
      sscanf(buffer,"%*i %lf %lf %lf", xIn, xIn+1, xIn+2);
      XYZ.push_back( xIn );
    }
    vecIn.CloseFile();
  // ----- Generate nvecs vectors
  } else {
    for (int i = 0; i < nvecs_; i++) {
      double phi = TWOPI * RNgen_.rn_gen();
      // XYZ[i+2] is cos(theta)
      double costheta = 1 - RNgen_.rn_gen();
      double theta = acos( costheta );
      double sintheta = sin( theta );
      XYZ.push_back( Vec3(sintheta*cos(phi), sintheta*sin(phi), costheta) );
    }
  }
  // Print vectors
  if (!randvecOut_.empty()) {
    CpptrajFile rvout;
    if (rvout.OpenWrite(randvecOut_)) {
      mprinterr("Error: Rotdif: Could not set up %s for writing.\n",randvecOut_.c_str());
    } else {
      int idx = 1;
      for (std::vector<Vec3>::iterator vec = XYZ.begin(); vec != XYZ.end(); ++vec)
        rvout.Printf("%6i  %15.8lf  %15.8lf  %15.8lf\n",
                     idx++, (*vec)[0], (*vec)[1], (*vec)[2]);
      rvout.CloseFile();
    }
  }

  return XYZ;
}

// Action_Rotdif::compute_corr()
/** Given a normalized vector that has been randomly rotated itotframes 
  * times, compute the time correlation functions of the vector.
  * \param rotated_vectors array of vector coords for each frame 
  * \param maxdat Maximum length to compute time correlation functions (units of 'frames')
  * \param p2 Will be set with values for correlation function, l=2
  * \param p1 Will be set with values for correlation function, l=1
  */
// TODO: Make rotated_vectors const&
int Action_Rotdif::compute_corr(DataSet_Vector& rotated_vectors, int maxdat,
                                std::vector<double>& p2, std::vector<double>& p1)
{
  // Initialize p1 and p2
  p1.assign(maxdat, 0.0);
  p2.assign(maxdat, 0.0);
  int itotframes = rotated_vectors.Size();
  // i loop:  each value of i is a value of delay (correlation function argument)
  for (int i = 0; i < maxdat; i++) {
    int jmax = itotframes - i;
    for (int j = 0; j < jmax; j++) {
      // Dot vector j with vector j+i 
      double dot = rotated_vectors[j] * rotated_vectors[j+i];
      //printf("DBG i=%6i j=%6i k=%i  %f\n",i,j,i+j,dot);
      p2[i] += ((1.5*dot*dot) - 0.5);
      p1[i] += dot;
    }
    double one_jmax = (double) jmax;
    one_jmax = 1 / one_jmax;
    p2[i] = p2[i] * one_jmax;
    p1[i] = p1[i] * one_jmax;
  }
 
  return 0; 
}

// Action_Rotdif::fft_compute_corr()
int Action_Rotdif::fft_compute_corr(DataSet_Vector& rotated_vectors, int nsteps, 
                                    std::vector<double>& p2, int order)
{
  int n_of_vecs = rotated_vectors.Size();
  // Zero output array
  p2.assign(nsteps, 0.0);

  // Calculate spherical harmonics for each vector.
  rotated_vectors.CalcSphericalHarmonics( order );

  // Calculate correlation fn
  CorrF_FFT pubfft( n_of_vecs );
  ComplexArray data1 = pubfft.Array();
  // Loop over m = -olegendre, ..., +olegendre
  for (int midx = -order; midx <= order; ++midx) { 
    rotated_vectors.FillData( data1, midx );
    // Pad with zeros at the end
    // TODO: Does this always have to be done or can it just be done once
    data1.PadWithZero( n_of_vecs );
    pubfft.AutoCorr(data1);
    // Sum into pX
    for (int i = 0; i < nsteps; ++i)
      p2[i] += data1[i*2];
  }
  // Normalize correlation fn
  // 4/3*PI and 4/5*PI due to spherical harmonics addition theorem
  double norm = 1.0;
  if (order==1)
    norm = FOURTHIRDSPI;
  else if (order==2)
    norm = FOURFIFTHSPI;
  for (int i = 0; i < nsteps; ++i) {
    p2[i] *= (norm / (n_of_vecs - i));
  }
     
  return 0;
}

// Action_Rotdif::calcEffectiveDiffusionConst()
/** computes effect diffusion constant for a vector using the integral over
  * its correlation function as input. Starting with definition:
  *
  *   6*D=integral[0,inf;C(t)] 
  *
  * C(t) has already been integrated from ti -> tf yielding F(ti,tf).
  * Iteratively solves the equation 
  *
  *   D(i+1)=[exp(6*D(i)*ti)-exp(6*D(i)*tf)]/[6*F(ti,tf)]
  *
  * (numerator obtained by integrating exp(6*D*t) from ti -> tf)
  *
  * Modified so that itsolv now solves
  *
  * F(ti,tf;C(t)]=integral[ti,tf;C(t)]
  *
  * D(i+1)={exp[l*(l+1)*D(i)*ti]-exp[l*(l+1)*D(i)*tf)]}/[l*(l+1)*F(ti,tf)]
  * /param f Integral of Cl(t) from ti to tf
  * /return Effective value of D
  */
double Action_Rotdif::calcEffectiveDiffusionConst(double f ) {
// Class variables used:
//   ti,tf: Integration limits.
//   itmax: Maximum number of iterations in subroutine itsolv.
//   delmin: convergence criterion used in subroutine itsolv;
//           maximum accepted fractional change in successive 
//           iterations
//   d0: initial guess for diffusion constant; accurate estimate not needed
//   olegendre: order of Legendre polynomial in the correlation function <P(l)>
//
// Solves the equation 6*D=[exp(-6*D*ti)-exp(-6*D*tf)]/F(ti,tf) iteratively, 
// by putting 6*D(i+1)=[exp(-6*D(i)*ti)-exp(-6*D(i)*tf)]/F(ti,tf)
// where F(ti,tf) is input (integral[dt*C(t)] from ti->tf).
  double l, d, del, fac, di; 
  int i;

  // Always use d0 as initial guess
  di = d0_;
  l = (double) olegendre_;
  fac = (l*(l+1));
  i=1;
  d = 0;
  del = DBL_MAX;
  while ( i<=itmax_ && del>delmin_) {
     d = ( exp(-fac*di*ti_) - exp(-fac*di*tf_) );
     d = d / (fac*f);
     del = (d-di)/di;
     if (del < 0) del = -del; // Abs value
     if (debug_>2)
       mprintf("ITSOLV: %6i  %15.8e  %15.8e  %15.8e\n", i,di,d,del);
     di = d;
     ++i;
  }
  if ( i>itmax_ && del>delmin_) {
     mprintf("\tWarning, itsolv did not converge: # iterations=%i, fractional change=%lf\n",
             i, del);
  } else {
    if (debug_>1) mprintf("\tITSOLV Converged: # iterations=%i\n",i);
  }

  return d; 
}

// PrintMatrix()
/*static void PrintMatrix(CpptrajFile &outfile, const char *Title, double *U, 
                        int mrows, int ncols) 
{
  outfile.Printf("    %s",Title);
  int usize = mrows * ncols;
  for (int i = 0; i < usize; i++) {
    if ( (i%ncols)==0 ) outfile.Printf("\n");
    outfile.Printf(" %10.5lf",U[i]);
  }
  outfile.Printf("\n");
}*/
void Action_Rotdif::PrintMatrix(CpptrajFile& outfile, const char* Title, Matrix_3x3 const& U)
{
  outfile.Printf("    %s\n",Title);
  outfile.Printf(" %10.5f %10.5f %10.5f\n %10.5f %10.5f %10.5f\n %10.5f %10.5f %10.5f\n",
                 U[0], U[1], U[2], U[3], U[4], U[5], U[6], U[7], U[8]);
}

void Action_Rotdif::PrintVector(CpptrajFile& outfile, const char* Title, Vec3 const& V)
{
  outfile.Printf("    %s\n",Title);
  outfile.Printf(" %10.5f %10.5f %10.5f\n", V[0], V[1], V[2]);
}

void Action_Rotdif::PrintVec6(CpptrajFile& outfile, const char* Title, Vec6 const& V)
{
  outfile.Printf("    %s\n",Title);
  outfile.Printf(" %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", 
                 V[0], V[1], V[2], V[3], V[4], V[5]);
}

// Action_Rotdif::calc_Asymmetric()
/** Computes tau(l=1) and tau(l=2) with full anisotropy given principal 
  * components of a rotational diffusion tensor and angular coordinates 
  * of a vector (in the PA frame, 
  * n*ez = cos(theta), n*ey = sin(theta)*sin(phi), n*ex = sin(theta)*cos(phi).
  * It is expected that the principal components Dxyz and principal axes
  * D_matrix are in ascending order, i.e. Dx <= Dy <= Dz. It is also expected
  * that the eigenvectors are orthonormal. This is the default behavior when 
  * the LAPACK routines are used. Also, since LAPACK routines are being used
  * from fortran it is expected the matrix of eigenvectors will be in column
  * major order, i.e.:
  *   x0 y0 z0
  *   x1 y1 z1
  *   x2 y2 z2
  * \param Dxyz Vector of size 3 containing principal components
  * \param matrix_D matrix of size 3x3 containing orthonormal principal vectors
  *        in COLUMN MAJOR order.
  */
int Action_Rotdif::calc_Asymmetric(Vec3 const& Dxyz, Matrix_3x3 const& matrix_D) 
{
  double lambda[8];
  double dx = Dxyz[0];
  double dy = Dxyz[1];
  double dz = Dxyz[2];
  double Dav, Dpr2, delta;

  // DEBUG
  //printMatrix("Dxyz",Dxyz,1,3);
  //printMatrix("Dvec",matrix_D,3,3);

  // tau = sum(m){amp(l,m)/lambda(l,m)}
  // See Korzhnev DM, Billeter M, Arseniev AS, Orekhov VY; 
  // Prog. Nuc. Mag. Res. Spec., 38, 197 (2001) for details, Table 3 in 
  // particular. Only weights need to be computed for each vector, decay 
  // constants can be computed once. First three decay constants correspond 
  // to l=1, m=-1,0,+1; next five correspond to l=2, m=-2,-1,0,+1,+2.

  // l=1, m=-1 term: lambda(1,-1) = Dy + Dz
  lambda[0] = dy + dz;
  // l=1, m=0 term:  lambda(1, 0) = Dx + Dy
  lambda[1] = dx + dy;
  // l=1, m=+1 term: lambda(1,+1) = Dx + Dz
  lambda[2] = dx + dz;
  // l=2, m=-2 term: lambda(2,-2) = Dx + Dy + 4*Dz
  lambda[3] = dx + dy + (4*dz);
  // l=2, m=-1 term: lambda(2,-1) = Dx + 4*Dy + Dz
  lambda[4] = dx + (4*dy) + dz;
  // l=2, m=0 term:  lambda(2, 0) = 6*[Dav - sqrt(Dav*Dav - Dpr*Dpr)]
  //                 Dav = (Dx + Dy + Dz)/3 
  //                 Dpr = sqrt[(Dx*Dy + Dy*Dz + Dx*Dz)/3 ]
  //                 Dpr2 = Dpr*Dpr = (Dx*Dy + Dy*Dz + Dx*Dz)/3
  Dav = (dx + dy + dz) / 3;
  Dpr2 = ((dx*dy) + (dy*dz) + (dx*dz)) / 3;
  if (Dpr2 < 0) {
    //mprinterr("Error: Rotdif::calc_Asymmetric: Cannot calculate Dpr (Dpr^2 < 0)\n");
    // NOTE: Original code just set Dpr to 0 and continued. 
    Dpr2 = 0;
  }
  delta = (Dav*Dav) - Dpr2;
  if (delta < 0) {
    mprinterr("Error: Rotdif::calc_Asymmetric: Cannot calculate lambda l=2, m=0\n");
    return 1;
  }
  double sqrt_Dav_Dpr = sqrt( delta );
  lambda[5] = 6 * (Dav - sqrt_Dav_Dpr);
  // l=2, m=+1 term: lambda(2,+1) = 4*Dx + Dy + Dz
  lambda[6] = (4*dx) + dy + dz;
  // l=2, m=+2 term: lambda(2,+2) = 6*[Dav + sqrt(Dav*Dav - Dpr*Dpr)]
  lambda[7] = 6 * (Dav + sqrt_Dav_Dpr);

  // Check all normalization factors for 0.0, which can lead to inf values
  // propogating throughout the calc. Set to SMALL instead; will get 
  // very large numbers but should still not overflow.
  for (int i = 0; i < 8; i++) 
    if (lambda[i] < SMALL) lambda[i] = SMALL;

  // Loop over all random vectors
  int nvec = 0; // index into tau1, tau2, sumc2
  for (std::vector<Vec3>::iterator randvec = random_vectors_.begin();
                                   randvec != random_vectors_.end(); ++randvec)
  {
    // Rotate vector i into D frame
    // This is an inverse rotation. Since matrix_D is in column major order
    // however, do a normal rotation instead.
    Vec3 rotated_vec = matrix_D * (*randvec);
    double dot1 = rotated_vec[0];
    double dot2 = rotated_vec[1];
    double dot3 = rotated_vec[2];
    //mprintf("DBG: dot1-3= %10.5lf%10.5lf%10.5lf\n",dot1,dot2,dot3);
    
    // assuming e(3)*n = cos(theta), e(1)*n = sin(theta)*cos(phi),
    // e(2)*n = sin(theta)*sin(phi), theta >= 0;
    // sin(theta) = sqrt(1 - cos(theta)^2) and theta = tan^-1[sin(theta)/cos(theta)];
    // phi = tan^-1[sin(phi)/cos(phi)] = tan^-1[(e(2)*n)/(e(1)*n)]
    double theta = atan2( sqrt( 1 - (dot3*dot3) ), dot3 );
    double phi = atan2( dot2, dot1 );

    // Pre-calculate sines and cosines
    double sintheta = sin(theta);
    double costheta = cos(theta);
    double sinphi = sin(phi);
    double cosphi = cos(phi);

    // ----- Compute correlation time for l=1
    // tau(l=1) = [sin(theta)^2*cos(phi)^2]/(Dy+Dz) +
    //            [sin(theta)^2*sin(phi)^2]/(Dx+Dz) +
    //            [cos(theta)^2           ]/(Dx+Dy)
    double sintheta2 = sintheta * sintheta;
    //     cth2=cth*cth
    //     sphi2=sphi*sphi
    //     cphi2=cphi*cphi
    tau1_[nvec] = ((sintheta2 * (cosphi*cosphi)) / lambda[0]) +
               ((sintheta2 * (sinphi*sinphi)) / lambda[2]) +
               ((costheta*costheta)           / lambda[1]);

    // ----- Compute correlation time for l=2
    // pre-calculate squares
    double dot1_2 = dot1*dot1;
    double dot2_2 = dot2*dot2;
    double dot3_2 = dot3*dot3;
    //     m=-2 term:
    //     lambda(2,-2) = Dx + Dy + 4*Dz
    //     amp(2,-2) = 0.75*sin(theta)^4*sin(2*phi)^2 = 3*(l^2)*(m^2)
    double m_m2 = 3*dot1_2*dot2_2;
    //     m=-1 term:
    //     lambda(2,-1) = Dx + 4*Dy + Dz
    //     amp(2,-1) = 0.75*sin(2*theta)^2*cos(phi)^2 = 3*(l^2)*(n^2)
    double m_m1 = 3*dot1_2*dot3_2;
    //     m=0 term:
    //     lambda(2,0) = 6*[Dav - sqrt(Dav*Dav - Dpr*Dpr)]
    //         Dav = (Dx + Dy + Dz)/3 
    //         Dpr = sqrt[(Dx*Dy + Dy*Dz + Dx*Dz)/3 ]
    //     amp(2,0) = (w/N)^2*0.25*[3*cos(theta)^2-1]^2 +
    //                (u/N)^2*0.75*sin(theta)^4*cos(2*phi)^2 -
    //                (u/delta)*[sqrt(3)/8]*[3*cos(theta)^2-1]*sin(theta)^2*cos(2*phi)
    //         u = sqrt(3)*(Dx - Dy)
    //         delta = 3*sqrt(Dav*Dav - Dpr*Dpr)
    //         w = 2*Dz - Dx - Dy + 2*delta
    //         N = 2*sqrt(delta*w)
    delta = 3 * sqrt_Dav_Dpr;
    double dot1_4 = dot1_2 * dot1_2;
    double dot2_4 = dot2_2 * dot2_2;
    double dot3_4 = dot3_2 * dot3_2;
    double da = 0.25 * ( 3*(dot1_4 + dot2_4 + dot3_4) - 1);
    double ea = 0;
    if ( delta > SMALL) { 
      double epsx = 3*(dx-Dav)/delta;
      double epsy = 3*(dy-Dav)/delta;
      double epsz = 3*(dz-Dav)/delta;
      double d2d3 = dot2*dot3;
      double d1d3 = dot1*dot3;
      double d1d2 = dot1*dot2;
      ea = epsx*(3*dot1_4 + 6*(d2d3*d2d3) - 1) + 
           epsy*(3*dot2_4 + 6*(d1d3*d1d3) - 1) + 
           epsz*(3*dot3_4 + 6*(d1d2*d1d2) - 1);
      ea /= 12;
    }
    double m_0 = da + ea;
    //     m=+1 term:
    //     lambda(2,+1) = 4*Dx + Dy + Dz
    //     amp(2,+1) = 0.75*sin(2*theta)^2*sin(phi)^2 
    double m_p1 = 3*dot2_2*dot3_2;
    //     m=+2 term:
    //     lambda(2,+2) = 6*[Dav + sqrt(Dav*Dav = Dpr*Dpr)]
    //     amp(2,+2) = (u/n)^2*0.25*[3*cos(theta)^2-1]^2 +
    //                 (w/n)^2*0.75*sin(theta)^4*cos(2*phi)^2 +
    //                 (u/delta)*[sqrt(3)/8]*[3*cos(theta)^2-1]*sin(theta)^2*cos(2*phi)
    double m_p2 = da - ea;
    // Sum decay constants and weights for l=2 m=-2...2
    tau2_[nvec] = (m_m2 / lambda[3]) + (m_m1 / lambda[4]) + (m_0 / lambda[5]) +
               (m_p1 / lambda[6]) + (m_p2 / lambda[7]);
    sumc2_[nvec] = m_m2 + m_m1 + m_0 + m_p1 + m_p2;
    ++nvec;
  }

  return 0;
}

// Action_Rotdif::chi_squared()
/** Given Q, calculate effective D values for random_vectors with full 
  * anisotropy and compare to the effective D values in D_eff. Sets
  * D_XYZ with principal components and D_tensor with principal vectors
  * in column major order.
  * \return Deviation of calculated D values from given D values.
  */
double Action_Rotdif::chi_squared(Action_Rotdif::Vec6 const& Qin) {
#ifdef NO_MATHLIB
  return -1;
#else
  int n_cols = 3;
  int info;
  double chisq;

  // Convert input Q to Diffusion tensor D
  Qin.Q_to_D(D_tensor_);
  // Diagonalize D; it is assumed workspace (work, lwork) set up prior 
  // to this call.
  // NOTE: Due to the fortran call, the eigenvectors are returned in COLUMN
  //       MAJOR order.
  dsyev_((char*)"Vectors",(char*)"Upper", n_cols, D_tensor_.Dptr(), n_cols, 
         D_XYZ_.Dptr(), work_, lwork_, info);
  // Check for convergence
  //if (info > 0) {
  //  mprinterr("The algorithm computing the eigenvalues/eigenvectors of D failed to converge.\n");
  //  return -1;
  //}
  calc_Asymmetric(D_XYZ_, D_tensor_);

  chisq = 0;
  for (int i = 0; i < nvecs_; i++) {
    double diff = D_eff_[i] - (*Tau_)[i];
    chisq += (diff * diff); 
  }

  return chisq;  
#endif
}

// calculate_D_properties()
/** Given the principal components of D, calculate the isotropic diffusion
  * constant (Dav = (Dx + Dy + Dz) / 3), anisotropy (Dan = 2Dz / (Dx + Dy)),
  * and rhombicity (Drh = (3/2)(Dy - Dx) / (Dz - 0.5(Dx + Dy)).
  */
static Vec3 calculate_D_properties(Vec3 const& Dxyz) {
  double Dx = Dxyz[0];
  double Dy = Dxyz[1];
  double Dz = Dxyz[2];
  return Vec3( (Dx + Dy + Dz) / 3,
               (2 * Dz) / (Dx + Dy),
               (1.5 * (Dy - Dx)) / (Dz - (0.5 * (Dx + Dy))) );
}
// ---------- SIMPLEX MINIMIZER ------------------------------------------------
#define SM_NP 6
#define SM_NP1 7
// Amotry()
double Action_Rotdif::Amotry(double xsmplx[SM_NP1][SM_NP], double *ysearch,
                      Vec6& psum, int ihi, double fac)
{
  double ytry, fac1, fac2;
  Vec6 ptry;

  fac1 = (1-fac)/SM_NP;
  fac2 = fac1 - fac;
  for (int j=0; j < SM_NP; j++) {
    ptry[j] = (psum[j] * fac1) - (xsmplx[ihi][j] * fac2);
    //mprintf("\t\tAmotry: %6i%10.5lf\n",j,ptry[j]);
  }
  ytry = chi_squared( ptry );
  if (ytry < ysearch[ihi]) {
    ysearch[ihi] = ytry;
    for (int j = 0; j < SM_NP; j++) {
      psum[j] = psum[j] - xsmplx[ihi][j] + ptry[j];
      xsmplx[ihi][j] = ptry[j];
      //mprintf("\t\tAmotryX: %6i%6i%10.5lf\n",ihi+1,j+1,xsmplx[ihi][j]);
    }
  }
  return ytry;
}

// Amoeba()
/** Main driver for the simplex method */
int Action_Rotdif::Amoeba(double xsmplx[SM_NP1][SM_NP], double *ysearch) {
  int iter;
  bool loop1;
  bool loop2;
  Vec6 psum;
  int ilo, ihi, inhi;
  double rtol, swap;

  iter = 0;
  loop1 = true;
  while (loop1) {
    //mprintf("Hit loop one %6i\n",iter);
    for (int n = 0; n < SM_NP; n++) {
      psum[n] = 0;
      for (int m = 0; m < SM_NP1; m++) { 
        //mprintf("Xsmplx %6i%6i%10.5lf\n",m,n,xsmplx[m][n]);
        //mprintf("Ysearch %6i%10.5lf\n",m,ysearch[m]);
        psum[n] += xsmplx[m][n];
      }
    }
    loop2 = true;
    while (loop2) {
      //mprintf("Hit loop two %6i\n",iter);
      //for (int n = 0; n < SM_NP; n++)
      //  mprintf("\tPsum %6i%10.5lf\n",n,psum[n]);
      ilo = 0;
      if (ysearch[0] > ysearch[1]) {
        ihi = 0;
        inhi = 1;
      } else {
        ihi = 1;
        inhi = 0;
      }
      for (int i = 0; i < SM_NP1; i++) {
        if (ysearch[i] <= ysearch[ilo]) ilo = i;
        if (ysearch[i] > ysearch[ihi]) {
          inhi = ihi;
          ihi = i;
        } else if ( ysearch[i] > ysearch[inhi]) {
          if (i != ihi) inhi = i;
        }
      }  
      //mprintf("Yihi Yilo = %10.5lf%10.5lf\n",ysearch[ihi],ysearch[ilo]);
      //mprintf("\tYihi Yilo = %6i%6i\n",ihi+1,ilo+1);
      double abs_yhi_ylo = ysearch[ihi] - ysearch[ilo];
      if (abs_yhi_ylo < 0) abs_yhi_ylo = -abs_yhi_ylo;
      double abs_yhi = ysearch[ihi];
      if (abs_yhi < 0) abs_yhi = -abs_yhi;
      double abs_ylo = ysearch[ilo];
      if (abs_ylo < 0) abs_ylo = -abs_ylo;
      //mprintf("Abs(yihi - yilo)=%lf, Abs(yihi)=%lf, Abs(yilo)=%lf\n",
      //        abs_yhi_ylo,abs_yhi,abs_ylo);
      rtol = 2.0 * (abs_yhi_ylo / (abs_yhi + abs_ylo));
      if (rtol < amoeba_ftol_) {
        swap = ysearch[0];
        ysearch[0] = ysearch[ilo];
        ysearch[ilo] = swap;
        for (int n = 0; n < SM_NP; n++) {
          swap = xsmplx[0][n];
          xsmplx[0][n] = xsmplx[ilo][n];
          xsmplx[ilo][n] = swap;
        }
        return iter;
      }
      //mprintf("\tIn amoeba, iter=%i, rtol=%15.6le\n",iter,rtol); 

      if (iter >= amoeba_itmax_) {
        mprintf("Max iterations (%i) exceeded in amoeba.\n",amoeba_itmax_);
        return iter;
      }
      iter += 2;

      double ytry = Amotry(xsmplx, ysearch, psum, ihi, -1.0);
      //mprintf("\tYtry %6i%10.5f\n",iter,ytry);
      if (ytry <= ysearch[ilo]) { 
        ytry = Amotry(xsmplx, ysearch, psum, ihi, 2.0);
        //mprintf("\tCase 1 %10.5lf\n",ytry);
      } else if (ytry >= ysearch[inhi]) {
        double ysave = ysearch[ihi];
        ytry = Amotry(xsmplx, ysearch, psum, ihi, 0.5);
        //mprintf("\tCase 2 %10.5lf\n",ytry);
        if (ytry >= ysave) {
          for (int i=0; i < SM_NP1; i++) {
            if (i != ilo) {
              for (int j = 0; j < SM_NP; j++) {
                psum[j] = 0.5 * (xsmplx[i][j] + xsmplx[ilo][j]);
                xsmplx[i][j] = psum[j];
              }
              ysearch[i] = chi_squared(psum);
            }
          }
          iter += SM_NP;
          // GO TO 1
          loop2 = false;
        }
      } else {
        //mprintf("\tCase 3\n");
        iter--;
      }
      // GO TO 2
    } // END LOOP 2
  } // END LOOP 1
  return iter;
}

// Average_vertices()
void Action_Rotdif::Average_vertices(Vec6& xsearch, double xsmplx[SM_NP1][SM_NP]) 
{
    for (int j=0; j < SM_NP; j++) {
      xsearch[j] = 0;
      for (int k = 0; k < SM_NP1; k++) 
        xsearch[j] += xsmplx[k][j];
      xsearch[j] /= SM_NP1;
    }
}

// Action_Rotdif::Simplex_min()
/** Main driver routine for Amoeba (downhill simplex) minimizer. In the 
  * simplex method, N+1 initial points (where N is the dimension of the 
  * search space) must be chosen; the SVD solution provides one of these. 
  * Initial points are stored in rows of xsmplx. Components of vector 
  * delq should be of the order of the characteristic "lengthscales" over 
  * which the Q tensor varies. delqfrac determines the size of variation
  * for each of the components of Q; the sign of the variation is randomly 
  * chosen.
  */
int Action_Rotdif::Simplex_min(Vec6& Q_vector) {
  Vec6 xsearch;
  double ysearch[SM_NP1];
  double xsmplx[SM_NP1][SM_NP];
  const int nsearch = 1;
  double sgn;
  //int test_seed = -3001796; // For tensorfit_ comparison

  mprintf("\tDetermining diffusion tensor with full anisotropy.\n");
  // Allocate tau1, tau2, and sumc2; used in calc_Asymmetric
  tau1_.resize(nvecs_);
  tau2_.resize(nvecs_);
  sumc2_.resize(nvecs_);
  // Set Tau to be used in chi_squared based on olegendre
  if (olegendre_ == 1)
    Tau_ = &tau1_;
  else if (olegendre_ == 2)
    Tau_ = &tau2_;
  else
    Tau_ = 0; // Should never get here
  // First, back-calculate with the SVD tensor, but with the full anisotropy
  // chi_squared performs diagonalization. The workspace for dsyev should
  // already have been set up in Tensor_Fit.
  outfile_.Printf("Same diffusion tensor, but full anisotropy:\n");
  outfile_.Printf("  chi_squared for SVD tensor is %15.5lf\n",chi_squared(Q_vector));
  outfile_.Printf("     taueff(obs) taueff(calc)\n");
  for (int i = 0; i < nvecs_; i++) 
    outfile_.Printf("%5i%10.5lf%10.5lf%10.5lf\n",i+1,D_eff_[i],(*Tau_)[i],sumc2_[i]);
  outfile_.Printf("\n");

  // Now execute the simplex search method with initial vertices,
  // xsimplx, and chi-squared values, ysearch.
  // We restart the minimization a pre-set number of times to avoid
  // an anomalous result.
  for (int n = 0; n < SM_NP; n++) xsearch[n] = Q_vector[n];

  // BEGIN SIMPLEX LOOP
  for (int i = 0; i < nsearch; i++) {
    for (int j = 0; j < SM_NP; j++) xsmplx[0][j] = xsearch[j];

    for (int j = 0; j < SM_NP; j++) {
      for (int k = 0; k < SM_NP; k++) {
        if (j==k) {
          sgn = RNgen_.rn_gen() - 0.5; // -0.5 <= sgn <= 0.5
          //sgn = random_(test_seed) - 0.5; // For tensorfit_ comparison
          if (sgn < 0)
            sgn = -1.0;
          else
            sgn = 1.0;
          xsmplx[j+1][k] = xsmplx[0][k] * (1+(sgn*delqfrac_));
        } else {
          xsmplx[j+1][k] = xsmplx[0][k];
        }
      }
    }

    // As to amoeba, chi-squared must be evaluated for all
    // vertices in the initial simplex.
    for (int j=0; j < SM_NP1; j++) {
      for (int k=0; k < SM_NP; k++) {
        xsearch[k] = xsmplx[j][k];
      }
      ysearch[j] = chi_squared(xsearch);
    }

    // Average the vertices and compute details of the average.
    /*for (int j=0; j < SM_NP; j++) {
      xsearch[j] = 0;
      for (int k = 0; k < SM_NP1; k++) 
        xsearch[j] += xsmplx[k][j];
      xsearch[j] /= SM_NP1;
    }*/
    Average_vertices( xsearch, xsmplx );
    sgn = chi_squared( xsearch );
    Vec3 d_props = calculate_D_properties(D_XYZ_);

    outfile_.Printf("Input to amoeba - average at cycle %i\n",i+1);
    outfile_.Printf("    Initial chisq = %15.5lf\n",sgn);
    PrintVector(outfile_,"Dav, aniostropy, rhombicity:",d_props);
    PrintVector(outfile_,"D tensor eigenvalues:",D_XYZ_);
    PrintMatrix(outfile_,"D tensor eigenvectors (in columns):",D_tensor_);

    Amoeba( xsmplx, ysearch );

    // Put amoeba results into xsearch
    for (int j=0; j < SM_NP1; j++) {
      for (int k=0; k < SM_NP; k++) {
        xsearch[k] = xsmplx[j][k];
      }
      ysearch[j] = chi_squared(xsearch);
    }

    // Average the vertices and compute details of the average.
    Average_vertices( xsearch, xsmplx );
    sgn = chi_squared( xsearch );
    d_props = calculate_D_properties(D_XYZ_);

    outfile_.Printf("Output from amoeba - average at cycle %i\n",i+1);
    outfile_.Printf("    Final chisq = %15.5lf\n",sgn);
    PrintVector(outfile_,"Dav, aniostropy, rhombicity:",d_props);
    PrintVector(outfile_,"D tensor eigenvalues:",D_XYZ_);
    PrintMatrix(outfile_,"D tensor eigenvectors (in columns):",D_tensor_);
    outfile_.Printf("     taueff(obs) taueff(calc)\n");
    for (int i = 0; i < nvecs_; i++)
      outfile_.Printf("%5i%10.5lf%10.5lf%10.5lf\n",i+1,D_eff_[i],(*Tau_)[i],sumc2_[i]);
    outfile_.Printf("\n");
   
    // cycle over main loop, but first reduce the size of delqfrac:
    delqfrac_ *= 0.750;
    if (debug_>0)mprintf("\tAmoeba: Setting delqfrac to %15.7lf\n",delqfrac_);
  }

  // Set q vector to the final average result from simpmin
  for (int n = 0; n < SM_NP; n++) Q_vector[n] = xsearch[n];

  return 0;
}
#undef SM_NP
#undef SM_NP1
// -----------------------------------------------------------------------------
// Action_Rotdif::Grid_search()
/** Given a Q tensor, search for a better Q tensor by simple
  * grid search on all 6 elements. If a better solution is found
  * store it.
  * \return 1 if grid search found a better solution.
  * \return 0 if no better solution was found.
  */
int Action_Rotdif::Grid_search(Vec6& Q_vector, int gridsize) {
  Vec6 xsearch;
  Vec6 best;
  double sgn, sgn0;
  bool success = false;
  if (gridsize < 1) {
    mprinterr("Error: Rotdif: Grid search: Grid < 1\n");
    return 0;
  }
  int gridmax = gridsize + 1;
  int gridmin = -gridsize;

  sgn0 = chi_squared(Q_vector);
  mprintf("Grid search: Starting chisq is %15.5lf\n",sgn0);
  //mprintf("_Grid q0 %10.5lf\n",Q_vector[0]);

  // Store initial solution
  for (int b = 0; b < 6; b++) best[b] = Q_vector[b];

  for (int i = gridmin; i < gridmax; i++) {
    xsearch[0] = Q_vector[0] + (i*delqfrac_/100.0);
    for (int j = gridmin; j < gridmax; j++) {
      xsearch[1] = Q_vector[1] + (j*delqfrac_/100.0);
      for (int k = gridmin; k < gridmax; k++) {
        xsearch[2] = Q_vector[2] + (k*delqfrac_/100.0);
        for (int l = gridmin; l < gridmax; l++) {
          xsearch[3] = Q_vector[3] + (l*delqfrac_/100.0);
          for (int m = gridmin; m < gridmax; m++) {
            xsearch[4] = Q_vector[4] + (m*delqfrac_/100.0);
            for (int n = gridmin; n < gridmax; n++) {
              xsearch[5] = Q_vector[5] + (n*delqfrac_/100.0);

              sgn = chi_squared( xsearch );
              //mprintf("_Grid x0 %10.5lf\n",xsearch[0]);
              //mprintf("_Grid %3i%3i%3i%3i%3i%3i%10.5lf%10.5lf\n",
              //        i,j,k,l,m,n,sgn,sgn0);
              if (sgn < sgn0) {
                for (int b = 0; b < 6; b++) best[b] = xsearch[b];
                sgn0 = sgn;
                success = true;
              }
            } // n
          } // m
        } // l
      } // k
    } // j
  } // i
      
  // Set best solution
  if (success) {
    mprintf("  Grid search succeeded.\n");
    for (int b = 0; b < 6; b++) Q_vector[b] = best[b];
    return 1;
  }
  return 0;
}

static void printMatrix(const char *Title, const double *U, int mrows, int ncols) {
  mprintf("    %s",Title);
  int usize = mrows * ncols;
  for (int i = 0; i < usize; i++) {
    if ( (i%ncols)==0 ) mprintf("\n");
    mprintf(" %10.5lf",U[i]);
  }
  mprintf("\n");
}

// Action_Rotdif::Tensor_Fit()
/** Based on random_vectors and effective diffusion constants D_eff previously
  * calculated, first find the tensor Q (and therefore D) in the small
  * anisotropic limit by solving:
  *   D_eff(n) = At(n) * Q
  * where At(n) is composed of D_eff vector components:
  *   { x^2, y^2, z^2, 2xy, 2yz, 2xz }
  * \param vector_q Will be set with Q tensor for small anisotropic limit.
  */
int Action_Rotdif::Tensor_Fit(Vec6& vector_q) {
#ifdef NO_MATHLIB
  return 1;
#else
  int info;
  double wkopt;
  //double cut_ratio = 0.000001; // threshold ratio for removing small singular values in SVD

  mprintf("\tDetermining diffusion tensor with small anisotropy.\n");
  // Generate matrix At
  // NOTE: The LAPACK fortran routines are COLUMN MAJOR, so m and n must be 
  //       flipped, i.e. matrix At must be tranposed before passing it in
  //       so we actually need to construct matrix A for SVD.
  // NOTE: LAPACK SVD routine destroys matrix A which is needed later, so 
  //       create a non-flipped version (i.e. matrix At) as well.
  int m_rows = nvecs_;
  int n_cols = 6;
  double *matrix_A = new double[ m_rows * n_cols ];
  double *matrix_At = new double[ m_rows * n_cols ];
  // Pre-compute array offsets of matrix_A for performing transpose
  int A1 = m_rows;
  int A2 = m_rows * 2;
  int A3 = m_rows * 3;
  int A4 = m_rows * 4;
  int A5 = m_rows * 5;
  double* At = matrix_At;
  int nvec = 0;
  for (std::vector<Vec3>::iterator randvec = random_vectors_.begin();
                                   randvec != random_vectors_.end(); ++randvec)
  {
    // Transpose of matrix A
    double *A = matrix_A + nvec;
    A[ 0] = (*randvec)[0] * (*randvec)[0];     // x^2
    A[A1] = (*randvec)[1] * (*randvec)[1];     // y^2
    A[A2] = (*randvec)[2] * (*randvec)[2];     // z^2
    A[A3] = 2*((*randvec)[0] * (*randvec)[1]); // 2xy
    A[A4] = 2*((*randvec)[1] * (*randvec)[2]); // 2yz
    A[A5] = 2*((*randvec)[0] * (*randvec)[2]); // 2xz
    // Matrix A
    At[0] = A[ 0];
    At[1] = A[A1];
    At[2] = A[A2];
    At[3] = A[A3];
    At[4] = A[A4];
    At[5] = A[A5];
    At += 6;
    ++nvec;
  }
  if (debug_>1) {
    printMatrix("matrix_A",matrix_A,n_cols,m_rows);
    printMatrix("matrix_At",matrix_At,m_rows,n_cols);
  }

  // Perform SVD on matrix At to generate U, Sigma, and Vt
  int lda = m_rows;
  int ldu = m_rows;
  int ldvt = n_cols;
  // NOTE: Final dimension of matrix_S is min(m,n)
  int s_dim = m_rows;
  if (n_cols < m_rows) s_dim = n_cols;
  double *matrix_S = new double[ s_dim ];
  double *matrix_U = new double[ m_rows * m_rows ];
  double *matrix_Vt = new double[ n_cols * n_cols ];
  lwork_ = -1;
  // Allocate Workspace
  dgesvd_((char*)"All",(char*)"All",m_rows, n_cols, matrix_A, lda, 
          matrix_S, matrix_U, ldu, matrix_Vt, ldvt, &wkopt, lwork_, info );
  lwork_ = (int)wkopt;
  work_ = new double[ lwork_ ];
  // Compute SVD
  dgesvd_((char*)"All",(char*)"All", m_rows, n_cols, matrix_A, lda, 
          matrix_S, matrix_U, ldu, matrix_Vt, ldvt, work_, lwork_, info );
  // matrix_A and work no longer needed
  delete[] matrix_A;
  delete[] work_;
  // DEBUG - Print Sigma
  if (debug_>0) {
    for (int i = 0; i < s_dim; i++) 
      mprintf("Sigma %6i%12.6lf\n",i+1,matrix_S[i]);
  }
  // Check for convergence
  if ( info > 0 ) {
    mprinterr( "The algorithm computing SVD of At failed to converge.\n" );
    delete[] matrix_At;
    delete[] matrix_U;
    delete[] matrix_S;
    delete[] matrix_Vt;
    return 1;
  }

  // DEBUG: Print U and Vt
  // NOTE: U and Vt are in column-major order from fortran routine
  //       so are currently implicitly transposed.
  if (debug_>1) {
    printMatrix("matrix_Ut",matrix_U,m_rows,m_rows);
    printMatrix("matrix_V",matrix_Vt,n_cols,n_cols);
  }

  // Remove small singular values from SVD
/*  double wmax = 0;
  for (int i=0; i < n; i++)
    if (matrix_S[i] > wmax) wmax=matrix_S[i];
  double wmin=wmax*cut_ratio;
  for (int i=0; i < n; i++)
    if (matrix_S[i] < wmin) matrix_S[i]=0;*/

  // Calculate x = V * Sigma^-1 * Ut * b
  // Take advantage of the fact that Vt and U have already been implicitly 
  // transposed in dgesvd to do everything in one step.
  // First invert Sigma.
  for (int s = 0; s < s_dim; s++) 
    if (matrix_S[s] > 0) matrix_S[s] = 1 / matrix_S[s];
  int v_idx = 0;
  for (int n = 0; n < 6; n++) {
    vector_q[n] = 0.0;
    for (int m = 0; m < m_rows; m++) {
      double vsu_sum = 0.0;
      //mprintf("Row %i (",m);
      for (int vs = 0; vs < s_dim; vs++) {
        vsu_sum += (matrix_Vt[v_idx+vs] * matrix_S[vs] * matrix_U[(vs*m_rows)+m]);
        //mprintf(" VSU v=%i s=%i u=%i",v_idx+vs,vs,(vs*m_rows)+m);
      }
      //mprintf(") m=%i\n",m);
      vector_q[n] += (vsu_sum * D_eff_[m]);
    }
    v_idx += 6;
  }
  // matrix_U, matrix_Vt and matrix_S no longer needed
  delete[] matrix_S;
  delete[] matrix_Vt;
  delete[] matrix_U;

  outfile_.Printf("Results of small anisotropy (SVD) analysis:\n");
  // Print Q vector
  PrintVec6(outfile_,"Qxx Qyy Qzz Qxy Qyz Qxz",vector_q);
  // ---------------------------------------------

  // Convert vector Q to diffusion tensor D
  vector_q.Q_to_D(D_tensor_);
  PrintMatrix(outfile_,"D_tensor",D_tensor_);

  // Save D for later use in back calculating deff from Q
  Matrix_3x3 matrix_D_local = D_tensor_;

  // Diagonalize D to find eigenvalues and eigenvectors
  // (principal components and axes)
  // Determine workspace. Do not delete work after set up since the 
  // Action_Rotdif::chi_squared routine will use dsyev_ for diagnolization.
  lwork_ = -1;
  n_cols = 3;
  dsyev_((char*)"Vectors",(char*)"Upper", n_cols, D_tensor_.Dptr(), n_cols, 
         D_XYZ_.Dptr(), &wkopt, lwork_, info);
  lwork_ = (int) wkopt;
  work_ = new double[ lwork_ ];
  // Diagonalize D_tensor
  dsyev_((char*)"Vectors",(char*)"Upper", n_cols, D_tensor_.Dptr(), n_cols, 
         D_XYZ_.Dptr(), work_, lwork_, info);
  // Check for convergence
  if (info > 0) {
    mprinterr("The algorithm computing the eigenvalues/eigenvectors of D failed to converge.\n");
    delete[] work_;
    delete[] matrix_At;
    return 1;
  }
  //delete[] work;
  // eigenvectors are stored in columns due to implicit transpose from fortran,
  // i.e. Ex = {D_tensor[0], D_tensor[3], D_tensor[6]} etc.
  // Transpose back.
  //matrix_transpose_3x3( D_tensor );

  // Print eigenvalues/eigenvectors
  PrintVector(outfile_,"D eigenvalues",D_XYZ_);
  PrintMatrix(outfile_,"D eigenvectors (in columns)",D_tensor_);

  // Calculate Dav, Daniso, Drhomb
  Vec3 d_props = calculate_D_properties(D_XYZ_);
  PrintVector(outfile_,"Dav, Daniso, Drhomb",d_props);

  // Back-calculate the local diffusion constants via At*Q=Deff
  // First convert original D back to Q
  Vec6 vector_q_local;
  vector_q_local.D_to_Q(matrix_D_local);
  if (debug_>0) {
    mprintf("    D_to_Q\n %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
            vector_q_local[0], vector_q_local[1], vector_q_local[2],
            vector_q_local[3], vector_q_local[4], vector_q_local[4]);
  }
  std::vector<double> deff_local;
  deff_local.reserve( nvecs_ );
  At = matrix_At;
  // At*Q
  for (int i=0; i < nvecs_; i++) {
    deff_local.push_back( (At[0] * vector_q_local[0]) +
                          (At[1] * vector_q_local[1]) +
                          (At[2] * vector_q_local[2]) +
                          (At[3] * vector_q_local[3]) +
                          (At[4] * vector_q_local[4]) +
                          (At[5] * vector_q_local[5])   );
    At += 6;
  }
  // Convert deff to tau, Output
  outfile_.Printf("     taueff(obs) taueff(calc)\n");
  double sgn = 0;
  for (int i = 0; i < nvecs_; i++) {
    // For the following chisq fits, convert deff to taueff
    D_eff_[i] = 1 / (6 * D_eff_[i]);
    deff_local[i] = 1 / (6 * deff_local[i]);
    outfile_.Printf("%5i%10.5lf%10.5lf\n", i+1, D_eff_[i], deff_local[i]);
    // NOTE: in rotdif code, sig is 1.0 for all nvecs 
    double diff = deff_local[i] - D_eff_[i];
    sgn += (diff * diff);
  }
  outfile_.Printf("  chisq for above is %15.5lf\n\n",sgn);

  // Cleanup
  delete[] matrix_At;

  // DEBUG -------------------------------------------------
/* 
  // CHECK: Ensure U is orthogonal (U * Ut) = I
  double *Ut = matrix_transpose(matrix_U,m_rows,m_rows);
  printMatrix("Ut",Ut,m_rows,m_rows);
  double *UtU = new double[ m_rows * m_rows ];
  dgemm_((char*)"N",(char*)"N",m_rows,m_rows,m_rows,alpha,Ut,m_rows,
         matrix_U,m_rows,beta,UtU,m_rows);
  printMatrix("UtU",UtU,m_rows,m_rows);
  delete[] UtU;
  delete[] Ut;
*/
/*
  // CHECK: Ensure that U * S * Vt = A
  double *diag_S = new double[ m_rows * n_cols ];
  int sidx=0;
  int diag_sidx=0;
  for (int m=0; m < m_rows; m++) {
    for (int n=0; n < n_cols; n++) {
      if (m==n) 
        diag_S[diag_sidx]=matrix_S[sidx++];
      else
        diag_S[diag_sidx]=0.0;
      ++diag_sidx;
    }
  }
  printMatrix("diag_S",diag_S,m_rows,n_cols);
  // Transpose diag_S for fortran; U and Vt are already transposed (result from fortran SVD)
  double *diag_St = matrix_transpose( diag_S, m_rows, n_cols );
  printMatrix("diag_St",diag_St,n_cols,m_rows);
  // U * S
  double *US = new double[ m_rows * n_cols ];
  dgemm_((char*)"N",(char*)"N",m_rows,n_cols,m_rows,alpha,matrix_U,m_rows,
         diag_St,m_rows,beta,US,m_rows);
  // US * Vt
  double *USVt = new double[ m_rows * n_cols ];
  dgemm_((char*)"N",(char*)"N",m_rows,n_cols,n_cols,alpha,US,m_rows,
         matrix_Vt,n_cols,beta,USVt,m_rows);
  // The result from fortran will actually be (USVt)t
  //double *USVtt = matrix_transpose( USVt, m_rows, n_cols );
  printMatrix("matrix_A",matrix_A,m_rows,n_cols);
  printMatrix("USVt",USVt,n_cols,m_rows);
  // ----- Cleanup
  delete[] diag_S;
  delete[] diag_St;
  delete[] US;
  delete[] USVt;
*/
/*  
  // ---------------------------------------------
  // CHECK: Calculate x = V * Sigma^-1 * Ut * b with BLAS routines
  double alpha = 1.0;
  double beta = 0.0;
  // First create transposed (for fortran) n x m sigma^-1
  double *matrix_St = new double[ n_cols * m_rows];
  int ndiag = n_cols+1;
  int sidx = 0;
  for (int i = 0; i < n_cols*m_rows; i++) {
    if ( (i%ndiag)==0 && sidx<s_dim) 
      matrix_St[i] = (1 / matrix_S[sidx++]);
    else
      matrix_St[i] = 0;
  }
  printMatrix("matrix_St",matrix_St,m_rows,n_cols);
  // 1) V * Sigma^-1
  //    Call dgemm with first arg "T" to indicate we want matrix_Vt to
  //    be transposed. 
  double *matrix_VS = new double[ n_cols * m_rows ];
  dgemm_((char*)"T",(char*)"N",n_cols,m_rows,n_cols,alpha,matrix_Vt,n_cols,
         matrix_St,n_cols,beta,matrix_VS,n_cols);
  // matrix_S, matrix_Vt, and matrix_St no longer needed
  delete[] matrix_S;
  delete[] matrix_Vt;
  delete[] matrix_St;
  // 2) VSigma^-1 * Ut
  //    Call dgemm with second arg "T" to indicate we want matrix_U to
  //    be transposed.
  double *matrix_VSUt = new double[ n_cols * m_rows ];
  dgemm_((char*)"N",(char*)"T",n_cols,m_rows,m_rows,alpha,matrix_VS,n_cols,
         matrix_U,m_rows,beta,matrix_VSUt,n_cols);
  // matrix_VS and matrix_U no longer needed
  delete[] matrix_VS;
  delete[] matrix_U;
  // 3) VSigma^-1*Ut * b
  int incx = 1;
  dgemv_((char*)"N",n_cols,m_rows,alpha,matrix_VSUt,n_cols,deff,incx,beta,vector_q,incx);
  // matrix_VSUt no longer needed
  delete[] matrix_VSUt;
  // ---------------------------------------------
*/
  // END DEBUG ---------------------------------------------

  return 0;
#endif
}

// Action_Rotdif::DetermineDeffs()
/** Calculate effective diffusion constant for each random vector. 
  * Vectors will be normalized during this phase. First rotate the vector 
  * by all rotation matrices, storing the resulting vectors. The first entry 
  * of rotated_vectors is the original random vector. Then calculate the time 
  * correlation function of the vector. Finally compute the area under the 
  * time correlation function curve and estimate the diffusion constant.
  * Sets D_Eff, normalizes random_vectors.
  */
// TODO: OpenMP Parallelize
int Action_Rotdif::DetermineDeffs() {
  int itotframes;                 // Total number of frames (rotation matrices) 
  DataSet_Vector rotated_vectors; // Hold vectors after rotation with Rmatrices
  int maxdat;                     // Length of C(t) 
  std::vector<double> pX;         // Hold X values of C(t)
  std::vector<double> p1;         // Hold Y values of C(t) for p1
  std::vector<double> p2;         // Hold Y values of C(t) for p2
  std::vector<double>* pY;        // Point to p1 or p2 depending on olegendre
  int meshSize;                   // Total mesh size, maxdat * NmeshPoints

  mprintf("\tDetermining local diffusion constants for each vector.\n");
  ProgressBar progress( nvecs_ );

  itotframes = (int) Rmatrices_.size();
  if (ncorr_ == 0) ncorr_ = itotframes;
  maxdat = ncorr_ + 1;
  // Allocate memory to hold calcd effective D values
  D_eff_.reserve( nvecs_ );
  // Allocate memory to hold rotated vectors. Need +1 since the original
  // vector is stored at position 0. 
  rotated_vectors.Allocate( itotframes + 1 );
  // Allocate memory for C(t)
  p1.reserve( maxdat );
  p2.reserve( maxdat );
  pX.reserve( maxdat );
  // Set X values of C(t) based on tfac
  for (int i = 0; i < maxdat; i++)
    pX.push_back( (double)i * tfac_ );
  // Allocate mesh to hold interpolated C(t)
  if (NmeshPoints_ < 1)
    meshSize = maxdat * 2;
  else
    meshSize = maxdat * NmeshPoints_;
  // spline will be used to interpolate C(t) for better intergration
  Interpolate spline( ti_, tf_, meshSize );
  // Check order of legendre polynomial to determine whether we use p1 or p2
  if (olegendre_==1)
    pY = &p1;
  else if (olegendre_==2)
    pY = &p2;
  else
    return 1; // Should never get here, olegendre is checked in init
  // LOOP OVER RANDOM VECTORS
  int nvec = 0;
  for (std::vector<Vec3>::iterator rndvec = random_vectors_.begin();
                                   rndvec != random_vectors_.end(); ++rndvec)
  {
    progress.Update( nvec );
    // Reset rotated_vectors to the beginning 
    rotated_vectors.reset();
    // Normalize vector
    (*rndvec).Normalize();
    // Assign normalized vector to rotated_vectors position 0
    rotated_vectors.AddVxyz( *rndvec );
   
    // Loop over rotation matrices
    for (std::vector<Matrix_3x3>::iterator rmatrix = Rmatrices_.begin();
                                           rmatrix != Rmatrices_.end();
                                           ++rmatrix)
    {
      // Rotate normalized vector
      rotated_vectors.AddVxyz( *rmatrix * (*rndvec) );
      // DEBUG
      //Vec3 current = rotated_vectors.CurrentVec();
      //mprintf("DBG:Rotated %6u: %15.8f%15.8f%15.8f\n", rmatrix - Rmatrices_.begin(),
      //        current[0], current[1], current[2]); 
    }
    // Calculate time correlation function for this vector
    if (usefft_)
      fft_compute_corr(rotated_vectors, maxdat, *pY, olegendre_);
    else
      compute_corr(rotated_vectors, maxdat, p2, p1);
    // Calculate mesh Y values
    spline.SetMesh_Y(pX, *pY);
    // Integrate
    double integral = spline.Integrate_Trapezoid();
    // Solve for deff
    D_eff_.push_back( calcEffectiveDiffusionConst(integral) );

    // DEBUG: Write out p1 and p2 ------------------------------------
    if (!corrOut_.empty() || debug_ > 3) {
        CpptrajFile outfile;
        std::string namebuffer;
        if (!corrOut_.empty())
          namebuffer = NumberFilename( corrOut_, nvec );
        else
          namebuffer = NumberFilename( "p1p2.dat", nvec );
        outfile.OpenWrite(namebuffer);
        for (int i = 0; i < maxdat; i++) 
          //outfile.Printf("%lf %lf %lf\n",pX[i], p2[i], p1[i]);
          outfile.Printf("%f %f\n",pX[i], (*pY)[i]);
        outfile.CloseFile();
        //    Write Mesh
        if (debug_>3) {
          namebuffer = NumberFilename( "mesh.dat", nvec );
          outfile.OpenWrite(namebuffer);
          for (int i=0; i < spline.Mesh_Size(); i++)
            outfile.Printf("%f %f\n", spline.X(i), spline.Y(i));
          outfile.CloseFile();
        }
        if (!corrOut_.empty()) 
          corrOut_.clear();
    }
    if (debug_ > 0) {
      mprintf("DBG: Vec %i Spline integral= %12.4lf\n",nvec,integral);
      mprintf("DBG: deff is %lf\n",D_eff_[nvec]);
    }
    // END DEBUG -----------------------------------------------------
    ++nvec;
  }

  return 0;
}

// Action_Rotdif::Print()
/** Main tensorfit calculation.
  * - Read/generate random vectors; analogous to e.g. N-H bond vectors.
  * - For each random vector:
  *   a) Normalize it.
  *   b) Rotate it by each of the rotation matrices.
  *   c) Calculate the time correlation function of the rotated vector.
  *   d) Integrate the time correlation function to obtain F.
  *   e) Given f, iteratively solve eq 18  in Wong&Case 2008 for
  *      effective value of diffusion constant (deff) for that vector
  * - Given deff for each vector, solve for Q assuming small anisotropy
  *   via SVD using eq 13 from Wong&Case 2008
  * - Based on Q from small anisotropic limit, use downhill simplex
  *   minimizer to optimize Q in full anisotropic limit
  */
void Action_Rotdif::Print() {

  mprintf("    ROTDIF:\n");
 
  // Read/Generate nvecs random vectors
  random_vectors_ = RandomVectors();
  if (random_vectors_.empty()) return;

  // If no rotation matrices generated, exit
  if (Rmatrices_.empty()) return;
  // HACK: To match results from rmscorr.f (where rotation matrices are
  //       implicitly transposed), transpose each rotation matrix.
  // NOTE: Is this actually correct? Want inverse rotation?
  for (std::vector<Matrix_3x3>::iterator rmatrix = Rmatrices_.begin();
                                      rmatrix != Rmatrices_.end(); rmatrix++) 
    (*rmatrix).Transpose();
  // Print rotation matrices
  if (!rmOut_.empty()) {
    CpptrajFile rmout;
    if (rmout.SetupWrite(rmOut_,debug_)) {
      mprinterr("    Error: Rotdif: Could not set up %s for writing.\n",rmOut_.c_str());
    } else {
      rmout.OpenFile();
      int rmframe=1;
      for (std::vector<Matrix_3x3>::iterator rmatrix = Rmatrices_.begin();
                                             rmatrix != Rmatrices_.end(); rmatrix++) 
      {
        rmout.Printf("%13i %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f\n",
             rmframe++,
            (*rmatrix)[0], (*rmatrix)[1], (*rmatrix)[2],
            (*rmatrix)[3], (*rmatrix)[4], (*rmatrix)[5],
            (*rmatrix)[6], (*rmatrix)[7], (*rmatrix)[8]);
      }
      rmout.CloseFile();
    }
  }
  mprintf("\t%i vectors, %u rotation matrices.\n",nvecs_,Rmatrices_.size());

  // Determine effective D for each vector
  DetermineDeffs( );
  // Print deffs
  if (!deffOut_.empty()) {
    CpptrajFile dout;
    if (dout.SetupWrite(deffOut_,debug_)) {
      mprinterr("    Error: Rotdif: Could not set up file %s\n",deffOut_.c_str());
    } else {
      dout.OpenFile();
      for (int vec = 0; vec < nvecs_; vec++)
        dout.Printf("%6i%15.8lf\n",vec+1,D_eff_[vec]);
      dout.CloseFile();
    }
  }

  // All remaining functions require LAPACK
# ifndef NO_MATHLIB
  Vec6 Q_isotropic;
  Tensor_Fit( Q_isotropic );

  // Using Q (small anisotropy) as a guess, calculate Q with
  // full anisotropy
  Vec6 Q_anisotropic;
  Q_anisotropic[0] = Q_isotropic[0];
  Q_anisotropic[1] = Q_isotropic[1];
  Q_anisotropic[2] = Q_isotropic[2];
  Q_anisotropic[3] = Q_isotropic[3];
  Q_anisotropic[4] = Q_isotropic[4];
  Q_anisotropic[5] = Q_isotropic[5];
  Simplex_min( Q_anisotropic );

  // Brute force grid search
  if (do_gridsearch_)
    Grid_search( Q_anisotropic, 5 );
# endif
}
