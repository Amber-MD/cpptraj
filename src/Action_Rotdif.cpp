// Rotdif
#include <cmath>
#include <cfloat> // DBL_MAX
#include <cstdio> //sscanf
#include "Action_Rotdif.h"
#include "CpptrajStdio.h"
#include "Constants.h" // TWOPI
#include "vectormath.h"
#include "Integrate.h"
#include "ProgressBar.h"

// Definition of Fortran subroutines in Rotdif.f called from this class
#ifndef NO_PTRAJ_ANALYZE
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
#endif

// CONSTRUCTOR
Rotdif::Rotdif() :
  rseed_( 1 ),
  nvecs_( 0 ),
  tfac_( 0.0 ),
  ti_( 0.0 ),
  tf_( 0.0 ),
  itmax_( 0 ),
  delmin_( 0.0 ),
  d0_( 0.0 ),
  olegendre_( 2 ),
  ncorr_( 0 ),
  delqfrac_( 0 ),
  amoeba_ftol_(0.0000001 ),
  amoeba_itmax_(10000 ),
  do_gridsearch_( false ),
  work_( NULL ),
  lwork_( 0 ),
  randvecOut_( NULL ),
  randvecIn_( NULL ),
  rmOut_( NULL ),
  deffOut_( NULL ),
  random_vectors_( NULL ),
  D_eff_( NULL ),
  Tau_( NULL )
{
  //fprintf(stderr,"Rotdif Con\n") ),
} 

// DESTRUCTOR
Rotdif::~Rotdif() {
  //fprintf(stderr,"Rotdif Destructor.\n");
  for (std::vector<double*>::iterator Rmatrix = Rmatrices_.begin();
                                      Rmatrix != Rmatrices_.end();
                                      Rmatrix++)
    delete[] *Rmatrix;
  if (work_!=NULL) delete[] work_;
  if (random_vectors_!=NULL) delete[] random_vectors_;
  if (D_eff_!=NULL) delete[] D_eff_;
  // Close output file
  outfile_.CloseFile();
}

// Rotdif::init()
/** Expected call: rotdif [rseed <rseed>] [nvecs <nvecs>]  
  *                       ref <refname> | refindex <refindex> | reference
  *                       [<refmask>] [ncorr <ncorr>] dt <tfac> [ti <ti>] tf <tf>
  *                       [itmax <itmax>] [tol <delmin>] [d0 <d0>] [order <olegendre>]
  *                       [delqfrac <delqfrac>] [rvecout <randvecOut>]
  *                       [rmout <rmOut>] [deffout <deffOut>] [outfile <outfilename>]
  *                       [corrout <corrOut>]
  *                       [rvecin <randvecIn>]
  *                       [gridsearch]
  */ 
// Dataset name will be the last arg checked for. Check order is:
//    1) Keywords
//    2) Masks
//    3) Dataset name
int Rotdif::init( ) {
  int refindex;
  char *referenceName, *mask0, *outfilename;
  double Trans[3]; // Dummy variable for CenteredRef routine
  // Get Keywords
  nvecs_ = actionArgs.getKeyInt("nvecs",1000);
  rseed_ = actionArgs.getKeyInt("rseed",80531);
  ncorr_ = actionArgs.getKeyInt("ncorr",0);
  tfac_ = actionArgs.getKeyDouble("dt",0);
  if (tfac_<=0) {
    mprinterr("    Error: Rotdif: dt (timestep) must be specified and > 0.\n");
    return 1;
  }
  ti_ = actionArgs.getKeyDouble("ti",0);
  tf_ = actionArgs.getKeyDouble("tf",0);
  if (tf_ <= ti_) {
    mprinterr("    Error: Rotdif: Initial time ti (%lf) must be < final time tf (%lf).\n",
              ti_, tf_);
    return 1;
  }
  itmax_ = actionArgs.getKeyInt("itmax",500);
  delmin_ = actionArgs.getKeyDouble("tol",0.000001);
  d0_ = actionArgs.getKeyDouble("d0",0.03);
  olegendre_ = actionArgs.getKeyInt("order",2);
  if (olegendre_!=1 && olegendre_!=2) {
    mprinterr("    Error: Rotdif: Order of legendre polynomial (%i) must be 1 or 2.\n",
              olegendre_);
    return 1;
  }
  delqfrac_ = actionArgs.getKeyDouble("delqfrac",0.5);
  randvecOut_ = actionArgs.getKeyString("rvecout",NULL);
  randvecIn_ = actionArgs.getKeyString("rvecin",NULL);
  rmOut_ = actionArgs.getKeyString("rmout",NULL);
  deffOut_ = actionArgs.getKeyString("deffout",NULL);
  outfilename = actionArgs.getKeyString("outfile",NULL);
  corrOut_ = actionArgs.GetStringKey("corrout");
  do_gridsearch_ = actionArgs.hasKey("gridsearch");

  referenceName=actionArgs.getKeyString("ref",NULL);
  refindex=actionArgs.getKeyInt("refindex",-1);
  if (actionArgs.hasKey("reference")) // For compatibility with ptraj
    refindex = 0;

  // Get Masks
  mask0 = actionArgs.getNextMask();
  AtomMask RefMask(mask0);
  TargetMask_.SetMaskString(mask0);

  // Dataset
  
  // Initialize random number generator
  RNgen_.rn_set( rseed_ );

  // Set up reference for RMSD
  // Attempt to get reference index by name
  if (referenceName!=NULL)
    refindex=FL->FindName(referenceName);
  // Get reference frame by index
  Frame* TempFrame=FL->GetFrame(refindex);
  if (TempFrame==NULL) {
    mprinterr("    Error: Rotdif::init: Could not get reference index %i\n",refindex);
    return 1;
  }
  //RefFrame = *TempFrame;
  // Set reference parm
  Topology* RefParm = FL->GetFrameParm(refindex);
  // Setup reference mask
  if (RefParm->SetupIntegerMask( RefMask )) return 1;
  if (RefMask.None()) {
    mprintf("    Error: Rotdif::init: No atoms in reference mask.\n");
    return 1;
  }
  // Allocate frame for selected reference atoms
  SelectedRef_.SetupFrameFromMask(RefMask, RefParm->Mass());
  // Set reference frame coordinates
  SelectedRef_.SetCoordinates(*TempFrame, RefMask);
  // Always fitting; Pre-center reference frame
  SelectedRef_.CenterReference(Trans, useMass_); 

  // Open output file. Defaults to stdout if no name specified
  if (outfile_.SetupWrite(outfilename,debug)) {
    mprinterr("Error setting up Rotdif output file.\n");
    return 1;
  }
  if (outfile_.OpenFile()) {
    mprinterr("Error opening Rotdif output file.\n");
    return 1;
  }

  mprintf("    ROTDIF: Random seed %i, # of random vectors to generate: %i\n",rseed_,nvecs_);
  mprintf("            Max length to compute time correlation functions:");
  if (ncorr_ == 0)
    mprintf(" Total # of frames\n");
  else
    mprintf(" %i\n",ncorr_);
  mprintf("            Timestep = %.4lf, T0 = %.4lf, TF = %.4lf\n",tfac_,ti_,tf_);
  mprintf("            Iterative solver: Max iterations = %i, tol = %lf, initial guess = %lf\n",
          itmax_, delmin_, d0_);
  mprintf("            Order of Legendre polynomial = %i\n",olegendre_);
  mprintf("            Simplex scaling factor=%.4lf\n",delqfrac_);
  if (do_gridsearch_)
    mprintf("            Grid search will be performed for Q with full anisotropy (time consuming)\n");
  if (randvecIn_!=NULL)
    mprintf("            Random vectors will be read from %s\n",randvecIn_);
  if (randvecOut_!=NULL)
    mprintf("            Random vectors will be written to %s\n",randvecOut_);
  if (rmOut_!=NULL)
    mprintf("            Rotation matrices will be written out to %s\n",rmOut_);
  if (deffOut_!=NULL)
    mprintf("            Deff will be written out to %s\n",deffOut_);
  if (!corrOut_.empty())
    mprintf("            Time correlation for l=1 and l=2 for vector 0 will be written to %s\n",
            corrOut_.c_str());
#ifdef NO_PTRAJ_ANALYZE
  mprintf("------------------------------------------------------\n");
  mprintf("Warning: Cpptraj was compiled with -DNO_PTRAJ_ANALYZE.\n");
  mprintf("         The final tensor fit cannot be performed.\n");
  mprintf("         Only Deffs will be calculated.\n");
  mprintf("------------------------------------------------------\n");
#else
  if (outfilename!=NULL)
    mprintf("            Diffusion constants and tau will be written to %s\n",outfilename);
  else
    mprintf("            Diffusion constants and tau will be written to STDOUT.\n");
#endif
  return 0;
}

// Rotdif::setup()
/** Determine what atoms each mask pertains to for the current parm file.
  * Also determine whether imaging should be performed.
  */
int Rotdif::setup() {

  if ( currentParm->SetupIntegerMask( TargetMask_ ) ) return 1;
  if ( TargetMask_.None() ) {
    mprintf("    Error: Rotdif::setup: No atoms in mask.\n");
    return 1;
  }
  // Allocate space for selected atoms in the frame. This will also put the
  // correct masses in based on the mask.
  SelectedTgt_.SetupFrameFromMask(TargetMask_, currentParm->Mass());
  // Check that num atoms in frame mask from this parm match ref parm mask
  if ( SelectedRef_.Natom() != TargetMask_.Nselected() ) {
    mprintf( "    Error: Number of atoms in RMS mask (%i) does not \n",TargetMask_.Nselected());
    mprintf( "           equal number of atoms in Ref mask (%i).\n",SelectedRef_.Natom());
    return 1;
  }
  
  // Print info for this parm
  mprintf("    ROTDIF: %i atoms selected for RMS fit.\n",TargetMask_.Nselected());
        
  return 0;  
}

// Rotdif::action()
/** Calculate and store the rotation matrix for frame to reference.
  */
int Rotdif::action() {
  double *U, Trans[6];

  // Set selected frame atoms. Masses have already been set.
  SelectedTgt_.SetCoordinates(*currentFrame, TargetMask_);

  U = new double[ 9 ];
  SelectedTgt_.RMSD_CenteredRef(SelectedRef_, U, Trans, useMass_);

  Rmatrices_.push_back( U );

  return 0;
} 

// ---------- ROTATIONAL DIFFUSION CALC ROUTINES -------------------------------
// Rotdif::randvec()
/** If no input file is specified by randvecIn, generate nvecs vectors of length 
  * 1.0 centered at the coordinate origin that are randomly oriented. The x,
  * y, and z components of each vector are generated from a polar coordinate 
  * system in which phi is randomly chosen in the range 0 to 2*PI and theta
  * is randomly chosen in the range 0 to PI/2 (single hemisphere).
  *   R = 1.0
  *   x = R * sin(theta) * cos(phi)
  *   y = R * sin(theta) * sin(phi)
  *   z = R * cos(theta)
  * \return Array containing vector XYZ coordinates, X0Y0Z0X1Y1Z1...
  * \return NULL on error
  */
// NOTE: Theta could also be generated in the same way as phi. Currently done
//       to be consistent with the original implementation in randvec.F90
double *Rotdif::randvec() {
  const size_t BUF_SIZE = 128;
  double *XYZ;
  int xyz_size = nvecs_ * 3;
  CpptrajFile vecIn;
  char buffer[BUF_SIZE];

  XYZ = new double[ xyz_size ];

  // ----- Read nvecs vectors from a file
  if (randvecIn_!=NULL) {
    if (vecIn.SetupRead(randvecIn_, debug)) {
      mprinterr("Error: Could not setup random vectors input file %s",randvecIn_);
      delete[] XYZ;
      return NULL;
    }
    if (vecIn.OpenFile()) {
      mprinterr("Error: Could not open random vectors input file %s",randvecIn_);
      delete[] XYZ;
      return NULL;
    }
    for (int i = 0; i < xyz_size; i+=3) {
      if (vecIn.IO->Gets(buffer, BUF_SIZE) ) {
        mprinterr("Error: Could not read vector %i from file %s\n",(i/3)+1,randvecIn_);
        delete[] XYZ;
        return NULL;
      }
      sscanf(buffer,"%*i %lf %lf %lf",XYZ+i,XYZ+i+1,XYZ+i+2);
    }
    vecIn.CloseFile();
  // ----- Generate nvecs vectors
  } else {
    for (int i = 0; i < xyz_size; i+=3) {
      double phi = TWOPI * RNgen_.rn_gen();
      // XYZ[i+2] is cos(theta)
      XYZ[i+2] = 1 - RNgen_.rn_gen();
      double theta = acos( XYZ[i+2] );
      double sintheta = sin( theta );
      XYZ[i  ] = sintheta * cos( phi );
      XYZ[i+1] = sintheta * sin( phi );
    }
  }
  // Print vectors
  if (randvecOut_!=NULL) {
    CpptrajFile rvout;
    if (rvout.SetupWrite(randvecOut_,debug)) {
      mprinterr("    Error: Rotdif: Could not set up %s for writing.\n",randvecOut_);
    } else {
      rvout.OpenFile();
      int idx = 0;
      for (int i = 1; i <= nvecs_; i++) {
        rvout.Printf("%6i  %15.8lf  %15.8lf  %15.8lf\n",i,XYZ[idx],XYZ[idx+1],XYZ[idx+2]);
        idx += 3;
      }
      rvout.CloseFile();
    }
  }

  return XYZ;
}

// Rotdif::compute_corr()
/** Given a normalized vector that has been randomly rotated itotframes 
  * times, compute the time correlation functions of the vector.
  * \param rotated_vectors array of vector coords for each frame, V0x,V0y,V0z,V1x,V1y,V1z 
  * \param maxdat Maximum length to compute time correlation functions (units of 'frames')
  * \param itotframes total number of frames provided
  * \param p2 Will be set wit values for correlation function, l=2
  * \param p1 Will be set wit values for correlation function, l=1
  */
int Rotdif::compute_corr(double *rotated_vectors, int maxdat, int itotframes, 
                         double *p2, double *p1)
{
  double *VJ, *VK;
  // Initialize p1 and p2
  for (int i = 0; i < maxdat; i++) {
    p1[i] = 0.0;
    p2[i] = 0.0;
  }

  // i loop:  each value of i is a value of delay (correlation function argument)
  // NOTE: Eventually optimize kidx and jidx?
  for (int i = 0; i < maxdat; i++) {
    int jmax = itotframes - i + 1;
    for (int j = 0; j < jmax; j++) {
      int jidx = j * 3;
      VJ = rotated_vectors + jidx;

      int k = j + i;
      int kidx = k * 3;
      VK = rotated_vectors + kidx;

      //mprintf("DBG i=%6i j=%6i k=%i\n",i,j,k);
      // Dot vector j with vector k
      double dot = dot_product( VJ, VK ); 
      p2[i] = p2[i] + (1.5*dot*dot) - 0.5;
      p1[i] = p1[i] + dot;
    }
    double one_jmax = (double) jmax;
    one_jmax = 1 / one_jmax;
    p2[i] = p2[i] * one_jmax;
    p1[i] = p1[i] * one_jmax;
  }
 
  return 0; 
}

// Rotdif::calcEffectiveDiffusionConst()
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
double Rotdif::calcEffectiveDiffusionConst(double f ) {
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
     if (debug>2)
       mprintf("ITSOLV: %6i  %15.8e  %15.8e  %15.8e\n", i,di,d,del);
     di = d;
     ++i;
  }
  if ( i>itmax_ && del>delmin_) {
     mprintf("\tWarning, itsolv did not converge: # iterations=%i, fractional change=%lf\n",
             i, del);
  } else {
    if (debug>1) mprintf("\tITSOLV Converged: # iterations=%i\n",i);
  }

  return d; 
}

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
static void Q_to_D(double *D, double *Q) {
  double tq = Q[0] + Q[1] + Q[2];
  
  D[0] = tq - (2 * Q[0]); // tq-2Qxx
  D[1] = -2 * Q[3];       // -2Qxy
  D[2] = -2 * Q[5];       // -2Qxz
  D[3] = D[1];            // -2Qyx
  D[4] = tq - (2 * Q[1]); // tq-2Qyy
  D[5] = -2 * Q[4];       // -2Qyz
  D[6] = D[2];            // -2Qzx
  D[7] = D[5];            // -2Qzy
  D[8] = tq - (2 * Q[2]); // tq-2Qzz
}

// static void D_to_Q()
/** Given diffusion tensor D[9], calculate Q[6] where Q is:
  *   Q = {Qxx, Qyy, Qzz, Qxy, Qyz, Qxz}
  * from
  *   Q = (3*Dav*I - D) / 2
  * where Dav = trace(D). See Q_to_D for more discussion.
  */
static void D_to_Q(double *Q, double *D) {
  double td = D[0] + D[4] + D[8];

  Q[0] = (td - D[0]) / 2; // Qxx
  Q[1] = (td - D[4]) / 2; // Qyy
  Q[2] = (td - D[8]) / 2; // Qzz
  Q[3] = -D[1] / 2; // Qxy
  Q[4] = -D[5] / 2; // Qyz
  Q[5] = -D[2] / 2; // Qxz
}

// PrintMatrix()
static void PrintMatrix(CpptrajFile &outfile, const char *Title, double *U, 
                        int mrows, int ncols) 
{
  outfile.Printf("    %s",Title);
  int usize = mrows * ncols;
  for (int i = 0; i < usize; i++) {
    if ( (i%ncols)==0 ) outfile.Printf("\n");
    outfile.Printf(" %10.5lf",U[i]);
  }
  outfile.Printf("\n");
}

// Rotdif::calc_Asymmetric()
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
int Rotdif::calc_Asymmetric(double *Dxyz, double *matrix_D) 
{
  double lambda[8];
  double dx = Dxyz[0];
  double dy = Dxyz[1];
  double dz = Dxyz[2];
  double Dav, Dpr2, delta;
  double rotated_vec[3];

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

  // Loop over all vectors
  double *randvec = random_vectors_;
  for (int i=0; i < nvecs_; i++) {
    // Rotate vector i into D frame
    // This is an inverse rotation. Since matrix_D is in column major order
    // however, do a normal rotation instead.
    // NOTE: Replace this with 3 separate dot products?
    matrix_times_vector(rotated_vec, matrix_D, randvec);
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
    tau1_[i] = ((sintheta2 * (cosphi*cosphi)) / lambda[0]) +
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
    tau2_[i] = (m_m2 / lambda[3]) + (m_m1 / lambda[4]) + (m_0 / lambda[5]) +
               (m_p1 / lambda[6]) + (m_p2 / lambda[7]);
    sumc2_[i] = m_m2 + m_m1 + m_0 + m_p1 + m_p2;

    randvec += 3; // Next random vector
  }

  return 0;
}

// Rotdif::chi_squared()
/** Given Q, calculate effective D values for random_vectors with full 
  * anisotropy and compare to the effective D values in D_eff. Sets
  * D_XYZ with principal components and D_tensor with principal vectors
  * in column major order.
  * \return Deviation of calculated D values from given D values.
  */
double Rotdif::chi_squared(double *Qin) {
#ifdef NO_PTRAJ_ANALYZE
  return -1;
#else
  int n_cols = 3;
  int info;
  double chisq;

  // Convert input Q to Diffusion tensor D
  Q_to_D(D_tensor_, Qin);
  // Diagonalize D; it is assumed workspace (work, lwork) set up prior 
  // to this call.
  // NOTE: Due to the fortran call, the eigenvectors are returned in COLUMN
  //       MAJOR order.
  dsyev_((char*)"Vectors",(char*)"Upper", n_cols, D_tensor_, n_cols, D_XYZ_, work_, lwork_, info);
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
static void calculate_D_properties(double Dxyz[3], double Dout[3]) {
  double Dx = Dxyz[0];
  double Dy = Dxyz[1];
  double Dz = Dxyz[2];
  Dout[0] = (Dx + Dy + Dz) / 3;
  Dout[1] = (2 * Dz) / (Dx + Dy);
  Dout[2] = (1.5 * (Dy - Dx)) / (Dz - (0.5 * (Dx + Dy)) );
}

#define SM_NP 6
#define SM_NP1 7
// Amotry()
double Rotdif::Amotry(double xsmplx[SM_NP1][SM_NP], double *ysearch,
                      double *psum, int ihi, double fac)
{
  double ytry, fac1, fac2, ptry[SM_NP];

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
int Rotdif::Amoeba(double xsmplx[SM_NP1][SM_NP], double *ysearch) {
  int iter;
  bool loop1;
  bool loop2;
  double psum[SM_NP];
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
static void Average_vertices(double *xsearch, double xsmplx[SM_NP1][SM_NP]) {
    for (int j=0; j < SM_NP; j++) {
      xsearch[j] = 0;
      for (int k = 0; k < SM_NP1; k++) 
        xsearch[j] += xsmplx[k][j];
      xsearch[j] /= SM_NP1;
    }
}

// Rotdif::Simplex_min()
/** Main driver routine for Amoeba (downhill simplex) minimizer. In the 
  * simplex method, N+1 initial points (where N is the dimension of the 
  * search space) must be chosen; the SVD solution provides one of these. 
  * Initial points are stored in rows of xsmplx. Components of vector 
  * delq should be of the order of the characteristic "lengthscales" over 
  * which the Q tensor varies. delqfrac determines the size of variation
  * for each of the components of Q; the sign of the variation is randomly 
  * chosen.
  */
int Rotdif::Simplex_min(double *Q_vector) {
  double xsearch[SM_NP];
  double ysearch[SM_NP1];
  double xsmplx[SM_NP1][SM_NP];
  const int nsearch = 1;
  double sgn;
  double d_props[3];
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
    Tau_ = NULL; // Should never get here
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
    calculate_D_properties(D_XYZ_, d_props);

    outfile_.Printf("Input to amoeba - average at cycle %i\n",i+1);
    outfile_.Printf("    Initial chisq = %15.5lf\n",sgn);
    PrintMatrix(outfile_,"Dav, aniostropy, rhombicity:",d_props,1,3);
    PrintMatrix(outfile_,"D tensor eigenvalues:",D_XYZ_,1,3);
    PrintMatrix(outfile_,"D tensor eigenvectors (in columns):",D_tensor_,3,3);

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
    calculate_D_properties(D_XYZ_, d_props);

    outfile_.Printf("Output from amoeba - average at cycle %i\n",i+1);
    outfile_.Printf("    Final chisq = %15.5lf\n",sgn);
    PrintMatrix(outfile_,"Dav, aniostropy, rhombicity:",d_props,1,3);
    PrintMatrix(outfile_,"D tensor eigenvalues:",D_XYZ_,1,3);
    PrintMatrix(outfile_,"D tensor eigenvectors (in columns):",D_tensor_,3,3);
    outfile_.Printf("     taueff(obs) taueff(calc)\n");
    for (int i = 0; i < nvecs_; i++)
      outfile_.Printf("%5i%10.5lf%10.5lf%10.5lf\n",i+1,D_eff_[i],(*Tau_)[i],sumc2_[i]);
    outfile_.Printf("\n");
   
    // cycle over main loop, but first reduce the size of delqfrac:
    delqfrac_ *= 0.750;
    if (debug>0)mprintf("\tAmoeba: Setting delqfrac to %15.7lf\n",delqfrac_);
  }

  // Set q vector to the final average result from simpmin
  for (int n = 0; n < SM_NP; n++) Q_vector[n] = xsearch[n];

  return 0;
}
#undef SM_NP
#undef SM_NP1

// Rotdif::Grid_search()
/** Given a Q tensor, search for a better Q tensor by simple
  * grid search on all 6 elements. If a better solution is found
  * store it.
  * \return 1 if grid search found a better solution.
  * \return 0 if no better solution was found.
  */
int Rotdif::Grid_search(double *Q_vector, int gridsize) {
  double xsearch[6];
  double best[6];
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

// Rotdif::Tensor_Fit()
/** Based on random_vectors and effective diffusion constants D_eff previously
  * calculated, first find the tensor Q (and therefore D) in the small
  * anisotropic limit by solving:
  *   D_eff(n) = At(n) * Q
  * where At(n) is composed of D_eff vector components:
  *   { x^2, y^2, z^2, 2xy, 2yz, 2xz }
  * \param vector_q Will be set with Q tensor for small anisotropic limit.
  */
int Rotdif::Tensor_Fit(double *vector_q) {
#ifdef NO_PTRAJ_ANALYZE
  return 1;
#else
  int info;
  double wkopt;
  double vector_q_local[6];
  double d_props[3];
  double matrix_D_local[9];
  double *deff_local;
  double *At; // Used to index into matrix_At
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
  double *rvecs = random_vectors_;
  At = matrix_At;
  for (int i = 0; i < nvecs_; i++) {
    // Transpose of matrix A
    double *A = matrix_A + i;
    A[ 0] = rvecs[0] * rvecs[0];     // x^2
    A[A1] = rvecs[1] * rvecs[1];     // y^2
    A[A2] = rvecs[2] * rvecs[2];     // z^2
    A[A3] = 2*(rvecs[0] * rvecs[1]); // 2xy
    A[A4] = 2*(rvecs[1] * rvecs[2]); // 2yz
    A[A5] = 2*(rvecs[0] * rvecs[2]); // 2xz
    // Matrix A
    At[0] = A[ 0];
    At[1] = A[A1];
    At[2] = A[A2];
    At[3] = A[A3];
    At[4] = A[A4];
    At[5] = A[A5];
    At += 6;
    rvecs += 3;
  }
  if (debug>1) {
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
  if (debug>0) {
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
  if (debug>1) {
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
  PrintMatrix(outfile_,"Qxx Qyy Qzz Qxy Qyz Qxz",vector_q,1,6);
  // ---------------------------------------------

  // Convert vector Q to diffusion tensor D
  Q_to_D(D_tensor_, vector_q);
  PrintMatrix(outfile_,"D_tensor",D_tensor_,3,3);

  // Save D for later use in back calculating deff from Q
  for (int i = 0; i < 9; i++)
    matrix_D_local[i] = D_tensor_[i];

  // Diagonalize D to find eigenvalues and eigenvectors
  // (principal components and axes)
  // Determine workspace. Do not delete work after set up since the 
  // Rotdif::chi_squared routine will use dsyev_ for diagnolization.
  lwork_ = -1;
  n_cols = 3;
  dsyev_((char*)"Vectors",(char*)"Upper", n_cols, D_tensor_, n_cols, D_XYZ_, &wkopt, lwork_, info);
  lwork_ = (int) wkopt;
  work_ = new double[ lwork_ ];
  // Diagonalize D_tensor
  dsyev_((char*)"Vectors",(char*)"Upper", n_cols, D_tensor_, n_cols, D_XYZ_, work_, lwork_, info);
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
  PrintMatrix(outfile_,"D eigenvalues",D_XYZ_,1,3);
  PrintMatrix(outfile_,"D eigenvectors (in columns)",D_tensor_,3,3);

  // Calculate Dav, Daniso, Drhomb
  calculate_D_properties(D_XYZ_, d_props);
  PrintMatrix(outfile_,"Dav, Daniso, Drhomb",d_props,1,3);

  // Back-calculate the local diffusion constants via At*Q=Deff
  // First convert original D back to Q
  D_to_Q(vector_q_local, matrix_D_local);
  if (debug>0) printMatrix("D_to_Q",vector_q_local,1,6);
  deff_local = new double[ nvecs_ ];
  At = matrix_At;
  // At*Q
  for (int i=0; i < nvecs_; i++) {
    deff_local[i]  = (At[0] * vector_q_local[0]);
    deff_local[i] += (At[1] * vector_q_local[1]);
    deff_local[i] += (At[2] * vector_q_local[2]);
    deff_local[i] += (At[3] * vector_q_local[3]);
    deff_local[i] += (At[4] * vector_q_local[4]);
    deff_local[i] += (At[5] * vector_q_local[5]);
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
  delete[] deff_local; 

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

// Rotdif::DetermineDeffs()
/** Calculate effective diffusion constant for each random vector. 
  * Vectors will be normalized during this phase. First rotate the vector 
  * by all rotation matrices, storing the resulting vectors. The first entry 
  * of rotated_vectors is the original random vector. Then calculate the time 
  * correlation function of the vector. Finally compute the area under the 
  * time correlation function curve and estimate the diffusion constant.
  * Sets D_Eff, normalizes random_vectors.
  */
// TODO: OpenMP Parallelize
int Rotdif::DetermineDeffs() {
  int itotframes;          // Total number of frames (rotation matrices) 
  double *rotated_vectors; // Hold vectors after rotation with Rmatrices
  int maxdat;              // Length of C(t), itotframes + 1 (the original vector)
  double *pX;              // Hold X values of C(t)
  double *p1;              // Hold Y values of C(t) for p1
  double *p2;              // Hold Y values of C(t) for p2
  double *pY;              // Point to p1 or p2 depending on olegendre
  Interpolate spline;      // Used to interpolate C(t)
  // DEBUG
  CpptrajFile outfile;
  // DEBUG

  mprintf("\tDetermining local diffusion constants for each vector.\n");
  ProgressBar progress( nvecs_ );

  itotframes = (int) Rmatrices_.size();
  if (ncorr_ == 0) ncorr_ = itotframes;
  maxdat = ncorr_ + 1;
  // Allocate memory to hold calcd effective D values
  D_eff_ = new double[ nvecs_ ];
  // Allocate memory to hold rotated vectors. Need +1 since the original
  // vector is stored at position 0. 
  rotated_vectors = new double[ 3 * (itotframes+1) ];
  // Allocate memory for C(t)
  p1 = new double[ maxdat ];
  p2 = new double[ maxdat ];
  pX = new double[ maxdat ];
  // Set X values of C(t) based on tfac
  for (int i = 0; i < maxdat; i++) {
    pX[i] = (double) i;
    pX[i] *= tfac_;
  }
  // Allocate mesh to hold interpolated C(t)
  spline.Set_meshX( ti_, tf_, maxdat * 2 );
  // Pointer to random vectors
  double *rndvec = random_vectors_;
  // Check order of legendre polynomial to determine whether we use p1 or p2
  if (olegendre_==1)
    pY = p1;
  else if (olegendre_==2)
    pY = p2;
  else
    return 1; // Should never get here, olegendre is checked in init
  // LOOP OVER RANDOM VECTORS
  for (int vec = 0; vec < nvecs_; vec++) {
    progress.Update( vec );
    // Pointer to soon-to-be rotated vectors
    double *rotvec = rotated_vectors; 
    // Normalize vector
    normalize( rndvec );
    // Assign normalized vector to rotated_vectors position 0
    rotvec[0] = rndvec[0];
    rotvec[1] = rndvec[1];
    rotvec[2] = rndvec[2];
    rotvec += 3;
    // Loop over rotation matrices
    for (std::vector<double*>::iterator rmatrix = Rmatrices_.begin();
                                        rmatrix != Rmatrices_.end();
                                        rmatrix++)
    {
      // Rotate normalized vector
      matrix_times_vector(rotvec, *rmatrix, rndvec);
      // DEBUG
      //mprintf("DBG:%6i%15.8lf%15.8lf%15.8lf\n",ridx,rotated_x[ridx],
      //        rotated_y[ridx],rotated_z[ridx]);
      rotvec += 3;
    }
    // Calculate time correlation function for this vector
    compute_corr(rotated_vectors,maxdat,itotframes,p2,p1);
    // Calculate spline coefficients
    spline.cubicSpline_coeff(pX, pY, maxdat);
    // Calculate mesh Y values
    spline.cubicSpline_eval( pX, pY, maxdat );
    // Integrate
    double integral = spline.Integrate_Trapezoid();
    // Solve for deff
    D_eff_[vec] = calcEffectiveDiffusionConst(integral);

    // DEBUG: Write out p1 and p2 ------------------------------------
    if (!corrOut_.empty() || debug > 3) {
        std::string namebuffer("p1p2.dat");
        if (!corrOut_.empty())
          namebuffer = NumberFilename( corrOut_, vec );
        else
          namebuffer = NumberFilename( namebuffer, vec );
        //NumberFilename(namebuffer, (char*)"p1p2.dat", vec);
        outfile.SetupWrite((char*)namebuffer.c_str(), debug);
        outfile.OpenFile();
        for (int i = 0; i < maxdat; i++) 
          outfile.Printf("%lf %lf %lf\n",pX[i], p2[i], p1[i]);
        outfile.CloseFile();
        //    Write Mesh
        if (debug>3) {
          namebuffer.assign("mesh.dat");
          namebuffer = NumberFilename( namebuffer, vec );
          //NumberFilename(namebuffer, (char*)"mesh.dat", vec);
          outfile.SetupWrite((char*)namebuffer.c_str(), debug);
          outfile.OpenFile();
          for (int i=0; i < spline.Mesh_Size(); i++)
            outfile.Printf("%lf %lf\n", spline.X(i), spline.Y(i));
          outfile.CloseFile();
        }
        if (!corrOut_.empty()) 
          corrOut_.clear();
    }
    if (debug > 0) {
      mprintf("DBG: Vec %i Spline integral= %12.4lf\n",vec,integral);
      mprintf("DBG: deff is %lf\n",D_eff_[vec]);
    }
    // END DEBUG -----------------------------------------------------

    rndvec += 3;
    //break;
  }
  progress.Update( nvecs_ );

  // Cleanup
  delete[] rotated_vectors;
  delete[] p1;
  delete[] p2;
  delete[] pX;
  return 0;
}

// Rotdif::print()
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
void Rotdif::print() {
  double Q_isotropic[6];
  double Q_anisotropic[6];

  mprintf("    ROTDIF:\n");
 
  // Read/Generate nvecs random vectors
  random_vectors_ = randvec();
  if (random_vectors_ == NULL) return;

  // If no rotation matrices generated, exit
  if (Rmatrices_.empty()) return;
  // HACK: To match results from rmscorr.f (where rotation matrices are
  //       implicitly transposed), transpose each rotation matrix.
  // NOTE: Is this actually correct? Want inverse rotation?
  for (std::vector<double*>::iterator rmatrix = Rmatrices_.begin();
                                      rmatrix != Rmatrices_.end();
                                      rmatrix++) {
    matrix_transpose_3x3(*rmatrix);
  }
  // Print rotation matrices
  if (rmOut_!=NULL) {
    CpptrajFile rmout;
    if (rmout.SetupWrite(rmOut_,debug)) {
      mprinterr("    Error: Rotdif: Could not set up %s for writing.\n",rmOut_);
    } else {
      rmout.OpenFile();
      int rmframe=1;
      for (std::vector<double*>::iterator rmatrix = Rmatrices_.begin();
                                          rmatrix != Rmatrices_.end();
                                          rmatrix++) {
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
  if (deffOut_!=NULL) {
    CpptrajFile dout;
    if (dout.SetupWrite(deffOut_,debug)) {
      mprinterr("    Error: Rotdif: Could not set up file %s\n",deffOut_);
    } else {
      dout.OpenFile();
      for (int vec = 0; vec < nvecs_; vec++)
        dout.Printf("%6i%15.8lf\n",vec+1,D_eff_[vec]);
      dout.CloseFile();
    }
  }

/*
  // Create temporary copy of deff for tensorfit_ call
  double *temp_deff = new double[ nvecs_ ];
  for (int i = 0; i < nvecs_; i++)
    temp_deff[i] = D_eff[i];
  // Create temporary value of delqfrac for tensorfit_ call
  double temp_delqfrac = delqfrac;
*/
  // All remaining functions require LAPACK
# ifndef NO_PTRAJ_ANALYZE
  Tensor_Fit( Q_isotropic );

  // Using Q (small anisotropy) as a guess, calculate Q with
  // full anisotropy
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
  //int lflag = olegendre;
  //tensorfit_(random_vectors_,nvecs_,temp_deff,nvecs_,lflag,temp_delqfrac_,debug);
  //delete[] temp_deff;
}
  
