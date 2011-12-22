// Rotdif
#include <cmath>
#include "Action_Rotdif.h"
#include "CpptrajStdio.h"
#include "Constants.h" // TWOPI
#include "vectormath.h"

// Definition of Fortran subroutines in Rotdif.f called from this class
extern "C" {
  void compute_corr_(double*, double*, double*, int&, int&, int&, double*, double*);
  void dlocint_(double&,double&,int&,double&,double&,int&,double&);
  double random_(int&);
  void tensorfit_(double*,int&,double*,int&,int&,double&,int&);
}

// Fortran common block used by rmscorr.f functions
#define MAXDAT 100000
typedef struct {
  //double *tdat;
  //double *p2;
  double tdat[MAXDAT];
  double p2[MAXDAT];
  int ndat;
} rmscorr_common;

// CONSTRUCTOR
Rotdif::Rotdif() {
  //fprintf(stderr,"Rotdif Con\n");
  rseed = 1;
  nvecs = 0;
  tfac = 0.0;
  ti = 0.0;
  tf = 0.0;
  itmax = 0;
  delmin = 0.0;
  d0 = 0.0;
  olegendre = 0;
  ncorr = 0;
  lflag = 0;
  delqfrac = 0;

  randvecOut = NULL;
  rmOut = NULL;
  deffOut = NULL;

  RefParm = NULL;

  //iff = 0;
  //inext = 0;
  //inextp = 0;
  //for (int i = 0; i < 55; i++)
  //  ma[i] = 0;
} 

// DESTRUCTOR
Rotdif::~Rotdif() {
  //fprintf(stderr,"Rotdif Destructor.\n");
  for (std::vector<double*>::iterator Rmatrix = Rmatrices.begin();
                                      Rmatrix != Rmatrices.end();
                                      Rmatrix++)
    delete[] *Rmatrix;
}

// Rotdif::init()
/** Expected call: rotdif [rseed <rseed>] [nvecs <nvecs>]  
  *                       ref <refname> | refindex <refindex> | reference
  *                       [<refmask>] [ncorr <ncorr>] dt <tfac> [ti <ti>] tf <tf>
  *                       [itmax <itmax>] [tol <delmin>] [d0 <d0>] [order <olegendre>]
  *                       [lflag <lflag>] [delqfrac <delqfrac>] [rvecout <randvecOut>]
  *                       [rmout <rmOut>] [deffout <deffOut>]
  */ 
// Dataset name will be the last arg checked for. Check order is:
//    1) Keywords
//    2) Masks
//    3) Dataset name
int Rotdif::init( ) {
  int refindex;
  char *referenceName, *mask0;
  // Get Keywords
  nvecs = actionArgs.getKeyInt("nvecs",1000);
  rseed = actionArgs.getKeyInt("rseed",80531);
  ncorr = actionArgs.getKeyInt("ncorr",0);
  tfac = actionArgs.getKeyDouble("dt",0);
  if (tfac<=0) {
    mprinterr("    Error: Rotdif: dt (timestep) must be specified and > 0.\n");
    return 1;
  }
  ti = actionArgs.getKeyDouble("ti",0);
  tf = actionArgs.getKeyDouble("tf",0);
  if (tf <= ti) {
    mprinterr("    Error: Rotdif: Initial time ti (%lf) must be < final time tf (%lf).\n",
              ti, tf);
    return 1;
  }
  itmax = actionArgs.getKeyInt("itmax",500);
  delmin = actionArgs.getKeyDouble("tol",0.000001);
  d0 = actionArgs.getKeyDouble("d0",0.03);
  olegendre = actionArgs.getKeyInt("order",2);
  lflag = actionArgs.getKeyInt("lflag",2);
  delqfrac = actionArgs.getKeyDouble("delqfrac",0.5);
  randvecOut = actionArgs.getKeyString("rvecout",NULL);
  rmOut = actionArgs.getKeyString("rmout",NULL);
  deffOut = actionArgs.getKeyString("deffout",NULL);

  referenceName=actionArgs.getKeyString("ref",NULL);
  refindex=actionArgs.getKeyInt("refindex",-1);
  if (actionArgs.hasKey("reference")) // For compatibility with ptraj
    refindex = 0;

  // Get Masks
  mask0 = actionArgs.getNextMask();
  RefMask.SetMaskString(mask0);
  TargetMask.SetMaskString(mask0);

  // Dataset

  // Set up reference for RMSD
  // Attempt to get reference index by name
  if (referenceName!=NULL)
    refindex=FL->GetFrameIndex(referenceName);
  // Get reference frame by index
  Frame *TempFrame=FL->GetFrame(refindex);
  if (TempFrame==NULL) {
    mprinterr("    Error: Rotdif::init: Could not get reference index %i\n",refindex);
    return 1;
  }
  RefFrame = *TempFrame;
  // Set reference parm
  RefParm=FL->GetFrameParm(refindex);
  // Setup reference mask
  if (RefParm->SetupIntegerMask( RefMask, activeReference )) return 1;
  if (RefMask.None()) {
    mprintf("    Error: Rotdif::init: No atoms in reference mask.\n");
    return 1;
  }
  // Allocate frame for selected reference atoms
  SelectedRef.SetupFrameFromMask(&RefMask, RefParm->mass); 

  mprintf("    ROTDIF: Random seed %i, # of random vectors to generate: %i\n",rseed,nvecs);
  mprintf("            Max length to compute time correlation functions:");
  if (ncorr == 0)
    mprintf(" Total # of frames\n");
  else
    mprintf(" %i\n",ncorr);
  mprintf("            Timestep = %.4lf, T0 = %.4lf, TF = %.4lf\n",tfac,ti,tf);
  mprintf("            Max iterations = %i, tol = %lf, initial guess = %lf\n",itmax,
          delmin,d0);
  mprintf("            Order of Legendre polynomial = %i\n",olegendre);
  mprintf("            lflag=%i, delqfrac=%.4lf\n",lflag,delqfrac);
  if (randvecOut!=NULL)
    mprintf("            Random vectors will be written to %s\n",randvecOut);
  if (rmOut!=NULL)
    mprintf("            Rotation matrices will be written out to %s\n",rmOut);
  if (deffOut!=NULL)
    mprintf("            Deff will be written out to %s\n",deffOut);

  return 0;
}

// Rotdif::setup()
/** Determine what atoms each mask pertains to for the current parm file.
  * Also determine whether imaging should be performed.
  */
int Rotdif::setup() {

  if ( currentParm->SetupIntegerMask( TargetMask, activeReference ) ) return 1;
  if ( TargetMask.None() ) {
    mprintf("    Error: Rotdif::setup: No atoms in mask.\n");
    return 1;
  }
  // Allocate space for selected atoms in the frame. This will also put the
  // correct masses in based on the mask.
  SelectedTarget.SetupFrameFromMask(&TargetMask, currentParm->mass);
  // Check that num atoms in frame mask from this parm match ref parm mask
  if ( RefMask.Nselected != TargetMask.Nselected ) {
    mprintf( "    Error: Number of atoms in RMS mask (%i) does not \n",TargetMask.Nselected);
    mprintf( "           equal number of atoms in Ref mask (%i).\n",RefMask.Nselected);
    return 1;
  }
  
  // Print info for this parm
  mprintf("    ROTDIF:\n");
        
  return 0;  
}

// Rotdif::action()
/** Calculate and store the rotation matrix for frame to reference.
  */
int Rotdif::action() {
  double R, Trans[6], *U;

  // Set selected reference atoms - always done since RMS fit modifies SelectedRef 
  SelectedRef.SetFrameCoordsFromMask(RefFrame.X, &RefMask);

  // Set selected frame atoms. Masses have already been set.
  SelectedTarget.SetFrameCoordsFromMask(currentFrame->X, &TargetMask);

  U = new double[ 9 ];
  R = SelectedTarget.RMSD(&SelectedRef, U, Trans, useMass);

  Rmatrices.push_back( U );

  return 0;
} 

// Rotdif::random_number()
/** Generate random number based on rseed.
  */
/*double Rotdif::random_number() {
  float mbig = 4000000;
  float mseed = 1618033;
  float mz = 0;
  float fac = 1 / mbig;
  float mj, mk, random, seed;
  double random_ret;
 
  seed = (float) rseed;
  if (seed < 0 || iff < 0) {
    iff = 1;
    mj = mseed - fabs(seed);
    mj = fmod(mj, mbig);
    for (int i=0; i < 55; i++)
      ma[i] = mj;
    mk = 1;
    for (int i=0; i < 54; i++) {
      int ii = (21*i) % 55;
      ma[ii] = mk;
      mk = mj - mk;
      if ( mk < mz) mk = mk + mbig;
      mj = ma[ii];
    } 
    for (int k = 0; k < 4; k++) { 
      for (int i = 0; i < 55; i++) {
        int idx = (i+30) % 55;
        ma[i] = ma[i] - ma[idx];
        if (ma[i] < mz) ma[i] = ma[i] + mbig;
      }
    }
    inext = 0;
    inextp = 31;
    seed = fabs(seed);
    rseed = (int) seed;
  } 
  ++inext;
  if (inext==55) inext = 0;
  ++inextp;
  if (inextp==55) inextp = 0;
  mj = ma[inext] - ma[inextp];
  if (mj < mz) mj = mj + mbig;
  ma[inext] = mj;
  random = mj * fac;
  random_ret = (double) random;
  return random_ret; 
}*/

// Rotdif::randvec()
/** Generate a set of points uniformly covering surface of unit sphere.
  *   u,v: uniformly distributed random deviates 0 < u, v < 1
  *   azimuthal angle phi: rho(phi)=1/2*pi -> phi=2*pi*u 0 < phi < 2*pi
  *   polar angle theta: rho(theta)=sin(theta) 
  *                      -> theta=arccos(1-v) 0 < theta < 0.5*pi
  *   only a single hemisphere is required
  *   rho(theta)=0.5d0*sin(theta) -> theta=arccos(1-2*v) 0 < theta < pi
  *   if both hemispheres are required
  *   see Numerical Recipes sec. 7.2
  *   http://mathworld.wolfram.com/SpherePoi.html
  */
double *Rotdif::randvec() {
  double *XYZ;
  int xyz_size = nvecs * 3;

  XYZ = new double[ xyz_size ];
  for (int i = 0; i < xyz_size; i+=3) {
    //phi=2d0*pi*random(seed)
    //double phi = TWOPI * random_number();
    double phi = TWOPI * random_(rseed);
    //theta=dacos(1d0-random(seed))
    //double theta = acos( 1 - random_number() );
    double theta = acos( 1 - random_(rseed) );
    double sintheta = sin( theta );
    //x=dsin(theta)*dcos(phi)
    XYZ[i  ] = sintheta * cos( phi );
    //y=dsin(theta)*dsin(phi)
    XYZ[i+1] = sintheta * sin( phi );
    //z=dcos(theta)
    XYZ[i+2] = cos( theta );
  }
  // Print vectors
  if (randvecOut!=NULL) {
    CpptrajFile rvout;
    if (rvout.SetupFile(randvecOut,WRITE,debug)) {
      mprinterr("    Error: Rotdif: Could not set up %s for writing.\n",randvecOut);
    } else {
      rvout.OpenFile();
      int idx = 0;
      for (int i = 1; i <= nvecs; i++) {
        rvout.IO->Printf("%6i  %15.8lf  %15.8lf  %15.8lf\n",i,XYZ[idx],XYZ[idx+1],XYZ[idx+2]);
        idx += 3;
      }
      rvout.CloseFile();
    }
  }

  return XYZ;
}
   
// Rotdif::print()
void Rotdif::print() {
  double *random_vectors;
  double vec_result[3];
  double *deff;
  double *rotated_x;
  double *rotated_y;
  double *rotated_z;
  double *p1;
  //double *p2;
  double deff_val;
  extern rmscorr_common dat_; // For common block in rmscorr
  int itotframes = (int) Rmatrices.size();
  int maxdat = itotframes + 1;
 
  // If no rotation matrices generated, exit
  if (Rmatrices.empty()) return;
 
  // Generate nvecs random vectors
  random_vectors = randvec();
  if (random_vectors == NULL) return;

  if (ncorr == 0) ncorr = itotframes;

  // Allocate space for rmscorr common block
  //dat_.tdat = new double[ maxdat ];
  //dat_.p2 = new double[ maxdat ];
  dat_.ndat = maxdat;
  for (int i = 0; i < maxdat; i++) {
    dat_.tdat[i] = (double) i;
    dat_.tdat[i] *= tfac;
  }

  // HACK: To match results from rmscorr.f (where rotation matrices are
  //       implicitly transposed), transpose each rotation matrix.
  for (std::vector<double*>::iterator rmatrix = Rmatrices.begin();
                                      rmatrix != Rmatrices.end();
                                      rmatrix++) {
    matrix_transpose(*rmatrix);
  }

  // Print rotation matrices
  if (rmOut!=NULL) {
    CpptrajFile rmout;
    if (rmout.SetupFile(rmOut,WRITE,debug)) {
      mprinterr("    Error: Rotdif: Could not set up %s for writing.\n",rmOut);
    } else {
      rmout.OpenFile();
      int rmframe=1;
      for (std::vector<double*>::iterator rmatrix = Rmatrices.begin();
                                          rmatrix != Rmatrices.end();
                                          rmatrix++) {
        rmout.IO->Printf("%13i %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f\n",
             rmframe++,
            (*rmatrix)[0], (*rmatrix)[1], (*rmatrix)[2],
            (*rmatrix)[3], (*rmatrix)[4], (*rmatrix)[5],
            (*rmatrix)[6], (*rmatrix)[7], (*rmatrix)[8]);
      }
      rmout.CloseFile();
    }
  }

  // For each random vector, rotate by all rotation matrices, storing the
  // resulting vectors. The first entry of rotated_ is the original
  // random vector.
  //double *rotated_vectors = new double[ 3 * itotframes ];
  deff = new double[ nvecs ];
  rotated_x = new double[ maxdat ];
  rotated_y = new double[ maxdat ];
  rotated_z = new double[ maxdat ];
  p1 = new double[ maxdat ];
  //p2 = new double[ maxdat ];
  double *rndvec = random_vectors; // Pointer to random vectors
  //double *vec_result = rotated_vectors; // Pointer to rotated vectors
  for (int vec = 0; vec < nvecs; vec++) {
    rotated_x[0] = rndvec[0];
    rotated_y[0] = rndvec[1];
    rotated_z[0] = rndvec[2];
    int ridx = 1;
    for (std::vector<double*>::iterator rmatrix = Rmatrices.begin();
                                        rmatrix != Rmatrices.end();
                                        rmatrix++)
    {
      matrix_times_vector(vec_result, *rmatrix, rndvec);
      rotated_x[ridx] = vec_result[0];
      rotated_y[ridx] = vec_result[1];
      rotated_z[ridx] = vec_result[2];
  
      // DEBUG
      //mprintf("DBG:%6i%15.8lf%15.8lf%15.8lf\n",ridx,rotated_x[ridx],
      //        rotated_y[ridx],rotated_z[ridx]);

      ++ridx;
      //vec_result += 3;
    }
    compute_corr_(rotated_x,rotated_y,rotated_z,ncorr,itotframes,maxdat,dat_.p2,p1);
    //for (int i = 0; i < maxdat; i++)
    //  mprintf("%15.8lf%15.8lf\n",dat_.tdat[i],dat_.p2[i]);
    //break;

    dlocint_(ti,tf,itmax,delmin,d0,olegendre,deff_val);

    deff[vec] = deff_val;

    rndvec += 3;
    //break;
  }

  // Print deff
  if (deffOut!=NULL) {
    CpptrajFile dout;
    if (dout.SetupFile(deffOut,WRITE,debug)) {
      mprinterr("    Error: Rotdif: Could not set up file %s\n",deffOut);
    } else {
      dout.OpenFile();
      for (int vec = 0; vec < nvecs; vec++)
        dout.IO->Printf("%6i%15.8lf\n",vec+1,deff[vec]);
      dout.CloseFile();
    }
  }

  tensorfit_(random_vectors,nvecs,deff,nvecs,lflag,delqfrac,debug);

  // Cleanup
  delete[] random_vectors;
  delete[] deff;
  delete[] rotated_x;
  delete[] rotated_y;
  delete[] rotated_z;
  delete[] p1;
  //delete[] dat_.tdat;
  //delete[] dat_.p2;
}
  
