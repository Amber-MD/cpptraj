// Rotdif
#include <cmath>
#include <cstdio> //sscanf
#include "Action_Rotdif.h"
#include "CpptrajStdio.h"
#include "Constants.h" // TWOPI
#include "vectormath.h"
#include "Integrate.h"

// Definition of Fortran subroutines in Rotdif.f called from this class
extern "C" {
  void tensorfit_(double*,int&,double*,int&,int&,double&,int&);
}

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
  randvecIn = NULL;
  rmOut = NULL;
  deffOut = NULL;

  RefParm = NULL;
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
  *                       [rvecin <randvecIn>]
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
  randvecIn = actionArgs.getKeyString("rvecin",NULL);
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

  // Initialize random number generator
  RNgen.rn_set( rseed );

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
  mprintf("            Max iterations = %i, tol = %lf, initial guess = %lf\n",
          itmax, delmin,d0);
  mprintf("            Order of Legendre polynomial = %i\n",olegendre);
  mprintf("            lflag=%i, delqfrac=%.4lf\n",lflag,delqfrac);
  if (randvecIn!=NULL)
    mprintf("            Random vectors will be read from %s\n",randvecIn);
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
  double *XYZ;
  int xyz_size = nvecs * 3;
  CpptrajFile vecIn;
  char buffer[BUFFER_SIZE];

  XYZ = new double[ xyz_size ];

  // ----- Read nvecs vectors from a file
  if (randvecIn!=NULL) {
    if (vecIn.SetupFile(randvecIn, READ, debug)) {
      mprinterr("Error: Could not setup random vectors input file %s",randvecIn);
      delete[] XYZ;
      return NULL;
    }
    if (vecIn.OpenFile()) {
      mprinterr("Error: Could not open random vectors input file %s",randvecIn);
      delete[] XYZ;
      return NULL;
    }
    for (int i = 0; i < xyz_size; i+=3) {
      if (vecIn.IO->Gets(buffer, BUFFER_SIZE) ) {
        mprinterr("Error: Could not read vector %i from file %s\n",(i/3)+1,randvecIn);
        delete[] XYZ;
        return NULL;
      }
      sscanf(buffer,"%*i %lf %lf %lf",XYZ+i,XYZ+i+1,XYZ+i+2);
    }
    vecIn.CloseFile();
  // ----- Generate nvecs vectors
  } else {
    for (int i = 0; i < xyz_size; i+=3) {
      double phi = TWOPI * RNgen.rn_gen();
      // XYZ[i+2] is cos(theta)
      XYZ[i+2] = 1 - RNgen.rn_gen();
      double theta = acos( XYZ[i+2] );
      double sintheta = sin( theta );
      XYZ[i  ] = sintheta * cos( phi );
      XYZ[i+1] = sintheta * sin( phi );
    }
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

// Rotdif::compute_corr()
/** Given a vector that has been randomly rotated itotframes times, compute
  * the time correlation function of the vector.
  * \param ncorr maximum length to compute time correlation functions - units of 'frames'
  * \param itotframes total number of frames provided
  * \param rotated_vectors array of vector coords for each frame, V0x,V0y,V0z,V1x,V1y,V1z 
  */
// rotated_vectors should be size itotframes+1
// p2 and p1 should be size ncorr+1
int Rotdif::compute_corr(double *rotated_vectors, int maxdat, int itotframes, 
                         double *p2, double *p1)
{
  double *VJ, *VK;
  // Initialize p1 and p2
  for (int i = 0; i < maxdat; i++) {
    p1[i] = 0.0;
    p2[i] = 0.0;
  }
  // Normalize all vectors
  //VJ = rotated_vectors;
  //for (int i = 0; i <= itotframes; i++) {
  //  normalize( VJ );
  //  VJ += 3;
  //}

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
/** computes effect diffusion constant for a vector using its correlation 
  * function as input. Starting with definition:
  *
  *   6*D=integral[0,inf;C(t)] 
  *
  * integrate C(t) from ti -> tf yielding F(ti,tf).
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
  */
// itmax:  maximum number of iterations in subroutine itsolv
// delmin:  convergence criterion used in subroutine itsolv;
//          maximum accepted fractional change in successive 
//          iterations
// d0: initial guess for diffusion constant; accurate estimate not
//     needed
// l: order of Legendre polynomial in the correlation function
//    <P(l)>
double Rotdif::calcEffectiveDiffusionConst(double f ) {
//       ti,tf:  Integration limits.
//    subroutine itsolv(itmax,delmin,l,d0,ti,tf,f,d,info)
//    solves the equation 6*D=[exp(-6*D*ti)-exp(-6*D*tf)]/F(ti,tf) iteratively, 
//    by putting 6*D(i+1)=[exp(-6*D(i)*ti)-exp(-6*D(i)*tf)]/F(ti,tf)
//    where F(ti,tf) is input (integral[dt*C(t)] from ti->tf)
  double d; 
  double del;
  int i;
  double fac;

  double l = (double) olegendre;
  fac = (l*(l+1));
  i=1;
  del=10000000000;
  while ( i<=itmax && del>delmin) {
     d = ( exp(-fac*d0*ti) - exp(-fac*d0*tf) );
     d = d / (fac*f);
     del = (d-d0)/d0;
     if (del < 0) del = -del;
     //del = abs( (d-d0)/d0 );
     //if (debug>0)
       mprintf("ITSOLV: %6i  %15.8e  %15.8e  %15.8e\n", i,d0,d,del);
     d0 = d;
     ++i;
  }
  if ( i>itmax && del>delmin) {
     mprintf("\tWarning, itsolv did not converge: # iterations=%i, fractional change=%lf\n",
             i, del);
  } else {
    mprintf("\tConverged: # iterations=%i\n",i);
  }

  return d; 
}

// Rotdif::print()
void Rotdif::print() {
  double *random_vectors;    // Hold nvecs random vectors
  double *deff;              // Hold calculated effective D
  double *rotated_vectors;   // Hold vectors after rotation with Rmatrices
  double *pX;                // Hold X values of C(t)
  double *p1;                // Hold Y values of C(t) for p1
  double *p2;                // Hold Y values of C(t) for p2
  double deff_val;           
  Interpolate spline;        // Used to interpolate C(t)
  int itotframes;            // Total number of frames (rotation matrices) 
  int maxdat;                // Length of C(t), itotframes + 1 (the original vector)
  // DEBUG
  CpptrajFile outfile;
  char namebuffer[32];
  // DEBUG
 
  // If no rotation matrices generated, exit
  if (Rmatrices.empty()) return;
 
  // Generate nvecs random vectors
  random_vectors = randvec();
  if (random_vectors == NULL) return;

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

  // Calculate effective diffusion constant deff for each random vector. 
  // Vectors will be normalized during this phase. First rotate the vector 
  // by all rotation matrices, storing the resulting vectors. The first entry 
  // of rotated_vectors is the original random vector. Then calculate the time 
  // correlation function of the vector. Finally compute the area under the 
  // time correlation function curve and estimate the diffusion constant.
  itotframes = (int) Rmatrices.size();
  if (ncorr == 0) ncorr = itotframes;
  maxdat = ncorr + 1;
  // Allocate memory to hold vectors, Deff, and C(t)
  rotated_vectors = new double[ 3 * maxdat ];
  deff = new double[ nvecs ];
  p1 = new double[ maxdat ];
  p2 = new double[ maxdat ];
  pX = new double[ maxdat ];
  // Set X values of C(t) based on tfac
  for (int i = 0; i < maxdat; i++) {
    pX[i] = (double) i;
    pX[i] *= tfac;
  }
  // Allocate mesh to hold interpolated C(t)
  int mesh_size = maxdat * 2;
  double *mesh_x = new double[ mesh_size ];
  double *mesh_y = new double[ mesh_size ];
  set_xvalues_range(mesh_x, ti, tf, mesh_size); 
  // Pointer to random vectors
  double *rndvec = random_vectors; 
  // LOOP OVER RANDOM VECTORS
  for (int vec = 0; vec < nvecs; vec++) {
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
    for (std::vector<double*>::iterator rmatrix = Rmatrices.begin();
                                        rmatrix != Rmatrices.end();
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
    //dlocint_(ti,tf,itmax,delmin,d0,olegendre,deff_val);
    // Calculate spline coefficients
    spline.cubicSpline_coeff(pX, p2, maxdat);
    // Calculate mesh Y values
    spline.cubicSpline_eval(mesh_x,mesh_y,mesh_size, pX, p2, maxdat);
    // Integrate
    double integral = integrate_trapezoid(mesh_x, mesh_y, mesh_size);
    // Solve for deff
    deff_val = calcEffectiveDiffusionConst(integral);
    deff[vec] = deff_val;

    // DEBUG: Write out p1 and p2 ------------------------------------
    if (debug > 1) {
      NumberFilename(namebuffer, (char*)"p1p2.dat", vec);
      outfile.SetupFile(namebuffer,WRITE,debug);
      outfile.OpenFile();
      for (int i = 0; i < maxdat; i++) 
        outfile.IO->Printf("%lf %lf %lf\n",pX[i], p2[i], p1[i]);
      outfile.CloseFile();
      //    Write Mesh
      NumberFilename(namebuffer, (char*)"mesh.dat", vec);
      outfile.SetupFile(namebuffer, WRITE, debug);
      outfile.OpenFile();
      for (int i=0; i < mesh_size; i++)
        outfile.IO->Printf("%lf %lf\n",mesh_x[i],mesh_y[i]);
      outfile.CloseFile(); 
      mprintf("DBG: Vec %i Spline integral= %12.4lf\n",vec,integral);
      mprintf("DBG: deff is %lf\n",deff_val);
    }
    // END DEBUG -----------------------------------------------------

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

  //tensorfit_(random_vectors,nvecs,deff,nvecs,lflag,delqfrac,debug);

  // Cleanup
  delete[] random_vectors;
  delete[] deff;
  delete[] rotated_vectors;
  delete[] p1;
  delete[] p2;
  delete[] pX;
  delete[] mesh_x;
  delete[] mesh_y;
}
  
