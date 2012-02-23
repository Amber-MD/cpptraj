/*! \file vectormath.cpp
 *
 * Simple vector/matrix math routines.
 * ROTATE, jacobi3, and diagEsort routines adapted from PTRAJ
 */
#include <cmath>
#include "vectormath.h"
#include "CpptrajStdio.h"
#include "Constants.h" // RADDEG, PI

// normalize()
/** Normalize vector in a[].  */
void normalize(double a[3]) {
  double r2;
  double b;

  r2 = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
  b = 1.0 / sqrt(r2);

  a[0] *= b;
  a[1] *= b;
  a[2] *= b;
}

// vector_norm()
/** Normalize vector in V.
  * \param V Vector to be normalized
  * \param r2 will be set to the length^2
  * \return vector length
  */
double vector_norm(double V[3], double *r2) {
  double b, r;

  *r2=(V[0]*V[0])+(V[1]*V[1])+(V[2]*V[2]);
  r = sqrt(*r2);
  b = 1.0 / r;
  V[0] *= b;
  V[1] *= b;
  V[2] *= b;
  return r;
}

// vector_sub()
/** V = U - W */
void vector_sub(double V[3], double U[3], double W[3]) {
  V[0] = U[0] - W[0];
  V[1] = U[1] - W[1];
  V[2] = U[2] - W[2];
}

// vector_sum()
/** V = U + W */
void vector_sum(double V[3], double U[3], double W[3]) {
  V[0] = U[0] + W[0];
  V[1] = U[1] + W[1];
  V[2] = U[2] + W[2];
}

// dot_product()
/** total = V . U */
double dot_product(double V[3], double U[3]) {
  double total;

  total = (V[0] * U[0]) + (V[1] * U[1]) + (V[2] * U[2]);

  return total;
}

// dot_product_angle()
/** Return the angle obtained from the dot product between vectors V 
  * and U. Only works correctly if V and U are normalized beforehand.
  */
double dot_product_angle(double V[3], double U[3]) {
  double total;

  total = V[0] * U[0];
  total += (V[1] * U[1]);
  total += (V[2] * U[2]);

  return acos(total);
}

// cross_product()
void cross_product(double V[3], double U[3], double W[3]) {
  V[0] = U[1]*W[2] - U[2]*W[1]; // UyWz - UzWy
  V[1] = U[2]*W[0] - U[0]*W[2]; // UzWx - UxWz
  V[2] = U[0]*W[1] - U[1]*W[0]; // UxWy - UyWx
}

// dot_product_sign()
/** Return the angle obtained from the dot product between vectors V
  * and U, with sign determined from (VxU) dot Z. Assumes V and U
  * are normalized.
  */
double dot_product_sign(double *V, double *U, double *Z) {
  double Vec[3];
  double dp = dot_product_angle(V, U);
  cross_product(Vec, V, U);
  double sign = dot_product(Vec, Z);
  if (sign < 0) dp = -dp;
  return dp;
}

// matrix_transpose()
/** M = Ut
  * Columns of U become rows of M and vice versa.
  */
void matrix_transpose(double M[9], double U[9]) {
  M[0] = U[0];
  M[1] = U[3];
  M[2] = U[6];
  M[3] = U[1];
  M[4] = U[4];
  M[5] = U[7];
  M[6] = U[2];
  M[7] = U[5];
  M[8] = U[8];
}

// matrix_transpose_3x3()
/** M = Mt
  * Columns of M become rows of M and vice versa.
  */
void matrix_transpose_3x3(double M[9]) {
  //double U0,U4,U8;
  double U1,U2,U3,U5,U6,U7;
  //U0 = M[0];
  U1 = M[1];
  U2 = M[2];

  U3 = M[3];
  //U4 = M[4];
  U5 = M[5];

  U6 = M[6];
  U7 = M[7];
  //U8 = M[8];

  //M[0] = U0;
  M[1] = U3;
  M[2] = U6;

  M[3] = U1;
  //M[4] = U4;
  M[5] = U7;

  M[6] = U2;
  M[7] = U5;
  //M[8] = U8;
}

// matrix_transpose()
double *matrix_transpose(double *M, int mrows, int ncols) {
  double *result = new double[ mrows * ncols];
  int midx = 0;
  for (int m = 0; m < mrows; m++) {
    int ridx = m;
    for (int n = 0; n < ncols; n++) {
      //mprintf("TRANSPOSE: %6i = %6i\n",ridx,midx);
      result[ridx] = M[midx];
      ++midx;
      ridx += mrows;
    }
  }
  return result;
}
  
// matrix_times_vector()
/** Multiply matrix R by vector V, store result in U
  */
void matrix_times_vector(double U[3], double R[9], double V[3]) {
  double x,y,z;
  // Store V so that U and V may overlap.
  x = V[0];
  y = V[1];
  z = V[2];
  U[0] = (R[0]*x) + (R[1]*y) + (R[2]*z);
  U[1] = (R[3]*x) + (R[4]*y) + (R[5]*z);
  U[2] = (R[6]*x) + (R[7]*y) + (R[8]*z);
}

// matrixT_times_vector()
/** Multiply transpose of matrix R by vector V, store result in U
  */
void matrixT_times_vector(double U[3], double R[9], double V[3]) {
  double x,y,z;
  // Store V so that U and V may overlap.
  x = V[0];
  y = V[1];
  z = V[2];
  U[0] = (R[0]*x) + (R[3]*y) + (R[6]*z);
  U[1] = (R[1]*x) + (R[4]*y) + (R[7]*z);
  U[2] = (R[2]*x) + (R[5]*y) + (R[8]*z);
}

// matrix_multiply_3x3()
/** Multiply 3x3 matrix Row by 3x3 matrix Col, store result in M
  */
void matrix_multiply_3x3(double M[9], double Row[9], double Col[9]) {
  M[0] = (Row[0] * Col[0]) + (Row[1] * Col[3]) + (Row[2] * Col[6]);
  M[1] = (Row[0] * Col[1]) + (Row[1] * Col[4]) + (Row[2] * Col[7]);
  M[2] = (Row[0] * Col[2]) + (Row[1] * Col[5]) + (Row[2] * Col[8]);
  M[3] = (Row[3] * Col[0]) + (Row[4] * Col[3]) + (Row[5] * Col[6]);
  M[4] = (Row[3] * Col[1]) + (Row[4] * Col[4]) + (Row[5] * Col[7]);
  M[5] = (Row[3] * Col[2]) + (Row[4] * Col[5]) + (Row[5] * Col[8]);
  M[6] = (Row[6] * Col[0]) + (Row[7] * Col[3]) + (Row[8] * Col[6]);
  M[7] = (Row[6] * Col[1]) + (Row[7] * Col[4]) + (Row[8] * Col[7]);
  M[8] = (Row[6] * Col[2]) + (Row[7] * Col[5]) + (Row[8] * Col[8]);
}

// matrix_multiply()
/** Multiply matrix M by matrix N.
  * \return matrix of size mrow*ncol containing MxN
  */
/*double *matrix_multiply(double *M, int mrow, int mcol,
                        double *N, int nrow, int ncol)
{
  int result_size = mrow * ncol;
  double *result = new double[ result_size ];

  for (int ridx = 0; ridx < result_size; ridx++) {
    double sum = 0;
    for (int */
  
 

// matrix_to_angle()
/** Return angle of rotation from rotation matrix according to
  * cos(t)=(trace(R)-1)/2
  * Equation taken from :
  *   3D game engine design: a practical approach to real-time Computer Graphics,
  *   Volume 385, By David H. Eberly, 2001, p. 16.
  */
double matrix_to_angle(double U[9]) {
  double trace;

  trace = U[0] + U[4] + U[8];

  trace = (trace - 1) / 2;

  return acos(trace);
}

// axis_of_rotation()
/** If theta is between 0 and pi extract axis of rotation from rotation matrix
  * U according to:
  *   R - Rt = (2 * sin(theta)) * S, where S is:
  *     0 -z  y
  *     z  0 -x
  *    -y  x  0
  * Place result in V.
  */
int axis_of_rotation(double V[3], double U[9], double theta) {
  double dx;
  if (theta>0 && theta<PI) {
    dx = 1 / (2 * sin(theta));
    V[0]=(U[5]-U[7]) * dx;
    V[1]=(U[6]-U[2]) * dx;
    V[2]=(U[1]-U[3]) * dx;
    normalize(V);
    return 0;
  } else {
    mprintf("Error: axis_of_rotation: Could not extract axis of rotation, angle is %lf\n",
            RADDEG*theta);
  }
  return 1;
}

// calcRotationMatrix()
/** Given an axis of rotation V and a magnitude (radians), calculate a 
  * rotation matrix and store it in T
  */
void calcRotationMatrix(double T[9], double V[3], double theta) {
  double ux2,uxuy,uxuz,uy2,uyuz,uz2,c,s,c1,uxs,uys,uzs;

  // Compute all prefactors
  ux2 =V[0]*V[0];
  uxuy=V[0]*V[1];
  uxuz=V[0]*V[2];
  uy2 =V[1]*V[1];
  uyuz=V[1]*V[2];
  uz2 =V[2]*V[2];
  c=cos(theta);
  s=sin(theta);
  c1=1-c;
  uxs=V[0]*s;
  uys=V[1]*s;
  uzs=V[2]*s;

  // Store rotation matrix elements
  T[0]=c + (ux2 * c1); 
  T[3]=(uxuy * c1) + uzs;
  T[6]=(uxuz * c1) - uys;

  T[1]=(uxuy * c1) - uzs;
  T[4]=c + (uy2 * c1); 
  T[7]=(uyuz * c1) + uxs;

  T[2]=(uxuz * c1) + uys;
  T[5]=(uyuz * c1) - uxs;
  T[8]=c + (uz2 * c1); 
}

// calcRotationMatrix()
/** Given rotations around the X, Y, and Z axes (radians), calculate a
  * rotation matrix and store it in T.
  */
void calcRotationMatrix(double T[9], double psiX, double psiY, double psiZ) {
  double Psi;
  double V[3];

  Psi = sqrt( (psiX*psiX) + (psiY*psiY) + (psiZ*psiZ) );
  //mprintf("\t\tcalcRotationMatrix(%.2lf,%.2lf,%.2lf) Psi=%lf\n",
  //        psiX*RADDEG,psiY*RADDEG,psiZ*RADDEG,Psi*RADDEG);
  V[0] = psiX / Psi;
  V[1] = psiY / Psi;
  V[2] = psiZ / Psi;

  calcRotationMatrix(T, V, Psi);
}

// ROTATE()
#define ROTATE(ARR,MAJ1,MIN1,MAJ2,MIN2) { \
  g = ARR[MAJ1 + MIN1]; \
  h = ARR[MAJ2 + MIN2]; \
  ARR[MAJ1 + MIN1] = g - s*(h+g*tau); \
  ARR[MAJ2 + MIN2] = h + s*(g-h*tau); }

// jacobi3()
/** Diagonalize 3x3 matrix with jacobi method.
  */
static int jacobi3(double *a, double *d, double *v, int *nrot) { 
/* n must be 3.  see b[3] and z[3] below */
  int  i, j, ip, iq, p3, j3;
  double  tresh, theta, tau, t, sm, s, h, g, c, b[3], z[3];

  for (ip=p3=0; ip<3; ip++,p3+=3) {
    /* initialize the identity matrix */
    for (iq=0; iq<3; iq++)
      v[p3 + iq] = 0.0;
    v[p3 + ip] = 1.0;
    /* initialize b and d to diagonal of a */
    b[ip] = d[ip] = a[p3 + ip];
    z[ip] = 0.0;
  }
  *nrot = 0;
  for (i=0; i<50; i++) {    /* 50 tries */

    sm = 0.0;
    for (ip=p3=0; ip<2; ip++,p3+=3) {
      for (iq=ip+1; iq<3; iq++)
        sm += fabs(a[p3 + iq]);
    }

    if (sm == 0.0) {
      return(1);
    }
    if (i < 3)
      tresh = sm * 0.2 / 9.0;   /* on 1st three sweeps... */
    else
      tresh = 0.0;      /* thereafter... */
    for (ip=p3=0; ip<2; ip++,p3+=3) {
      for (iq=ip+1; iq<3; iq++) {
        g = 100.0 * fabs(a[p3 + iq]);

        if ( i > 3  &&  fabs(d[ip])+g == fabs(d[ip])
        && fabs(d[iq])+g == fabs(d[iq])) {
          a[p3 + iq] = 0.0;
        } else if (fabs(a[p3 + iq]) > tresh) {
          h = d[iq]-d[ip];
          if (fabs(h)+g==fabs(h))
            t = a[p3 + iq] / h;
          else {
            theta = 0.5 * h / a[p3 + iq];
            t = 1.0 / (fabs(theta)+
            (double)sqrt(1.0+theta*theta));
            if (theta < 0.0)
              t = -t;
          }
          c = 1.0 / (double)sqrt(1+t*t);
          s = t * c;
          tau = s / (1.0+c);
          h = t * a[p3 + iq];
          z[ip] -= h;
          z[iq] += h;
          d[ip] -= h;
          d[iq] += h;
          a[p3 + iq] = 0.0;
          for (j=j3=0; j<=ip-1; j++,j3+=3)
            ROTATE(a,j3,ip,j3,iq)
          for (j=ip+1; j<=iq-1; j++)
            ROTATE(a,p3,j,j*3,iq)
          for (j=iq+1; j<3; j++)
            ROTATE(a,p3,j,iq*3,j)

          for (j3=0; j3<9; j3+=3)
            ROTATE(v,j3,ip,j3,iq)

          ++(*nrot);
        }
      }
    }
    for (ip=0; ip<3; ip++) {
      b[ip] += z[ip];
      d[ip] = b[ip];
      z[ip] = 0.0;
    }
  }
  mprintf("Too many iterations in routine JACOBI\n");
  return(0);
}

// diagEsort()
/** Diagonalize 3x3 matrix, sort eigenvalues/eigenvectors.
  */
int diagEsort(double *mat, double *Emat, double *Evec[], double *Eigenvalue) {
  int njrot;
  int i, j, k, i3;
  double eigenvector[9], *eA, v;

  if (!jacobi3(mat, Eigenvalue, eigenvector, &njrot)) {
    mprintf("convergence failed\n");
    return(0);
  }

  for (i=i3=0; i<3; i++, i3+=3)
    for (j=0; j<3; j++)
      Emat[i3+j] = eigenvector[j*3+i];

  for (i=0; i<3; i++)
    Evec[i] = (double *) &Emat[i*3];

  for (i=0; i<2; i++) {
    v = Eigenvalue[k=i];
    for (j=i+1; j<3; j++)
      if (Eigenvalue[j] > v)
        v = Eigenvalue[k=j];
    if (k != i) {

      Eigenvalue[k] = Eigenvalue[i];
      Eigenvalue[i] = v;
      eA = Evec[i];
      Evec[i] = Evec[k];
      Evec[k] = eA;
    }
  }
  return(1);
}

// printVector()
void printVector(const char *Name, double V[3]) {
  mprintf("    %s: %8.4lf %8.4lf %8.4lf\n",Name,V[0], V[1], V[2]);
}

// printMatrix_3x3()
void printMatrix_3x3(const char *Title, double U[9]) {
  mprintf("    %s\n",Title);
  mprintf("     %8.4lf %8.4lf %8.4lf\n", U[0], U[1], U[2]);
  mprintf("     %8.4lf %8.4lf %8.4lf\n", U[3], U[4], U[5]);
  mprintf("     %8.4lf %8.4lf %8.4lf\n", U[6], U[7], U[8]);
}

// printMatrix()
void printMatrix(const char *Title, double *U, int mrows, int ncols) {
  mprintf("    %s",Title);
  int usize = mrows * ncols;
  for (int i = 0; i < usize; i++) {
    if ( (i%ncols)==0 ) mprintf("\n");
    mprintf(" %10.5lf",U[i]);
  }
  mprintf("\n");
}
    
// printRotTransInfo()
void printRotTransInfo(double U[9], double trans[6]) {

  printMatrix_3x3("Rotation matrix follows",U);
  printVector("Translation 1",trans);
  printVector("Translation 2",trans+3);
}

