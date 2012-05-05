#include <cmath>
#include "Matrix_3x3.h"
#include "CpptrajStdio.h"
#include "vectormath.h"

Matrix_3x3::Matrix_3x3() {
  M_[0] = 0;
  M_[1] = 0;
  M_[2] = 0;
  M_[3] = 0;
  M_[4] = 0;
  M_[5] = 0;
  M_[6] = 0;
  M_[7] = 0;
  M_[8] = 0;
};

Matrix_3x3::Matrix_3x3(double *Min) {
  M_[0] = Min[0]; 
  M_[1] = Min[1]; 
  M_[2] = Min[2]; 
  M_[3] = Min[3]; 
  M_[4] = Min[4]; 
  M_[5] = Min[5]; 
  M_[6] = Min[6]; 
  M_[7] = Min[7]; 
  M_[8] = Min[8]; 
}

const int Matrix_3x3::MAX_ITERATIONS = 50;

#define ROTATE(ARR,MAJ1,MIN1,MAJ2,MIN2) { \
  dg = ARR[MAJ1 + MIN1]; \
  dh = ARR[MAJ2 + MIN2]; \
  ARR[MAJ1 + MIN1] = dg - ds*(dh+dg*tau); \
  ARR[MAJ2 + MIN2] = dh + ds*(dg-dh*tau); }

int Matrix_3x3::Diagonalize( double *vecD, double* matrixV ) {
  // Create identity matrix
  //double matrixV[9];
  matrixV[0] = 1;
  matrixV[1] = 0;
  matrixV[2] = 0;
  matrixV[3] = 0;
  matrixV[4] = 1;
  matrixV[5] = 0;
  matrixV[6] = 0;
  matrixV[7] = 0;
  matrixV[8] = 1;
  // Set vectors B and D equal to diagonal of M_. vector Z is 0.
  double vecB[3], vecZ[3];
  vecB[0] = vecD[0] = M_[0]; 
  vecB[1] = vecD[1] = M_[4]; 
  vecB[2] = vecD[2] = M_[8];
  vecZ[0] = 0; 
  vecZ[1] = 0; 
  vecZ[2] = 0;
  // MAIN LOOP
  double tresh = 0;
  int nrot;
  for (int i = 0; i < MAX_ITERATIONS; ++i) {
    // sm = SUM of UPPER RIGHT TRIANGLE
    double sm = fabs(M_[1]) + fabs(M_[2]) + fabs(M_[5]);
    if (sm == 0) return 0;
    
    if (i < 3)
      tresh = 0.2 * sm / 9;
    else
      tresh = 0;
    // BEGIN INNER LOOP OVER UPPER RIGHT TRIANGLE
    double dt;
    //int p3 = 0;
    int ip, p3;
    for ( ip = p3 = 0; ip < 2; ++ip, p3+=3) {
      for ( int iq = ip + 1; iq < 3; ++iq ) {
        int midx = p3 + iq;
        double dg = 100.0 * fabs(M_[midx]);
        if ( i > 3 && fabs(vecD[ip]) + dg == fabs(vecD[ip]) && 
                      fabs(vecD[iq]) + dg == fabs(vecD[iq]) )
        {
          M_[midx] = 0;
        } else if ( fabs(M_[midx]) > tresh) {
          double dh = vecD[iq] - vecD[ip];
          if (fabs(dh) + dg == fabs(dh))
            dt = M_[p3 + iq] / dh;
          else {
            double theta = 0.5 * dh / M_[midx];
            dt = 1.0 / (fabs(theta)+(double)sqrt(1.0+theta*theta));
            if (theta < 0.0)
              dt = -dt;
          }
          double dc = 1.0 / (double)sqrt(1+dt*dt);
          double ds = dt * dc;
          double tau = ds / (1.0+dc);
          dh = dt * M_[midx];
          vecZ[ip] -= dh;
          vecZ[iq] += dh;
          vecD[ip] -= dh;
          vecD[iq] += dh;
          M_[midx] = 0;
          int j, j3;
          for (j=j3=0; j<=ip-1; j++,j3+=3)
            ROTATE(M_,j3,ip,j3,iq)
          for (int j=ip+1; j<=iq-1; j++)
            ROTATE(M_,p3,j,j*3,iq)
          for (int j=iq+1; j<3; j++)
            ROTATE(M_,p3,j,iq*3,j)

          for (j3=0; j3<9; j3+=3)
            ROTATE(matrixV,j3,ip,j3,iq)

          ++nrot;
        }
      }
    } // END INNER LOOP OVER UPPER RIGHT TRIANGLE
    vecB[0] += vecZ[0];
    vecD[0] = vecB[0];
    vecZ[0] = 0;
    vecB[1] += vecZ[1];
    vecD[1] = vecB[1];
    vecZ[1] = 0;
    vecB[2] += vecZ[2];
    vecD[2] = vecB[2];
    vecZ[2] = 0;
  }
  mprintf("Too many iterations in routine!\n");
  return 1;
};

//int Matrix_3x3::Diagonalize_Sort(double (&EvecOut)[3][3], double *EvalOut) 
int Matrix_3x3::Diagonalize_Sort(double *EvecOut, double *EvalOut) 
{
  double Evec[9], Eval[3];

  if ( Diagonalize( Eval, Evec ) ) {
    mprintf("Convergence failed.\n");
    return 0; // 0 to be consistent with old diagEsort behavior
  }
  //printMatrix_3x3("Eigenvector Matrix", Evec);

  int i1,i2,i3;

  if (Eval[0] > Eval[1] && Eval[0] > Eval[2]) { // 0 is max
    if (Eval[1] > Eval[2]) {
      i1 = 0; i2 = 1; i3 = 2;
    } else {
      i1 = 0; i2 = 2; i3 = 1;
    }
  } else if (Eval[1] > Eval[0] && Eval[1] > Eval[2]) { // 1 is max
    if (Eval[0] > Eval[2]) {
      i1 = 1; i2 = 0; i3 = 2;
    } else {
      i1 = 1; i2 = 2; i3 = 0;
    }
  } else if (Eval[0] > Eval[1]) { // 2 is max
    i1 = 2; i2 = 0; i3 = 1;
  } else {
    i1 = 2; i2 = 1; i3 = 0;
  }
  //mprintf("EIGENVALUE ORDER (0=high, 3=med, 6=low): %i %i %i\n",i1,i2,i3);

  // Swap Eigenvectors - place them in rows
  EvecOut[0] = Evec[i1];
  EvecOut[1] = Evec[i1+3];
  EvecOut[2] = Evec[i1+6];

  EvecOut[3] = Evec[i2];
  EvecOut[4] = Evec[i2+3];
  EvecOut[5] = Evec[i2+6];

  EvecOut[6] = Evec[i3];
  EvecOut[7] = Evec[i3+3];
  EvecOut[8] = Evec[i3+6];

  // Swap eigenvalues
  EvalOut[0] = Eval[i1];
  EvalOut[1] = Eval[i2];
  EvalOut[2] = Eval[i3];

  return 1;
}
 
void Matrix_3x3::Print(const char* Title) {
  mprintf("    %s\n",Title);
  mprintf("     %8.4lf %8.4lf %8.4lf\n", M_[0], M_[1], M_[2]);
  mprintf("     %8.4lf %8.4lf %8.4lf\n", M_[3], M_[4], M_[5]);
  mprintf("     %8.4lf %8.4lf %8.4lf\n", M_[6], M_[7], M_[8]);
} 
