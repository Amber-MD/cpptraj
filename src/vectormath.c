/*
 * vectormath.c
 * Simple vector/matrix math routines.
 * ROTATE, jacobi3, and diagEsort routines adapted from PTRAJ
 */
#include <stdio.h>
#include <math.h>

/*
 * normalize()
 * Normalize vector in a[].
 */
void normalize(double a[3]) {
  double b;

  b = 1.0/sqrt((double)(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]));

  a[0] *= b;
  a[1] *= b;
  a[2] *= b;
}

/*
 * ROTATE()
 */
#define ROTATE(ARR,MAJ1,MIN1,MAJ2,MIN2) { \
  g = ARR[MAJ1 + MIN1]; \
  h = ARR[MAJ2 + MIN2]; \
  ARR[MAJ1 + MIN1] = g - s*(h+g*tau); \
  ARR[MAJ2 + MIN2] = h + s*(g-h*tau); }

/*
 * jacobi3()
 */
int jacobi3(double *a, double *d, double *v, int *nrot) { 
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
  printf("Too many iterations in routine JACOBI\n");
  return(0);
}

/*
 * diagEsort()
 */
int diagEsort(double *mat, double *Emat, double *Evec[], double *Eigenvalue) {
  int njrot;
  int i, j, k, i3;
  double eigenvector[9], *eA, v;

  if (!jacobi3(mat, Eigenvalue, eigenvector, &njrot)) {
    printf("convergence failed\n");
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

/*
 * printRotTransInfo()
 */
void printRotTransInfo(double U[9], double trans[6]) {

  fprintf(stdout, "    Rotation matrix follows\n");
  fprintf(stdout, "     %10.8f %10.8f %10.8f\n",
          U[0], U[1], U[2]);
  fprintf(stdout, "     %10.8f %10.8f %10.8f\n",
          U[3], U[4], U[5]);
  fprintf(stdout, "     %10.8f %10.8f %10.8f\n",
          U[6], U[7], U[8]);
  fprintf(stdout, "    Translation 1 is %10.8f %10.8f %10.8f\n",
          trans[0], trans[1], trans[2]);
  fprintf(stdout, "    Translation 2 is %10.8f %10.8f %10.8f\n",
          trans[3], trans[4], trans[5]);

  return;
}

