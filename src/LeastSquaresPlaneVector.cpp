#include "LeastSquaresPlaneVector.h"
#include "Frame.h"
#include "AtomMask.h"
#include "Constants.h"
#include <cmath>
#include <algorithm> //sort

/** CONSTRUCTOR */
LeastSquaresPlaneVector::LeastSquaresPlaneVector()
{}

/** Reserve space for selected number of atoms. */
void LeastSquaresPlaneVector::ReserveForNumAtoms(unsigned int nselected) {
  vcorr_.reserve( nselected );
}

/** Solves a cubic equation: ax^3 + bx^2 + cx + d = 0 using "Cardan's formula" 
  * (see: Bronstein, S.131f)
  */
double LeastSquaresPlaneVector::solve_cubic_eq(double a, double b, double c, double d) {
  const double one3 = 1.0 / 3.0;
  const double one27 = 1.0 / 27.0;
  double droot = 0;

  double r, s, t;
  double p, q, rho, phi;
  double D, u, v;
  std::vector<double> dtmp(3);

  /* Coeff. for normal form x^3 + rx^2 + sx + t = 0 */
  r = b / a;
  s = c / a;
  t = d / a;

  /* Coeff. for red. eq. y^3 + py + q = 0 with y = x + r/3 bzw. (x = y - r/3) */
  p = s - r * r * one3;
  q = 2.0 * r * r * r * one27 - r * s * one3 + t;

  /* Dummy variables */
  rho = sqrt(-p * p * p * one27);
  phi = acos(-q / (2.0 * rho));

  /* Discriminante(?) */
  D = pow((p * one3),3) + q * q * 0.25;

  if(D > 0){ /* x real -> one real solution */
    u = pow(-q * 0.5 + sqrt(D), one3);
    v = -p / u * one3;
    droot = (u + v) - r * one3;
  } else { // D <= 0
  /* three real solutions (d < 0) | one real solution + one real double solution or 
                                                     one real triple solution (d = 0) */
    dtmp[0] = 2.0 * pow(rho, one3) * cos(phi * one3) - r * one3;
    dtmp[1] = 2.0 * pow(rho, one3) * cos((phi + Constants::TWOPI ) * one3) - r * one3;
    dtmp[2] = 2.0 * pow(rho, one3) * cos((phi + Constants::FOURPI) * one3) - r * one3;

    sort(dtmp.begin(), dtmp.end());

    //qsort((void *) dtmp, (size_t) 3, sizeof(double), cmpdouble);
    droot = dtmp[0];
  }
  return droot;
}

/** Calcs (least-squares best) plane through a series of points
  * relative to their center of geom. (the latter has to be done outside this routine), 
  * returns (normalized) coeff. for plane eq. ax + by + cz = 0
  * following: Crystal Structure Analysis for Chem. and Biol.,
  * Glusker, Lewis, Rossi, S. 460ff
  */
Vec3 LeastSquaresPlaneVector::leastSquaresPlane(int n, const double* vcorr) {
  double Xout, Yout, Zout;
  if (n == 9) {  // Special case, only 3 coords
    double x1 = vcorr[3] - vcorr[0];
    double y1 = vcorr[4] - vcorr[1];
    double z1 = vcorr[5] - vcorr[2];
    double x2 = vcorr[6] - vcorr[3];
    double y2 = vcorr[7] - vcorr[4];
    double z2 = vcorr[8] - vcorr[5];

    Xout = y1 * z2 - z1 * y2;
    Yout = z1 * x2 - x1 * z2;
    Zout = x1 * y2 - y1 * x2;
  } else { // General case
    // Calc Var.
    double dSumXX = 0.0;
    double dSumYY = 0.0;
    double dSumZZ = 0.0;
    double dSumXY = 0.0;
    double dSumXZ = 0.0;
    double dSumYZ = 0.0;

    for (int i = 0; i < n; i+=3) {
      dSumXX += vcorr[i  ] * vcorr[i  ];
      dSumYY += vcorr[i+1] * vcorr[i+1];
      dSumZZ += vcorr[i+2] * vcorr[i+2];

      dSumXY += vcorr[i  ] * vcorr[i+1];
      dSumXZ += vcorr[i  ] * vcorr[i+2];
      dSumYZ += vcorr[i+1] * vcorr[i+2];
    }

    // Calc coeff. for -l^3 + o * l^2 + p * l + q = 0
    double o = dSumXX + dSumYY + dSumZZ;
    double p = pow(dSumXY,2) + pow(dSumXZ,2) + pow(dSumYZ,2) -
               (dSumXX * dSumYY + dSumXX * dSumZZ + dSumYY * dSumZZ);
    double q = dSumXX * dSumYY * dSumZZ + 2.0 * dSumXY * dSumXZ * dSumYZ -
               (dSumXX * dSumYZ * dSumYZ + dSumYY * dSumXZ * dSumXZ + dSumZZ * dSumXY * dSumXY);

    // Solve cubic eq.
    double root = solve_cubic_eq(-1.0, o, p, q);

    // Calc determinents
    Xout = (dSumYY - root) * dSumXZ - dSumXY * dSumYZ;
    Yout = (dSumXX - root) * dSumYZ - dSumXY * dSumXZ;
    Zout =  dSumXY         * dSumXY - (dSumYY - root) * (dSumXX - root);
  }
  // Normalize
  double dnorm = 1.0 / sqrt((Xout * Xout) + (Yout * Yout) + (Zout * Zout));
  return Vec3(Xout * dnorm, Yout * dnorm, Zout * dnorm);
}

