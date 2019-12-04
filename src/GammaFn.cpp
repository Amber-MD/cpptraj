#include "GammaFn.h"
#include <cmath>
#ifndef C11_SUPPORT
#include <cfloat> // DBL_MAX
#include "Constants.h"
#include "CpptrajStdio.h"
  // The code below is adapted from the Gamma code found here:
  //   https://www.johndcook.com/blog/stand_alone_code/
  // It is only intended for use when C++11 support is not available since
  // the errors from the "real" gamma function increase for large X
  // (about X=14 and above).
  // When C++11 support is available, the functions from cmath are used.
#endif

/** \return the Gamma function value of xIn.
  * The Gamma function is defined as:
  *   GammaFn(x)   = Integral(0, inf)[ exp(-t) * t^(x-1) dt ] for x > 0
  *   GammaFn(n+1) = n!                                       for n = 1,2,3...
  *   GammaFn(.5)  = sqrt(PI)
  */
double Cpptraj::Math::GammaFn(double xIn) {
# ifdef C11_SUPPORT
  return tgamma(xIn);
# else
  if (xIn <= 0) {
    mprinterr("Error: GammaFn argument is <= 0 (%g)\n", xIn);
    return 0;
  }

  // ---------------------------------------------
  if (xIn < 0.001) {
    // For small x, 1/Gamma(x) has power series x + gamma x^2  - ...
    // So in this range, 1/Gamma(x) = x + gamma x^2 with error on the order of x^3.
    // The relative error over this interval is less than 6e-7.
    return 1.0 / (xIn * (1.0 + Constants::EULER_MASC * xIn));
  }

  // ---------------------------------------------
  if (xIn < 12.0) {
    // The algorithm directly approximates gamma over (1,2) and uses
    // reduction identities to reduce other arguments to this interval.
    double y = xIn;
    int n = 0;
    bool lt_one = (xIn < 1.0);

    // Add or subtract integers as necessary to bring y into 1 < y < 2
    if (lt_one)
      y += 1.0;
    else {
      n = (int)(floor(y)) - 1;
      y -= n;
    }

    // numerator coefficients for approximation over 1 < y < 2
    static const double p[] =
    {
            -1.71618513886549492533811E+0,
             2.47656508055759199108314E+1,
            -3.79804256470945635097577E+2,
             6.29331155312818442661052E+2,
             8.66966202790413211295064E+2,
            -3.14512729688483675254357E+4,
            -3.61444134186911729807069E+4,
             6.64561438202405440627855E+4
    };

    // denominator coefficients for approximation over the interval (1,2)
    static const double q[] =
    {
            -3.08402300119738975254353E+1,
             3.15350626979604161529144E+2,
            -1.01515636749021914166146E+3,
            -3.10777167157231109440444E+3,
             2.25381184209801510330112E+4,
             4.75584627752788110767815E+3,
            -1.34659959864969306392456E+5,
            -1.15132259675553483497211E+5
    };
    double num = 0.0;
    double den = 1.0;
    int i;

    double z = y - 1;
    for (i = 0; i < 8; i++)
    {
      num = (num + p[i])*z;
      den = den*z + q[i];
    }
    double result = num/den + 1.0;

    // Apply correction if argument was not initially in (1,2)
    if (lt_one) {
      // Use identity gamma(z) = gamma(z+1)/z
      // The variable "result" now holds gamma of the original y + 1
      // Thus we use y-1 to get back the orginal y.
      result /= (y-1.0);
    } else {
      // Use the identity gamma(z+n) = z*(z+1)* ... *(z+n-1)*gamma(z)
      for (i = 0; i < n; i++)
        result *= y++;
    }

    return result;
  } // END if xIn < 12

  // ---------------------------------------------
  if (xIn > 171.624) {
    // Answer is out of range. Return +inf
    double temp = DBL_MAX;
    return 2 * temp;
  }

  // ---------------------------------------------
  return exp(LogGammaFn(xIn));
# endif /* C11_SUPPORT */
}

/** \return the natural log-Gamma function value of xIn. */
double Cpptraj::Math::LogGammaFn(double xIn) {
# ifdef C11_SUPPORT
  return lgamma( xIn );
# else
  if (xIn <= 0) {
    mprinterr("Error: GammaFn argument is <= 0 (%g)\n", xIn);
    return 0;
  }

  if (xIn < 12.0) {
    return log(fabs(GammaFn(xIn)));
  }

  // Abramowitz and Stegun 6.1.41
  // Asymptotic series should be good to at least 11 or 12 figures
  // For error analysis, see Whittiker and Watson
  // A Course in Modern Analysis (1927), page 252

  static const double c[8] =
  {
        1.0/12.0,
       -1.0/360.0,
        1.0/1260.0,
       -1.0/1680.0,
        1.0/1188.0,
     -691.0/360360.0,
        1.0/156.0,
    -3617.0/122400.0
  };
  double z = 1.0/(xIn*xIn);
  double sum = c[7];
  for (int i=6; i >= 0; i--)
  {
    sum *= z;
    sum += c[i];
  }
  double series = sum/xIn;

  static const double halfLogTwoPi = 0.91893853320467274178032973640562;
  double logGamma = (xIn - 0.5)*log(xIn) - xIn + halfLogTwoPi + series;    
  return logGamma;
# endif /* C11_SUPPORT */
}

