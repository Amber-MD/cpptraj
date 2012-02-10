#include "Integrate.h"
#include "CpptrajStdio.h"
/*! \file Integrate.cpp
    \brief Routines for interpolation and integration.
 */

// cubicSpline_coeff
/** Given a set of x and y values of size n, compute the b, c, and d
  * coefficients for n interpolating cubic splines of form:
  *
  *   Si(t) = y[i] + b[i]*dt + c[i]*dt^2 + d[i]*dt^3 
  *
  * where dt = t - x[i]. Si(t), and first and second derivatives Si'(t)
  * and Si"(t) must be continuous over interval dt, and both derivatives 
  * for adjacent points must be equal so that adjacent line segments become 
  * continuous.
  *
  *   Si'(t) = b[i] + 2*c[i]*dt + 3*d[i]*dt^2
  *   Si"(t) = 2*c[i] + 6*d[i]*dt
  *
  * \param x Input X values
  * \param y Corresponding Y values
  * \param n Total number of X/Y values
  * \param b Output b coefficients
  * \param c Output c coefficients
  * \param d Output d coefficients
  */
void Interpolate::cubicSpline_coeff(double *x, double *y, int n) 
{
  // NOTE: CHECK FOR NULLS!
  // No point if 1 or less values
  if (n < 2) return;

  b.resize(n,0);
  c.resize(n,0);
  d.resize(n,0);

  int n_minus1 = n - 1;

  if ( n > 2 ) {
    // Generate Tri-diagonal matrix
    d[0] = x[1] - x[0];
    c[1] = (y[1] - y[0]) / d[0];
    for (int i = 1; i < n_minus1; i++) {
      d[i] = x[i + 1] - x[i];
      b[i] = 2.0 * (d[i - 1] + d[i]);
      c[i+1] = (y[i + 1] - y[i]) / d[i];
      c[i] = c[i+1] - c[i];
      //mprintf("TRIDBG: %i b=%lf c=%lf d=%lf\n",i,b[i],c[i],d[i]);
    }
    
    // Set up boundary  conditions
    b[0]        = -d[0];
    b[n_minus1] = -d[n - 2];
    c[0]        = 0.0;
    c[n_minus1] = 0.0;
    if (n > 3) {
      c[0]        = c[2] / (x[3] - x[1]) - c[1] / (x[2] - x[0]);
      c[n_minus1] = c[n - 2] / (x[n_minus1] - x[n - 3]) - c[n - 3] / (x[n - 2] - x[n - 4]);
      c[0]        = c[0] * d[0] * d[0] / (x[3] - x[0]);
      c[n_minus1] = -c[n_minus1] * d[n - 2] * d[n - 2] / (x[n_minus1] - x[n - 4]);
    }

    // Forward elimination
    for (int i = 1; i < n; i++) {
        double t = d[i - 1] / b[i - 1];
        b[i]     = b[i] - t * d[i - 1];
        c[i]     = c[i] - t * c[i - 1];
        //mprintf("FWDDBG: %i b=%lf c=%lf t=%lf\n",i,b[i],c[i],t);
    }

    // Back substitution
    c[n_minus1] = c[n_minus1] / b[n_minus1];
    for (int i = n - 2; i > -1; i--) {
        c[i] = (c[i] - d[i] * c[i + 1]) / b[i];
        //mprintf("BAKDBG: %i c=%lf\n",i,c[i]);
    }

    // Calculate the polynomial coefficients
    b[n_minus1] = (y[n_minus1] - y[n - 2]) / d[n - 2] + d[n - 2] * (c[n - 2] + 2.0 * c[n_minus1]);
    for (int i = 0; i < n_minus1; i++) {
        b[i] = (y[i + 1] - y[i]) / d[i] - d[i] * (c[i + 1] + 2.0 * c[i]);
        d[i] = (c[i + 1] - c[i]) / d[i];
        c[i] = 3.0 * c[i];
        //mprintf("POLYDBG: %i b=%lf c=%lf d=%lf\n",i,b[i],c[i],d[i]);
    }
    c[n_minus1] = 3.0 * c[n_minus1];
    d[n_minus1] = d[n - 2];

  // Special case; n == 2
  } else {
    b[0] = (y[1] - y[0]) / (x[1] - x[0]);
    c[0] = 0.0;
    d[0] = 0.0;
    b[1] = b[0];
    c[1] = 0.0; 
    d[1] = 0.0;
  }
}

// cubicSpline_eval()
/** Evaluate cubic spline function with pre-calcd coefficients in b, c, and
  * d from coordinates x/y for all points in u.
  * \param u Points at which spline function will be evaluated.
  * \param v Evaluated points.
  * \param ulen Number of points to be evaluated
  * \param x Input X coordinates
  * \param y Input Y coordinates
  * \param b Precalcd b coefficients
  * \param c Precalcd c coefficients
  * \param d Precalcd d coefficients
  * \param n number of data points in x,y,b,c,d
  * \return 0 on success, 1 on error.
  */
int Interpolate::cubicSpline_eval(double *u, double *v, int ulen, double *x, double *y, int n)
{
  int xidx;  

  if (n < 2) {
    mprinterr("Error: cubicSpline_eval: Data length is < 2.\n");
    return 1;
  }

  for (int uidx = 0; uidx < ulen; uidx++) {
    double U = u[uidx];

    // Search for U in x
    if (U < x[0])
      xidx = 0;
    else if (U > x[n-1])
      xidx = n - 1;
    else {
      int i0 = 0;
      int i1 = n - 1;
      while (i0 <= i1) {
        xidx = (i0 + i1) / 2;
        if ( U < x[xidx] )
          i1 = xidx - 1;
        else if ( U > x[xidx+1] )
          i0 = xidx + 1;
        else
          break;
      }
    }
    
    // Evaluate v for this u
    double dx = U - x[xidx];
    v[uidx] = y[xidx] + dx*(b[xidx] + dx*(c[xidx] + dx*d[xidx])); 
  }
  return 0;
}

// integrate_trapezoid()
/** Integrate the set x,y of size n using the trapezoid rule.
  */
double integrate_trapezoid(double *x, double *y, int n) {
  double sum = 0.0;

  if (n < 2) return 0.0;
  
  for (int i = 1; i < n; i++) {
      double b_minus_a = (x[i] - x[i - 1]);
      sum += (b_minus_a * (y[i - 1] + y[i]) * 0.5);
  }
  return sum;
}

// set_xvalues_range()
/** Given a start, stop, and size, generate x values.
  */
void set_xvalues_range(double *x_out, double ti, double tf, int n) {
  double s = (ti + tf)/2;
  double d = (tf - ti)/2;
  for (int i = 0; i < n; i++)
    x_out[i] = s + d*((double) (2*i + 1 - n)/(n - 1));
}

