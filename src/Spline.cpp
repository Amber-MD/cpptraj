#include "Spline.h"
#include "CpptrajStdio.h"

// Spline::cubicSpline_coeff()
/** Given a set of x and y values of size n, compute the b, c, and d
  * coefficients for n interpolating cubic splines of form:
  *
  *   Si(t) = y[i] + b_[i]*dt + c_[i]*dt^2 + d_[i]*dt^3 
  *
  * where dt = t - x[i]. Si(t), and first and second derivatives Si'(t)
  * and Si"(t) must be continuous over interval dt, and both derivatives 
  * for adjacent points must be equal so that adjacent line segments become 
  * continuous.
  *
  *   Si'(t) = b_[i] + 2*c_[i]*dt + 3*d_[i]*dt^2
  *   Si"(t) = 2*c_[i] + 6*d_[i]*dt
  *
  * \param xIn Input X values
  * \param yIn Y values corresponding to input X values.
  */
void Spline::cubicSpline_coeff(Darray const& xIn, Darray const& yIn) 
{
  if (xIn.size() < 2) {
    mprinterr("Error: Cannot spline with less than 2 values.\n");
    return;
  }
  int n = (int)xIn.size();

  b_.resize(n, 0.0);
  c_.resize(n, 0.0);
  d_.resize(n, 0.0);

  int n_minus1 = n - 1;

  if ( n > 2 ) {
    // Generate Tri-diagonal matrix
    d_[0] = xIn[1] - xIn[0];
    c_[1] = (yIn[1] - yIn[0]) / d_[0];
    for (int i = 1; i < n_minus1; i++) {
      d_[i] = xIn[i + 1] - xIn[i];
      b_[i] = 2.0 * (d_[i - 1] + d_[i]);
      c_[i+1] = (yIn[i + 1] - yIn[i]) / d_[i];
      c_[i] = c_[i+1] - c_[i];
      //mprintf("TRIDBG: %i b=%lf c=%lf d=%lf\n",i,b_[i],c_[i],d_[i]);
    }
    
    // Set up boundary  conditions
    b_[0]        = -d_[0];
    b_[n_minus1] = -d_[n-2];
    c_[0]        = 0.0;
    c_[n_minus1] = 0.0;
    if (n > 3) {
      c_[0]        = c_[2] / (xIn[3] - xIn[1]) - c_[1] / (xIn[2] - xIn[0]);
      c_[n_minus1] = c_[n-2] / (xIn[n_minus1] - xIn[n-3]) - c_[n-3] / (xIn[n-2] - xIn[n-4]);
      c_[0]        = c_[0] * d_[0] * d_[0] / (xIn[3] - xIn[0]);
      c_[n_minus1] = -c_[n_minus1] * d_[n-2] * d_[n-2] / (xIn[n_minus1] - xIn[n-4]);
    }

    // Forward elimination
    for (int i = 1; i < n; i++) {
        double t = d_[i - 1] / b_[i - 1];
        b_[i]     = b_[i] - t * d_[i - 1];
        c_[i]     = c_[i] - t * c_[i - 1];
        //mprintf("FWDDBG: %i b=%lf c=%lf t=%lf\n",i,b_[i],c_[i],t);
    }

    // Back substitution
    c_[n_minus1] = c_[n_minus1] / b_[n_minus1];
    for (int i = n - 2; i > -1; i--) {
        c_[i] = (c_[i] - d_[i] * c_[i + 1]) / b_[i];
        //mprintf("BAKDBG: %i c=%lf\n",i,c_[i]);
    }

    // Calculate the polynomial coefficients
    b_[n_minus1] = (yIn[n_minus1] - yIn[n-2]) / d_[n-2] + d_[n-2] * (c_[n-2] + 2.0 * c_[n_minus1]);
    for (int i = 0; i < n_minus1; i++) {
        b_[i] = (yIn[i+1] - yIn[i]) / d_[i] - d_[i] * (c_[i+1] + 2.0 * c_[i]);
        d_[i] = (c_[i+1] - c_[i]) / d_[i];
        c_[i] = 3.0 * c_[i];
        //mprintf("POLYDBG: %i b=%lf c=%lf d=%lf\n",i,b_[i],c_[i],d_[i]);
    }
    c_[n_minus1] = 3.0 * c_[n_minus1];
    d_[n_minus1] = d_[n-2];

  // Special case; n == 2
  } else {
    b_[0] = (yIn[1] - yIn[0]) / (xIn[1] - xIn[0]);
    c_[0] = 0.0;
    d_[0] = 0.0;
    b_[1] = b_[0];
    c_[1] = 0.0; 
    d_[1] = 0.0;
  }
}

// Spline::cubicSpline_eval() 
/** Evaluate cubic spline function with pre-calcd coefficients in b, c, and
  * d from coordinates x/y for all points in given mesh.
  * \param xIn Original input X coordinates.
  * \param yIn Original input Y coordinates.
  * \param mesh_x New input X coordinates.
  * \return New output Y coordinates.
  */
Spline::Darray Spline::cubicSpline_eval(Darray const& xIn, Darray const& yIn, 
                                        Darray const& mesh_x) const
{
  int xidx = 0;
  int n = (int)xIn.size();
  int mesh_size = (int)mesh_x.size();
  Darray mesh_y;
  mesh_y.reserve( mesh_size );

  for (int uidx = 0; uidx < mesh_size; uidx++) {
    double U = mesh_x[uidx];
    // Search for U in x
    if (U < xIn[0])
      xidx = 0;
    else if (U > xIn[n-1])
      xidx = n - 1;
    else {
      int i0 = 0;
      int i1 = n - 1;
      while (i0 <= i1) {
        xidx = (i0 + i1) / 2;
        if ( U < xIn[xidx] )
          i1 = xidx - 1;
        else if ( U > xIn[xidx+1] )
          i0 = xidx + 1;
        else
          break;
      }
    }
    // Evaluate v for this u
    double dx = U - xIn[xidx];
    mesh_y.push_back( yIn[xidx] + dx*(b_[xidx] + dx*(c_[xidx] + dx*d_[xidx])) ); 
  }
  return mesh_y;
}
