#include "Integrate.h"
#include "CpptrajStdio.h"

/// CONSTRUCTOR - Set up X mesh
Interpolate::Interpolate(double ti, double tf, int n) :
  mesh_size_(n),
  mesh_x_(n),
  mesh_y_(n)
{
  CalcMeshX(ti, tf);
}

// Interpolate::CalcMeshX()
void Interpolate::CalcMeshX(double ti, double tf) {
  double s = (ti + tf)/2;
  double d = (tf - ti)/2;
  for (int i = 0; i < mesh_size_; i++)
    mesh_x_[i] = s + d*((double) (2*i + 1 - mesh_size_)/(mesh_size_ - 1));
}

// Interpolate::cubicSpline_coeff()
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
  */
void Interpolate::cubicSpline_coeff(std::vector<double> const& x, std::vector<double> const& y) 
{
  int n = (int)x.size();

  b.resize(n, 0.0);
  c.resize(n, 0.0);
  d.resize(n, 0.0);

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

// Interpolate::cubicSpline_eval() 
/** Evaluate cubic spline function with pre-calcd coefficients in b, c, and
  * d from coordinates x/y for all points in mesh.
  * \param x Input X coordinates
  * \param y Input Y coordinates
  */
void Interpolate::cubicSpline_eval(std::vector<double> const& x, std::vector<double> const& y)
{
  int xidx;  
  int n = (int)x.size();

  for (int uidx = 0; uidx < mesh_size_; uidx++) {
    double U = mesh_x_[uidx];
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
    mesh_y_[uidx] = y[xidx] + dx*(b[xidx] + dx*(c[xidx] + dx*d[xidx])); 
  }
}

// Interpolate::SetMesh_X() 
/** Given a start, stop, and size, set up mesh x values. */
void Interpolate::SetMesh_X(double ti, double tf, int n) {
  mesh_size_ = n;
  mesh_x_.resize( mesh_size_ );
  mesh_y_.resize( mesh_size_ );
  CalcMeshX(ti, tf);
}

// Interpolate::SetMesh_Y()
/** Set mesh Y values based on input X and Y and previously set up
  * mesh X.
  */
int Interpolate::SetMesh_Y(std::vector<double> const& x, std::vector<double> const& y) {
  if (x.size() != y.size()) {
    mprinterr("Error: Interpolate::SetMesh: X size (%u) != Y size (%u)\n",
              x.size(), y.size());
    return 1;
  }
  // No point if 1 or less values
  if (x.size() < 2) {
    mprinterr("Error: Interpolate::SetMesh: Requires > 1 values (%u specified).\n",
              x.size());
    return 1;
  }
  cubicSpline_coeff(x, y); 
  cubicSpline_eval(x, y);
  return 0;
}

// Interpolate::Integrate_Trapezoid() 
/** Integrate the mesh using the trapezoid rule. */
double Interpolate::Integrate_Trapezoid() {
  double sum = 0.0;

  if (mesh_size_ < 2) return 0.0;
  
  for (int i = 1; i < mesh_size_; i++) {
      double b_minus_a = (mesh_x_[i] - mesh_x_[i - 1]);
      sum += (b_minus_a * (mesh_y_[i - 1] + mesh_y_[i]) * 0.5);
  }
  return sum;
}
