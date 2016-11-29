#include "Spline.h"
#include "CpptrajStdio.h"
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
  * \param Xin Input X values
  * \param Yin Corresponding Y values
  */
int Spline::CalcCoeff_Cubic(Darray const& Xin, Darray const& Yin) {
  if (Xin.size() != Yin.size()) {
    mprinterr("Error: Spline::CalcCoeff_Cubic: X size (%zu) != Y size (%zu)\n", Xin.size(), Yin.size());
    return 1;
  }
  if (Xin.size() < 2) {
    mprinterr("Error: Spline::CalcCoeff_Cubic: Cannot spline fewer than 2 values.\n");
    return 1;
  }
  int n = (int)Xin.size();

  B_.resize(n, 0.0);
  C_.resize(n, 0.0);
  D_.resize(n, 0.0);

  int n_minus1 = n - 1;

  if ( n > 2 ) {
    // Generate Tri-diagonal matrix
    D_[0] = Xin[1] - Xin[0];
    C_[1] = (Yin[1] - Yin[0]) / D_[0];
    for (int i = 1; i < n_minus1; i++) {
      D_[i]   = Xin[i + 1] - Xin[i];
      B_[i]   = 2.0 * (D_[i - 1] + D_[i]);
      C_[i+1] = (Yin[i + 1] - Yin[i]) / D_[i];
      C_[i]   = C_[i+1] - C_[i];
      //mprintf("TRIDBG: %i b=%lf c=%lf d=%lf\n",i,b[i],c[i],d[i]);
    }

    // Set up boundary  conditions
    B_[0]        = -D_[0];
    B_[n_minus1] = -D_[n - 2];
    C_[0]        = 0.0;
    C_[n_minus1] = 0.0;
    if (n > 3) {
      C_[0]        = C_[2] / (Xin[3] - Xin[1]) - C_[1] / (Xin[2] - Xin[0]);
      C_[n_minus1] = C_[n-2] / (Xin[n_minus1] - Xin[n-3]) - C_[n-3] / (Xin[n-2] - Xin[n-4]);
      C_[0]        = C_[0] * D_[0] * D_[0] / (Xin[3] - Xin[0]);
      C_[n_minus1] = -C_[n_minus1] * D_[n-2] * D_[n-2] / (Xin[n_minus1] - Xin[n-4]);
    }

    // Forward elimination
    for (int i = 1; i < n; i++) {
      double t  = D_[i - 1] / B_[i - 1];
      B_[i]     = B_[i] - t * D_[i - 1];
      C_[i]     = C_[i] - t * C_[i - 1];
      //mprintf("FWDDBG: %i b=%lf c=%lf t=%lf\n",i,b[i],c[i],t);
    }

    // Back substitution
    C_[n_minus1] = C_[n_minus1] / B_[n_minus1];
    for (int i = n - 2; i > -1; i--) {
      C_[i] = (C_[i] - D_[i] * C_[i + 1]) / B_[i];
      //mprintf("BAKDBG: %i c=%lf\n",i,c[i]);
    }
    // Calculate the polynomial coefficients
    B_[n_minus1] = (Yin[n_minus1] - Yin[n-2]) / D_[n-2] + D_[n-2] * (C_[n-2] + 2.0 * C_[n_minus1]);
    for (int i = 0; i < n_minus1; i++) {
        B_[i] = (Yin[i+1] - Yin[i]) / D_[i] - D_[i] * (C_[i+1] + 2.0 * C_[i]);
        D_[i] = (C_[i+1] - C_[i]) / D_[i];
        C_[i] = 3.0 * C_[i];
        //mprintf("POLYDBG: %i b=%lf c=%lf d=%lf\n",i,B_[i],C_[i],D_[i]);
    }
    C_[n_minus1] = 3.0 * C_[n_minus1];
    D_[n_minus1] = D_[n - 2];
  } else {
    // Special case; n == 2
    B_[0] = (Yin[1] - Yin[0]) / (Xin[1] - Xin[0]);
    C_[0] = 0.0;
    D_[0] = 0.0;
    B_[1] = B_[0];
    C_[1] = 0.0;
    D_[1] = 0.0;
  }
  return 0;
}

double Spline::SplinedY(double Uin, Darray const& Xin, Darray const& Yin) const {
  int xidx = 0;
  int n = (int)Xin.size();
  // Search for Uin in x
  if (Uin < Xin[0])
    xidx = 0;
  else if (Uin > Xin[n-1])
    xidx = n - 1;
  else {
    int i0 = 0;
    int i1 = n - 1;
    while (i0 <= i1) {
      xidx = (i0 + i1) / 2;
      if ( Uin < Xin[xidx] )
        i1 = xidx - 1;
      else if ( Uin > Xin[xidx+1] )
        i0 = xidx + 1;
      else
        break;
    }
  }
  // Evaluate v for this u
  double dx = Uin - Xin[xidx];
  return  (Yin[xidx] + dx*(B_[xidx] + dx*(C_[xidx] + dx*D_[xidx])));
}

Spline::Darray Spline::SplinedYvals(Darray const& mesh_x, Darray const& Xin, Darray const& Yin) const
{
  Darray mesh_y(mesh_x.size(), 0.0);
  for (unsigned int uidx = 0; uidx < mesh_x.size(); uidx++)
    mesh_y[uidx] = SplinedY(mesh_x[uidx], Xin, Yin);
  return mesh_y;
}
