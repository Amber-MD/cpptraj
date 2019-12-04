#include "ChiSqDist.h"
#include <cmath>
#include "CpptrajStdio.h"
#include "GammaFn.h"

double Cpptraj::Math::ChiSqDist(double chisq, int dof) {
  if (chisq < 0) {
    mprinterr("Error: ChiSqDist: Chi^2 value is < 0 (%g).\n", chisq);
    return 0;
  }
  if (dof < 2) {
    mprintf("Warning: ChiSqDist: Less than 2 degrees of freedom (%i)\n", dof);
  }

  double n_over_2 = (double)dof / 2.0;

  double two_no2 = pow(2.0, n_over_2);
  double chisq_no21 = pow(chisq, n_over_2 - 1.0);
  double exp_chisq = exp(-chisq / 2.0);
  double gamma_dof = GammaFn(n_over_2);

  double pval = (chisq_no21 * exp_chisq) / (two_no2 * gamma_dof);
  return pval;
}
