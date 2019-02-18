#include <cmath> // fabs, atan2, sin, cos
#include "Metric_Data.h"
#include "../Constants.h" // RADDEG, DEGRAD

/// Calculate smallest difference between two angles (in degrees).
double Cpptraj::Cluster::Metric_Data::DistCalc_Dih(double d1, double d2) {
  double diff = fabs(d1 - d2);
  if (diff > 180.0)
    return (360.0 - diff);
  else
    return diff;
}

/// Calculate basic difference.
double Cpptraj::Cluster::Metric_Data::DistCalc_Std(double d1, double d2) {
  return fabs(d1 - d2);
}

/** Calculate unambiguous average dihedral angle (in degrees) by converting to 
  * cartesian coords using x = cos(theta), y = sin(theta), and:
  *   tan(avgtheta) = avgy / avgx = SUM[sin(theta)] / SUM[cos(theta)]
  * See Eq. 2 from Altis et al., J. Chem. Phys., 126 p. 244111 (2007).
  */
double Cpptraj::Cluster::Metric_Data::AvgCalc_Dih(DataSet_1D const& dsIn, 
                                                  Cframes const& cframesIn,
                                                  double& sumx, double& sumy)
{
  sumy = 0.0;
  sumx = 0.0;
  // TODO: Convert angles to radians prior to this call?
  for (Cframes::const_iterator frm = cframesIn.begin(); frm != cframesIn.end(); ++frm) {
    double theta = dsIn.Dval( *frm ) * Constants::DEGRAD;
    sumy += sin( theta );
    sumx += cos( theta );
  }
  return atan2(sumy, sumx) * Constants::RADDEG; 
}

double Cpptraj::Cluster::Metric_Data::AvgCalc_Std(DataSet_1D const& dsIn,
                                                  Cframes const& cframesIn)
{
  double val = 0.0;
  for (Cframes::const_iterator frm = cframesIn.begin(); frm != cframesIn.end(); ++frm)
    val += dsIn.Dval( *frm );
  return (val / (double)cframesIn.size());
}

/* Update centroid value for adding/removing a frame.
 * \param fval value of frame being added/removed.
 * \param cval current centroid value.
 * \param isTorsion data is periodic.
 * \param oldSize Previous size of the centroid.
 * \param OP Operation being performed.
 */
double Cpptraj::Cluster::Metric_Data::DistCalc_FrameCentroid(double fval, double cval,
                                     bool isTorsion,
                                     double oldSize, CentOpType OP,
                                     double& sumx, double& sumy)
{
  double newcval;
  if (isTorsion) {
    double ftheta = fval * Constants::DEGRAD;
    if (OP == ADDFRAME) {
      sumy += sin( ftheta );
      sumx += cos( ftheta );
    } else { // SUBTRACTFRAME
      sumy -= sin( ftheta );
      sumx -= cos( ftheta );
    }
    newcval = atan2(sumy, sumx) * Constants::RADDEG;
  } else {
    newcval = cval * oldSize;
    if (OP == ADDFRAME) {
      newcval += fval;
      newcval /= ( oldSize + 1 );
    } else { // SUBTRACTFRAME
      newcval -= fval;
      newcval /= ( oldSize - 1 );
    }
  }
  return newcval;
}
