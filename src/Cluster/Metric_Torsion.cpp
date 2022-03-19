#include "Metric_Torsion.h"
#include "Centroid_Num.h"
#include "Cframes.h"
#include "../Constants.h"
#include "../CpptrajStdio.h"
#include "../DataSet_1D.h"
#include <cmath> // fabs, atan2

/** CONSTRUCTOR */
Cpptraj::Cluster::Metric_Torsion::Metric_Torsion() :
  Metric(TORSION),
  data_(0)
{}

/** Initialize. */
int Cpptraj::Cluster::Metric_Torsion::Init(DataSet_1D* dataIn) {
  if (dataIn == 0) {
    mprinterr("Internal Error: Metric_Torsion::Init called with null set.\n");
    return 1;
  }
  data_ = dataIn;
  return 0;
}

/** Set up. Not much to do for a DataSet. */
int Cpptraj::Cluster::Metric_Torsion::Setup() {
  return 0;
}

/// \return smallest difference between two angles (in degrees).
static inline double DistCalc_Dih(double d1, double d2) {
  double diff = fabs(d1 - d2);
  if (diff > 180.0)
    return (360.0 - diff);
  else
    return diff;
}

/** \return Shortest absolute difference between torsions */
double Cpptraj::Cluster::Metric_Torsion::FrameDist(int f1, int f2) {
  return DistCalc_Dih( data_->Dval(f1), data_->Dval(f2) );
}

/** \return Absolute difference between centroids. */
double Cpptraj::Cluster::Metric_Torsion::CentroidDist(Centroid* c1, Centroid* c2) {
  return DistCalc_Dih( ((Centroid_Num*)c1)->Cval(), ((Centroid_Num*)c2)->Cval() );
}

/** \return Absolute difference between point and centroid. */
double  Cpptraj::Cluster::Metric_Torsion::FrameCentroidDist(int f1, Centroid* c1) {
  return DistCalc_Dih( data_->Dval(f1), ((Centroid_Num*)c1)->Cval() );
}

/** Calculate centroid from specified frames.
  * Calculate unambiguous average dihedral angle (in degrees) by converting to 
  * cartesian coords using x = cos(theta), y = sin(theta), and:
  *   tan(avgtheta) = avgy / avgx = SUM[sin(theta)] / SUM[cos(theta)]
  * See Eq. 2 from Altis et al., J. Chem. Phys., 126 p. 244111 (2007).
  */
void Cpptraj::Cluster::Metric_Torsion::CalculateCentroid(Centroid* centIn, Cframes const& cframesIn)
{
  double sumy = 0.0;
  double sumx = 0.0;
  // TODO: Convert angles to radians prior to this call?
  for (Cframes::const_iterator frm = cframesIn.begin(); frm != cframesIn.end(); ++frm) {
    double theta = data_->Dval( *frm ) * Constants::DEGRAD;
    sumy += sin( theta );
    sumx += cos( theta );
  }
  ((Centroid_Num*)centIn)->SetCval( atan2(sumy, sumx) * Constants::RADDEG );
  ((Centroid_Num*)centIn)->SetPeriodicSums( sumx, sumy );
}

/** \return New centroid from specified frames. */
Cpptraj::Cluster::Centroid*
  Cpptraj::Cluster::Metric_Torsion::NewCentroid(Cframes const& cframesIn)
{
  Centroid_Num* cent = new Centroid_Num();
  CalculateCentroid(cent, cframesIn);
  return cent;
}

/** Perform given operation between frame and centroid. */
void Cpptraj::Cluster::Metric_Torsion::FrameOpCentroid(int frame, Centroid* centIn,
                                                       double oldSize, CentOpType OP)
{
  Centroid_Num* cent = (Centroid_Num*)centIn;

  double sumx = cent->SumX();
  double sumy = cent->SumY();
  double ftheta = data_->Dval(frame) * Constants::DEGRAD;
  if (OP == ADDFRAME) {
    sumy += sin( ftheta );
    sumx += cos( ftheta );
  } else { // SUBTRACTFRAME
    sumy -= sin( ftheta );
    sumx -= cos( ftheta );
  }
  double newcval = atan2(sumy, sumx) * Constants::RADDEG;

  cent->SetCval( newcval );
  cent->SetPeriodicSums( sumx, sumy );
}

/** \return 1 line description */
std::string Cpptraj::Cluster::Metric_Torsion::Description() const {
  return "Torsion data set " + data_->Meta().Legend();
}

/** Print info to STDOUT. */
void Cpptraj::Cluster::Metric_Torsion::Info() const {
  mprintf("\tMetric: Torsion set type: '%s'\n", data_->description());
}

/** \return DataSet size */
unsigned int Cpptraj::Cluster::Metric_Torsion::Ntotal() const { return data_->Size(); }
