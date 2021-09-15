#include "Metric_Scalar.h"
#include "Centroid_Num.h"
#include "Cframes.h"
#include "../CpptrajStdio.h"
#include "../DataSet_1D.h"
#include <cmath> // fabs

/** CONSTRUCTOR */
Cpptraj::Cluster::Metric_Scalar::Metric_Scalar() :
  Metric(SCALAR),
  data_(0)
{}

/** Initialize. */
int Cpptraj::Cluster::Metric_Scalar::Init(DataSet_1D* dataIn) {
  if (dataIn == 0) {
    mprinterr("Internal Error: Metric_Scalar::Init called with null set.\n");
    return 1;
  }
  data_ = dataIn;
  return 0;
}

/** Set up. Not much to do for a DataSet. */
int Cpptraj::Cluster::Metric_Scalar::Setup() {
  return 0;
}

/** \return Absolute difference between points */
double Cpptraj::Cluster::Metric_Scalar::FrameDist(int f1, int f2) {
  return fabs( data_->Dval(f1) - data_->Dval(f2) );
}

/** \return Absolute difference between centroids. */
double Cpptraj::Cluster::Metric_Scalar::CentroidDist(Centroid* c1, Centroid* c2) {
  return fabs( ((Centroid_Num*)c1)->Cval() - ((Centroid_Num*)c2)->Cval() );
}

/** \return Absolute difference between point and centroid. */
double  Cpptraj::Cluster::Metric_Scalar::FrameCentroidDist(int f1, Centroid* c1) {
  return fabs( data_->Dval(f1) - ((Centroid_Num*)c1)->Cval() );
}

/** Calculate centroid from specified frames. */
void Cpptraj::Cluster::Metric_Scalar::CalculateCentroid(Centroid* centIn, Cframes const& cframesIn)
{
  double val = 0.0;
  for (Cframes::const_iterator frm = cframesIn.begin(); frm != cframesIn.end(); ++frm)
    val += data_->Dval( *frm );
  ((Centroid_Num*)centIn)->SetCval( val / (double)cframesIn.size() );
}

/** \return New centroid from specified frames. */
Cpptraj::Cluster::Centroid*
  Cpptraj::Cluster::Metric_Scalar::NewCentroid(Cframes const& cframesIn)
{
  Centroid_Num* cent = new Centroid_Num();
  CalculateCentroid(cent, cframesIn);
  return cent;
}

/** Perform given operation between frame and centroid. */
void Cpptraj::Cluster::Metric_Scalar::FrameOpCentroid(int frame, Centroid* centIn,
                                                      double oldSize, CentOpType OP)
{
  Centroid_Num* cent = (Centroid_Num*)centIn;

  double newcval = cent->Cval() * oldSize;
  double fval = data_->Dval(frame);
  if (OP == ADDFRAME) {
    newcval += fval;
    newcval /= ( oldSize + 1 );
  } else { // SUBTRACTFRAME
    newcval -= fval;
    newcval /= ( oldSize - 1 );
  }

  cent->SetCval( newcval );
}

/** \return 1 line description */
std::string Cpptraj::Cluster::Metric_Scalar::Description() const {
  return "Data set " + data_->Meta().Legend();
}

/** Print info to STDOUT. */
void Cpptraj::Cluster::Metric_Scalar::Info() const {
  mprintf("\tMetric: Data set type: '%s'\n", data_->description());
}

/** \return DataSet size */
unsigned int Cpptraj::Cluster::Metric_Scalar::Ntotal() const { return data_->Size(); }
