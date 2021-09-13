#include "Metric_Scalar.h"
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
