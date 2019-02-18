#include <cmath> // fabs 
#include "Metric_Data_Manhattan.h"
#include "Centroid_Multi.h"
#include "../CpptrajStdio.h"

double Cpptraj::Cluster::Metric_Data_Manhattan::FrameDist(int f1, int f2) {
  double dist = 0.0;
  DcArray::const_iterator dcalc = dcalcs_.begin();
  for (D1Array::iterator ds = dsets_.begin(); ds != dsets_.end(); ++ds, ++dcalc) {
    dist += fabs( (*dcalc)((*ds)->Dval(f1), (*ds)->Dval(f2)) );
  }
  return dist;
}

double Cpptraj::Cluster::Metric_Data_Manhattan::CentroidDist(Centroid* c1, Centroid* c2) {
  double dist = 0.0;
  Centroid_Multi::Darray::const_iterator c2val = ((Centroid_Multi*)c2)->Cvals().begin();
  DcArray::const_iterator dcalc = dcalcs_.begin();
  for (Centroid_Multi::Darray::const_iterator c1val = ((Centroid_Multi*)c1)->Cvals().begin();
                                              c1val != ((Centroid_Multi*)c1)->Cvals().end();
                                            ++c1val, ++dcalc)
  {
    dist += fabs( (*dcalc)(*c1val, *(c2val++)) );
  }
  return dist;
}

double Cpptraj::Cluster::Metric_Data_Manhattan::FrameCentroidDist(int f1, Centroid* c1) {
  double dist = 0.0;
  Centroid_Multi::Darray::const_iterator c1val = ((Centroid_Multi*)c1)->Cvals().begin();
  DcArray::const_iterator dcalc = dcalcs_.begin();
  for (D1Array::const_iterator ds = dsets_.begin(); ds != dsets_.end(); ++ds) {
    dist += fabs( (*dcalc)((*ds)->Dval(f1), *(c1val++)) );
  }
  return dist;
}

std::string Cpptraj::Cluster::Metric_Data_Manhattan::Description() const {
  return SetNames("data (Manhattan) ");
}

void Cpptraj::Cluster::Metric_Data_Manhattan::Info() const {
  mprintf("\tMetric: Manhattan distance (%zu sets)\n", dsets_.size());
}
