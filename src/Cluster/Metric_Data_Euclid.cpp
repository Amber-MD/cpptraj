#include <cmath> // sqrt
#include "Metric_Data_Euclid.h"
#include "Centroid_Multi.h"

double Cpptraj::Cluster::Metric_Data_Euclid::FrameDist(int f1, int f2) {
  double dist = 0.0;
  DcArray::const_iterator dcalc = dcalcs_.begin();
  for (D1Array::iterator ds = dsets_.begin(); ds != dsets_.end(); ++ds, ++dcalc) {
    double diff = (*dcalc)((*ds)->Dval(f1), (*ds)->Dval(f2));
    dist += (diff * diff);
  }
  return sqrt(dist);
}

double Cpptraj::Cluster::Metric_Data_Euclid::CentroidDist(Centroid* c1, Centroid* c2) {
  double dist = 0.0;
  Centroid_Multi::Darray::const_iterator c2val = ((Centroid_Multi*)c2)->Cvals().begin();
  DcArray::const_iterator dcalc = dcalcs_.begin();
  for (Centroid_Multi::Darray::const_iterator c1val = ((Centroid_Multi*)c1)->Cvals().begin();
                                              c1val != ((Centroid_Multi*)c1)->Cvals().end();
                                            ++c1val, ++dcalc)
  {
    double diff = (*dcalc)(*c1val, *(c2val++));
    dist += (diff * diff);
  }
  return sqrt(dist);
}

double Cpptraj::Cluster::Metric_Data_Euclid::FrameCentroidDist(int f1, Centroid* c1) {
  double dist = 0.0;
  Centroid_Multi::Darray::const_iterator c1val = ((Centroid_Multi*)c1)->Cvals().begin();
  DcArray::const_iterator dcalc = dcalcs_.begin();
  for (D1Array::const_iterator ds = dsets_.begin(); ds != dsets_.end(); ++ds) {
    double diff = (*dcalc)((*ds)->Dval(f1), *(c1val++));
    dist += (diff * diff);
  }
  return sqrt(dist);
}

std::string Cpptraj::Cluster::Metric_Data_Euclid::Description() const {
  return SetNames("data (Euclidean) ");
}
