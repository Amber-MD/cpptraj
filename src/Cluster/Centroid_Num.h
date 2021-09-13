#ifndef INC_CLUSTER_CENTROID_NUM_H
#define INC_CLUSTER_CENTROID_NUM_H
#include "Centroid.h"
namespace Cpptraj {
namespace Cluster {

/// Cluster centroid for generic DataSet
class Centroid_Num : public Centroid {
  public:
    Centroid_Num()                               : cval_(0.0), sumx_(0.0), sumy_(0.0) {}
    Centroid_Num(double val, double x, double y) : cval_(val), sumx_(x), sumy_(y) {}

    Centroid* Copy() { return (Centroid*)new Centroid_Num(cval_, sumx_, sumy_); }
    /// \return Current value of average.
    double Cval() const { return cval_; }

    /// Set average value.
    void SetCval(double c) { cval_ = c; }
  private:
    double cval_;
    double sumx_; // For storing periodic average
    double sumy_; // For storing periodic average
};

}
}
#endif
