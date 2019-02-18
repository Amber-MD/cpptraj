#ifndef INC_CLUSTER_CENTROID_MULTI_H
#define INC_CLUSTER_CENTROID_MULTI_H
#include <vector>
namespace Cpptraj {
namespace Cluster {

/// Cluster centroid for multiple DataSets
class Centroid_Multi : public Centroid {
  public:
    typedef std::vector<double> Darray;
    Centroid_Multi() {}
    Centroid_Multi(Darray const& val, Darray const& x, Darray const& y) :
      cvals_(val), Sumx_(x), Sumy_(y) {}

    Centroid* Copy() { return (Centroid*)new Centroid_Multi(cvals_, Sumx_, Sumy_); }

    Darray const& Cvals() const { return cvals_; }

    Darray& Cvals() { return cvals_; }
    Darray& SumX()  { return Sumx_;  }
    Darray& SumY()  { return Sumy_;  }
  private:
    Darray cvals_;
    Darray Sumx_; // For storing periodic average
    Darray Sumy_; // For storing periodic average
};


}
}
#endif
