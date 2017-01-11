#ifndef INC_EWALD_H
#define INC_EWALD_H
#include "Topology.h"
class Ewald {
  public:
    Ewald();
    void CalcSumQ(Topology const&, AtomMask const&);
    void FindEwaldCoefficient(double,double);
  private:
    double erfc_func(double) const;

    static double INVSQRTPI_;
    double sumq_;
    double sumq2_;
    double ew_coeff_;
};
#endif
