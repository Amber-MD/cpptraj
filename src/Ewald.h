#ifndef INC_EWALD_H
#define INC_EWALD_H
#include "Topology.h"
class Ewald {
  public:
    Ewald();
    int SetupParams(double, double, double);
    void CalcSumQ(Topology const&, AtomMask const&);
    double CalcEnergy(Frame const&, Topology const&, AtomMask const&);
    double Self(double);
  private:
    static double erfc_func(double);
    static double FindEwaldCoefficient(double,double);

    static double INVSQRTPI_;
    double sumq_; ///< Sum of charges
    double sumq2_; ///< Sum of charges squared
    double ew_coeff_; ///< Ewald coefficient
    double cutoff_; ///< Direct space cutoff
    double dsumTol_; ///< Direct space sum tolerance.
    bool needSumQ_; ///< True if sum over charges needs to be calcd. (TODO)
};
#endif
