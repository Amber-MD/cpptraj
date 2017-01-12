#ifndef INC_EWALD_H
#define INC_EWALD_H
#include "Topology.h"
class Ewald {
  public:
    Ewald();
    int SetupParams(Box const&, double, double, double, double, double, const int*);
    void CalcSumQ(Topology const&, AtomMask const&);
    double CalcEnergy(Frame const&, Topology const&, AtomMask const&);
    double Self(double);
  private:
    static double erfc_func(double);
    static double FindEwaldCoefficient(double,double);
    static double FindMaxexpFromMlim(const int*, Matrix_3x3 const&);
    static double FindMaxexpFromTol(double, double);
    static void GetMlimits(int*, double, double, Vec3 const&, Matrix_3x3 const&);

    Matrix_3x3 ucell_;
    Matrix_3x3 recip_;

    static double INVSQRTPI_;
    double sumq_; ///< Sum of charges
    double sumq2_; ///< Sum of charges squared
    double ew_coeff_; ///< Ewald coefficient
    double maxexp_;
    double cutoff_; ///< Direct space cutoff
    double dsumTol_; ///< Direct space sum tolerance.
    double rsumTol_; ///< Reciprocal space sum tolerance.
    int mlimit_[3];
    bool needSumQ_; ///< True if sum over charges needs to be calcd. (TODO)
};
#endif
