#ifndef INC_EWALD_H
#define INC_EWALD_H
#include "Topology.h"
#include "Timer.h"
#include "PairList.h"
class Ewald {
  public:
    Ewald();
    int EwaldInit(Box const&, double, double, double, double, double, const int*);
    void EwaldSetup(Topology const&, AtomMask const&);
    double CalcEnergy(Frame const&, Topology const&, AtomMask const&);
    void Timing(double) const;
  private:
    static double erfc_func(double);
    static double FindEwaldCoefficient(double,double);
    static double FindMaxexpFromMlim(const int*, Matrix_3x3 const&);
    static double FindMaxexpFromTol(double, double);
    static void GetMlimits(int*, double, double, Vec3 const&, Matrix_3x3 const&);

    void MapCoords(Frame const&, Matrix_3x3 const&,Matrix_3x3 const&, AtomMask const&);
    double Self(double);
    double Recip_Regular(Matrix_3x3 const&, double);
    double Direct(Matrix_3x3 const&, Topology const&, AtomMask const&);

    typedef std::vector<double> Darray;
    typedef std::vector<Vec3> Varray;

    Varray Cells_;  ///< Hold fractional translations to neighbor cells.
    Varray Frac_;   ///< Hold fractional coords back in primary cell.
    Varray Image_;  ///< Hold Cartesian coords back in primary cell.
    Darray Charge_; ///< Hold atomic charges converted to Amber units.
    Darray cosf1_;
    Darray cosf2_;
    Darray cosf3_;
    Darray sinf1_;
    Darray sinf2_;
    Darray sinf3_;

    PairList pairList_;

    static const double INVSQRTPI_;
    double sumq_; ///< Sum of charges
    double sumq2_; ///< Sum of charges squared
    double ew_coeff_; ///< Ewald coefficient
    double maxexp_;
    double cutoff_; ///< Direct space cutoff
    double dsumTol_; ///< Direct space sum tolerance.
    double rsumTol_; ///< Reciprocal space sum tolerance.
    int mlimit_[3];
    int maxmlim_;
    bool needSumQ_; ///< True if sum over charges needs to be calcd. (TODO)
    Timer t_total_;
    Timer t_map_;
    Timer t_self_;
    Timer t_recip_;
    Timer t_direct_;
    static Timer t_erfc_; // DEBUG
};
#endif
