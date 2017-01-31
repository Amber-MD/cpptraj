#ifndef INC_EWALD_H
#define INC_EWALD_H
#include "Topology.h"
#include "Timer.h"
#include "PairList.h"
#include "Spline.h"
class Ewald {
  public:
    Ewald();
    int EwaldInit(Box const&, double, double, double, double, double, double,
                  double, int, const int*);
    void EwaldSetup(Topology const&, AtomMask const&);
    double CalcEnergy_NoPairList(Frame const&, Topology const&, AtomMask const&);
    double CalcEnergy(Frame const&, AtomMask const&);
    void Timing(double) const;
  private:
    static double erfc_func(double);
    static double FindEwaldCoefficient(double,double);
    static double FindMaxexpFromMlim(const int*, Matrix_3x3 const&);
    static double FindMaxexpFromTol(double, double);
    static void GetMlimits(int*, double, double, Vec3 const&, Matrix_3x3 const&);
    void FillErfcTable(double,double);
    inline double ERFC(double) const;

    double Self(double);
    double Recip_Regular(Matrix_3x3 const&, double);
    double Direct(Matrix_3x3 const&, Topology const&, AtomMask const&);
    double Direct(PairList const&);
    inline void Adjust(double,double,double);

    typedef std::vector<double> Darray;
    typedef std::vector<Vec3> Varray;

    Varray Cells_;  ///< Hold fractional translations to neighbor cells (non-pairlist only)
    Darray Charge_; ///< Hold atomic charges converted to Amber units.
    Darray cosf1_;
    Darray cosf2_;
    Darray cosf3_;
    Darray sinf1_;
    Darray sinf2_;
    Darray sinf3_;

    PairList pairList_;

    Spline cspline_;
    Darray erfc_table_Y_;

    typedef std::vector< std::set<int> > Iarray2D;
    Iarray2D Excluded_; ///< Full exclusion list for each atom.

    static const double INVSQRTPI_;
    double sumq_; ///< Sum of charges
    double sumq2_; ///< Sum of charges squared
    double ew_coeff_; ///< Ewald coefficient
    double maxexp_;
    double cutoff_; ///< Direct space cutoff
    double dsumTol_; ///< Direct space sum tolerance.
    double rsumTol_; ///< Reciprocal space sum tolerance.
    double e_adjust_; ///< Adjustment for excluded pairs.
    double erfcTableDx_;
    double one_over_Dx_;
    int mlimit_[3];
    int maxmlim_;
    Timer t_total_;
    Timer t_self_;
    Timer t_recip_;
    Timer t_trig_tables_;
    Timer t_direct_;
    Timer t_erfc_;
    Timer t_adjust_;
};
#endif
