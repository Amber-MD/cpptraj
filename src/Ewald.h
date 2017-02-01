#ifndef INC_EWALD_H
#define INC_EWALD_H
#include "Topology.h"
#include "Timer.h"
#include "PairList.h"
/// Class for calculating electrostatics using Ewald summation.
class Ewald {
  public:
    Ewald();
    /// Initialize Ewald parameters.
    int EwaldInit(Box const&, double, double, double, double, double, double,
                  double, int, const int*);
    /// Set up for given topology and mask.
    void EwaldSetup(Topology const&, AtomMask const&);
    /// Calculate electrostatic energy via Ewald summation.
    double CalcEnergy(Frame const&, AtomMask const&);
    /// Report timings.
    void Timing(double) const;
    /// Slow non-pairlist version of energy calc. For debug only.
    double CalcEnergy_NoPairList(Frame const&, Topology const&, AtomMask const&);
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
    inline double Adjust(double,double,double) const;

    typedef std::vector<double> Darray;
    typedef std::vector<Vec3> Varray;
    typedef std::vector< std::set<int> > Iarray2D;

    Varray Cells_;  ///< Hold fractional translations to neighbor cells (non-pairlist only)
    Darray Charge_; ///< Hold atomic charges converted to Amber units.
    Darray cosf1_;
    Darray cosf2_;
    Darray cosf3_;
    Darray sinf1_;
    Darray sinf2_;
    Darray sinf3_;
    Darray c12_;
    Darray s12_;
    Darray c3_;
    Darray s3_;
    PairList pairList_;   ///< Atom pair list for direct sum.
    Darray erfc_table_;   ///< Hold Erfc cubic spline Y values and coefficients (Y B C D).
    Iarray2D Excluded_;   ///< Full exclusion list for each atom.
    static const double INVSQRTPI_;
    double sumq_;         ///< Sum of charges
    double sumq2_;        ///< Sum of charges squared
    double ew_coeff_;     ///< Ewald coefficient
    double maxexp_;       ///< Determines how far out recip vectors go? FIXME check!
    double cutoff_;       ///< Direct space cutoff
    double dsumTol_;      ///< Direct space sum tolerance.
    double rsumTol_;      ///< Reciprocal space sum tolerance.
    double e_adjust_;     ///< Adjustment for excluded pairs, calcd during direct sum.
    double erfcTableDx_;  ///< Spacing of X values in Erfc table.
    double one_over_Dx_;  ///< One over erfcTableDx_.
    int mlimit_[3];       ///< Number of units in each direction to calc recip. sum.
    int maxmlim_;         ///< The max of the three mlimit_ values.
    Timer t_total_;
    Timer t_self_;
    Timer t_recip_;
    Timer t_trig_tables_;
    Timer t_direct_;
    Timer t_erfc_;
    Timer t_adjust_;
};
#endif
