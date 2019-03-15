#ifndef INC_EWALD_H
#define INC_EWALD_H
#include "Topology.h"
#include "Timer.h"
#include "PairList.h"
/// Base class for calculating electrostatics using Ewald methods.
class Ewald {
  public:
    Ewald();
    // ----- Virtual functions -------------------
    virtual ~Ewald() {}
    virtual int Setup(Topology const&, AtomMask const&) = 0;
    virtual double CalcEnergy(Frame const&, AtomMask const&, double&) = 0; // TODO const?
    // -------------------------------------------
    /// Report timings.
    void Timing(double) const;
#   ifdef DEBUG_EWALD
    /// Slow non-pairlist version of energy calc. For debug only.
    double CalcEnergy_NoPairList(Frame const&, Topology const&, AtomMask const&);
#   endif
  protected:
    typedef std::vector<double> Darray;
    typedef std::vector<int> Iarray;
    typedef std::vector<Vec3> Varray;
    typedef std::set<int> Iset;
    typedef std::vector<Iset> Iarray2D;

    static inline double DABS(double xIn) { if (xIn < 0.0) return -xIn; else return xIn; }
    /// Complimentary error function, erfc.
    static double erfc_func(double);
    /// Determine Ewald coefficient from cutoff and direct sum tolerance.
    static double FindEwaldCoefficient(double,double);

    /// Fill erfc lookup table using cubic spline interpolation.
    void FillErfcTable(double,double);
    /// \return erfc value from erfc lookup table.
    inline double ERFC(double) const;
    /// Ewald "self" energy
    double Self(double);
    /// Ewald "self" energy for C6 term
    double Self6();
    /// Get analytical estimate of energy due to dispersion interactions > cutoff
    double Vdw_Correction(double);
    /// Box, debug, cutoff, dsum tol, ew coeff, lj coeff, switch window, erfc dx, nb skin
    int CheckInput(Box const&, int, double, double, double, double, double, double, double);
    /// Set up pair list
    int Setup_Pairlist(Box const&, Vec3 const&, double);
    /// Calculate sum q, sum q^2. Calls setup for vdw correction
    void CalculateCharges(Topology const&, AtomMask const&);
    /// Calculate VDW C6 parameters for LJ PME
    void CalculateC6params(Topology const&, AtomMask const&);
    /// Setup main excluded atom list
    void SetupExcluded(Topology const&, AtomMask const&);
    /// Setup VDW correction for selected atom types
    void Setup_VDW_Correction(Topology const&, AtomMask const&);

#   ifdef DEBUG_EWALD
    /// Slow version of direct space energy, no pairlist.
    double Direct(Matrix_3x3 const&, Topology const&, AtomMask const&);
#   endif
    /// Fast version of direct space energy using a pairlist
    double Direct(PairList const&, double&);
    /// \return adjusted energy for excluded atom pair
#   ifdef _OPENMP
    inline double Adjust(double,double,double) const;
#   else
    inline double Adjust(double,double,double); // Cannot be const bc timers
#   endif

    // TODO make variables private
#   ifdef DEBUG_EWALD
    Varray Cells_;  ///< Hold fractional translations to neighbor cells (non-pairlist only)
#   endif
    Darray Charge_;       ///< Hold selected atomic charges converted to Amber units.
    Darray Cparam_;       ///< Hold selected atomic C6 coefficients for LJ PME
    PairList pairList_;   ///< Atom pair list for direct sum.
    Darray erfc_table_;   ///< Hold Erfc cubic spline Y values and coefficients (Y B C D).
    Iarray2D Excluded_;   ///< Full exclusion list for each selected atom.
    Iarray TypeIndices_;  ///< Hold atom type indices for selected atoms
    NonbondParmType const* NB_; ///< Pointer to nonbonded parameters

    static const double INVSQRTPI_;
    double sumq_;         ///< Sum of charges
    double sumq2_;        ///< Sum of charges squared
    double ew_coeff_;     ///< Ewald coefficient for electrostatics
    double lw_coeff_;     ///< Ewald coefficient for LJ
    double switch_width_; ///< Switching window size for LJ switch if active
    double cutoff_;       ///< Direct space cutoff
    double cut2_;         ///< Direct space cutoff squared.
    double cut2_0_;       ///< Direct space cutoff minus switch width, squared.
    double dsumTol_;      ///< Direct space sum tolerance.
    double erfcTableDx_;  ///< Spacing of X values in Erfc table.
    double one_over_Dx_;  ///< One over erfcTableDx_.
    double Vdw_Recip_term_; ///< VDW recip correction term from # types and B parameters
    int debug_;
    Timer t_total_; // TODO make timing external
    Timer t_self_;
    Timer t_recip_;
    Timer t_trig_tables_;
    Timer t_direct_;
    Timer t_erfc_;
    Timer t_adjust_;
  private:
    double Direct_VDW_LongRangeCorrection(PairList const&, double&);
    double Direct_VDW_LJPME(PairList const&, double&);
};
#endif
