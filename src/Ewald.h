#ifndef INC_EWALD_H
#define INC_EWALD_H
#include "Topology.h"
#include "Timer.h"
#include "PairList.h"
#ifdef LIBPME
#  include <memory> // unique_ptr (libpme_standalone.h)
#  include "libpme_standalone.h"
#endif
/// Class for calculating electrostatics using Ewald summation.
class Ewald {
  public:
    Ewald();
    // ----- Virtual functions -------------------
    virtual ~Ewald() {}
    virtual int Setup(Topology const&, AtomMask const&) = 0;
    virtual double CalcEnergy(Frame const&, AtomMask const&) = 0; // TODO const?
    // -------------------------------------------
    /// Report timings.
    void Timing(double) const;
#   ifdef DEBUG_EWALD
    /// Slow non-pairlist version of energy calc. For debug only.
    double CalcEnergy_NoPairList(Frame const&, Topology const&, AtomMask const&);
#   endif
  protected:
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

    int CheckInput(Box const&, int, double, double, double, double, double);
    int Setup_Pairlist(Box const&, Vec3 const&, double);
    void CalculateCharges(Topology const&, AtomMask const&);
    void SetupExcluded(Topology const&); // TODO fix for atom mask

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
    typedef std::vector<double> Darray;
    typedef std::vector<Vec3> Varray;
    typedef std::set<int> Iset;
    typedef std::vector<Iset> Iarray2D;
#   ifdef DEBUG_EWALD
    Varray Cells_;  ///< Hold fractional translations to neighbor cells (non-pairlist only)
#   endif
    Darray Charge_; ///< Hold atomic charges converted to Amber units.
    PairList pairList_;   ///< Atom pair list for direct sum.
    Darray erfc_table_;   ///< Hold Erfc cubic spline Y values and coefficients (Y B C D).
    Iarray2D Excluded_;   ///< Full exclusion list for each atom.

    static const double INVSQRTPI_;
    double sumq_;         ///< Sum of charges
    double sumq2_;        ///< Sum of charges squared
    double ew_coeff_;     ///< Ewald coefficient
    double cutoff_;       ///< Direct space cutoff
    double dsumTol_;      ///< Direct space sum tolerance.
    double erfcTableDx_;  ///< Spacing of X values in Erfc table.
    double one_over_Dx_;  ///< One over erfcTableDx_.
    int debug_;
    Timer t_total_;
    Timer t_self_;
    Timer t_recip_;
    Timer t_trig_tables_;
    Timer t_direct_;
    Timer t_erfc_;
    Timer t_adjust_;
};
#endif
