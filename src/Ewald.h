#ifndef INC_EWALD_H
#define INC_EWALD_H
class Topology;
class AtomMask;
class Frame;
class NonbondParmType;
class EwaldOptions;
#include "Timer.h"
#include "PairList.h"
#include "ExclusionArray.h"
#include "SplineFxnTable.h"
/// Base class for calculating non-bonded energy using Ewald methods.
class Ewald {
  public:
    Ewald();
    // ----- Virtual functions -------------------
    virtual ~Ewald() {}
    virtual int Init(Box const&, EwaldOptions const&, int) = 0;
    virtual int Setup(Topology const&, AtomMask const&) = 0;
    /// Calculate electrostatic and van der Waals energy
    virtual int CalcNonbondEnergy(Frame const&, AtomMask const&, double&, double&) = 0;
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

    static inline double DABS(double xIn) { if (xIn < 0.0) return -xIn; else return xIn; }
    /// Complimentary error function, erfc.
    static double erfc_func(double);

    /// Ewald "self" energy
    double Self(double);
    /// Ewald "self" energy for C6 term
    double Self6();
    /// Get analytical estimate of energy due to dispersion interactions > cutoff
    double Vdw_Correction(double);
    /// Box, debug, cutoff, dsum tol, ew coeff, lj coeff, switch window, erfc dx, nb skin
    int CheckInput(Box const&, int, double, double, double, double, double, double, double);
    /// Set up pair list for given box and NB "skin" size
    int Setup_Pairlist(Box const&, double);
    /// Calculate sum q, sum q^2. Calls setup for vdw correction
    void CalculateCharges(Topology const&, AtomMask const&);
    /// Calculate VDW C6 parameters for LJ PME
    void CalculateC6params(Topology const&, AtomMask const&);
    /// Setup main excluded atom list
    void SetupExclusionList(Topology const&, AtomMask const&);

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

    /// \return sum of charges squared
    double SumQ2() const { return sumq2_; }
    /// \return sum of charges
    double SumQ()  const { return sumq_; }
    /// \return VDW recip correction term from # types and B parameters
    double Vdw_Recip_Term() const { return Vdw_Recip_term_; }
    /// \return Atom exclusion array
    ExclusionArray const& Excluded() const { return Excluded_; }
    /// \return Value of Erfc at given value
    double ErfcFxn(double) const;
    /// \return Nonbond parameters
    NonbondParmType const& NB() const { return *NB_; }
    /// \return Type index for given atom
    int TypeIdx(unsigned int idx) const { return TypeIndices_[idx]; }
    /// \return Value of LJ switching function
    static double SwitchFxn(double, double, double);
    /// \return Value of Ewald adjustment
#   ifdef _OPENMP
    double AdjustFxn(double,double,double) const;
#   else
    double AdjustFxn(double,double,double);
#   endif

    // TODO make variables private
    Darray Charge_;       ///< Hold selected atomic charges converted to Amber units.
    Darray Cparam_;       ///< Hold selected atomic C6 coefficients for LJ PME
    PairList pairList_;   ///< Atom pair list for direct sum.

    Iarray vdw_type_;              ///< Store nonbond vdw type for each atom
    Iarray N_vdw_type_;            ///< Total number of atoms for each vdw type
    Darray atype_vdw_recip_terms_; ///< Nonbond PME interaction correction for each vdw type

    static const double INVSQRTPI_;
    double ew_coeff_;     ///< Ewald coefficient for electrostatics
    double lw_coeff_;     ///< Ewald coefficient for LJ
    double switch_width_; ///< Switching window size for LJ switch if active
    double cutoff_;       ///< Direct space cutoff
    double cut2_;         ///< Direct space cutoff squared.
    double cut2_0_;       ///< Direct space cutoff minus switch width, squared.
    double dsumTol_;      ///< Direct space sum tolerance.
    int debug_;
    Timer t_total_; // TODO make timing external
    Timer t_self_;
    Timer t_recip_;
    Timer t_trig_tables_;
    Timer t_direct_;
    Timer t_erfc_;
    Timer t_adjust_;
  private:
    /// \return erfc value from erfc lookup table.
    inline double ERFC(double) const;
    /// Determine Ewald coefficient from cutoff and direct sum tolerance.
    static double FindEwaldCoefficient(double,double);

    /// Setup VDW correction for selected atom types
    void Setup_VDW_Correction(Topology const&, AtomMask const&);
    /// Direct-space energy with VDW long range corrected energy
    double Direct_VDW_LongRangeCorrection(PairList const&, double&);
    /// Direct-space energy with VDW handled via PME
    double Direct_VDW_LJPME(PairList const&, double&);

    SplineFxnTable table_; ///< Hold spline interpolation for erfc
#   ifdef DEBUG_EWALD
    Varray Cells_;  ///< Hold fractional translations to neighbor cells (non-pairlist only)
#   endif
    double sumq_;         ///< Sum of charges
    double sumq2_;        ///< Sum of charges squared
    double Vdw_Recip_term_; ///< VDW recip correction term from # types and B parameters
    // TODO should Exlcusions be passed in?
    ExclusionArray Excluded_;   ///< Full exclusion list for each selected atom.
    Iarray TypeIndices_;  ///< Hold atom type indices for selected atoms
    NonbondParmType const* NB_; ///< Pointer to nonbonded parameters

};
#endif
