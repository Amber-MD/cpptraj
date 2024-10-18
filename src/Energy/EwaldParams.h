#ifndef INC_ENERGY_EWALDPARAMS_H
#define INC_ENERGY_EWALDPARAMS_H
#include "ErfcFxn.h"
#include "../ParameterTypes.h" // NonbondParmType
class AtomMask;
class Box;
class EwaldOptions;
class Frame;
class Topology;
namespace Cpptraj {
namespace Energy {
/// Parameters common to all Ewald methods
class EwaldParams {
    static const double INVSQRTPI_;
  public:
    EwaldParams();
    // Virtual since inherited
    virtual ~EwaldParams() {}
    // -------------------------------------------
    virtual int InitEwald(Box const&, EwaldOptions const&, int);
    virtual int SetupEwald(Topology const&, AtomMask const&);
    // -------------------------------------------

    ///  Fill recip coords with XYZ coords of selected atoms.
    void FillRecipCoords(Frame const&, AtomMask const&);
    /// \return ERFC value of given distance times the Ewald coefficient
    double ErfcEW(double rIn) const { return erfc_.ErfcInterpolated( ew_coeff_*rIn ); }
    /// \return LJ switch fn value
    double Switch_Fn(double rij2) const {
      double cut2_1 = cut2_;
      if (rij2 <= cut2_0_)
        return 1.0;
      else if (rij2 > cut2_1)
        return 0.0;
      else {
        double xoff_m_x = cut2_1 - rij2;
        double fac = 1.0 / (cut2_1 - cut2_0_);
        return (xoff_m_x*xoff_m_x) * (cut2_1 + 2.0*rij2 - 3.0*cut2_0_) * (fac*fac*fac);
      }
    }
    /// \return Self energy for the given volume
    double SelfEnergy(double) const;
    /// \return Self energy for given volume. Set self energy for each atom.
    double DecomposedSelfEnergy(std::vector<double>&, double) const;

    /// \return Direct space cutoff (in Ang squared)
    double Cut2() const { return cut2_; }
    /// \return Charge for given atom index
    double Charge(int idx) const { return Charge_[idx]; }
    /// \return Nonbonded index for given atom indices
    int NbIndex(int idx0, int idx1) const { return NB_->GetLJindex(TypeIndices_[idx0], TypeIndices_[idx1]); }
    /// \return Nonbonded parameter at nonobonded parameter index
    NonbondType const& GetLJ(int nbindex) const { return NB_->NBarray()[ nbindex ]; }
    /// \return Number of atoms (current size of charge array)
    unsigned int Natom() const { return Charge_.size(); }

    /// \return Debug level
    int Debug() const { return debug_; }
    /// \return Direct space cutoff (in Ang)
    double Cutoff() const { return cutoff_; }
    /// \return Direct sum tolerance
    double DirectSumTol() const { return dsumTol_; }
    /// \return Ewald coefficient
    double EwaldCoeff() const { return ew_coeff_; }
    /// \return LJ switch width (in Ang.)
    double LJ_SwitchWidth() const { return switch_width_; }
    /// \return 1 / sqrt(PI)
    static const double INVSQRTPI() { return INVSQRTPI_; }
    /// \return Charge array
    std::vector<double> const& Charge() const { return Charge_; }

    // FIXME do not return const because helPME needs the array to be non-const. Should be fixed
    std::vector<double>& SelectedCharges() { return Charge_; }
    // FIXME do not return const because helPME needs the array to be non-const. Should be fixed
    std::vector<double>& SelectedCoords() { return coordsD_; }

    // FIXME these can probably go away after GIST_PME is converted
    double SumQ() const { return sumq_; }
    double SumQ2() const { return sumq2_; }
  protected:
    typedef std::vector<double> Darray;
    typedef std::vector<int> Iarray;

    static inline double DABS(double xIn) { if (xIn < 0.0) return -xIn; else return xIn; }
    /// Set Ewald parametsr, check them and set defaults if needed.
    int CheckInput(Box const&, int, double, double,
                   double, double, double, double);
  private:
    double FindEwaldCoefficient(double, double);
    /// Reserve space for selected coordinates for recip calcs
    void reserveRecipCoords(AtomMask const&);

    double ew_coeff_;                  ///< Ewald coefficient
    double switch_width_; ///< Switching window size for LJ switch if active
    double cutoff_;       ///< Direct space cutoff
    double cut2_;         ///< Direct space cutoff squared.
    double cut2_0_;       ///< Direct space cutoff minus switch width, squared.
    double dsumTol_;      ///< Direct space sum tolerance.
    int debug_;

    ErfcFxn erfc_;              ///< Hold spline interpolation for erfc

    Darray coordsD_;            ///< Hold selected coords for recip calcs
    Darray Charge_;             ///< Hold charges for selected atoms
    Iarray TypeIndices_;        ///< Hold atom type indices for selected atoms
    NonbondParmType const* NB_; ///< Pointer to nonbonded parameters
    double sumq_;         ///< Sum of charges
    double sumq2_;        ///< Sum of charges squared
};
}
}
#endif
