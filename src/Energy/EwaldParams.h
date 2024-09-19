#ifndef INC_ENERGY_EWALDPARAMS_H
#define INC_ENERGY_EWALDPARAMS_H
#include "ErfcFxn.h"
#include "../ParameterTypes.h" // NonbondParmType
class AtomMask;
class Box;
class EwaldOptions;
class Topology;
namespace Cpptraj {
namespace Energy {
class EwaldParams {
  public:
    EwaldParams();
    // Virtual since inherited
    virtual ~EwaldParams() {}
    // -------------------------------------------
    virtual int InitEwald(Box const&, EwaldOptions const&, int) = 0;
    virtual int SetupEwald(Topology const&, AtomMask const&) = 0;
    // -------------------------------------------

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

    /// \return Direct space cutoff (in Ang squared)
    double Cut2() const { return cut2_; }
    /// \return Charge for given atom index
    double Charge(int idx) { return Charge_[idx]; }
    /// \return Nonbonded index for given atom indices
    int NbIndex(int idx0, int idx1) { return NB_->GetLJindex(TypeIndices_[idx0], TypeIndices_[idx1]); }
    /// \return Nonbonded parameter at nonobonded parameter index
    NonbondType const& GetLJ(int nbindex) const { return NB_->NBarray()[ nbindex ]; }

    /// \return Direct space cutoff (in Ang)
    double Cutoff() const { return cutoff_; }
    /// \return Direct sum tolerance
    double DirectSumTol() const { return dsumTol_; }
    /// \return Ewald coefficient
    double EwaldCoeff() const { return ew_coeff_; }
    /// \return LJ Ewald coefficient
    double LJ_EwaldCoeff() const { return lw_coeff_; }
    /// \return LJ switch width (in Ang.)
    double LJ_SwitchWidth() const { return switch_width_; }
  protected:
    typedef std::vector<int> Iarray;
    typedef std::vector<double> Darray;

    /// Set Ewald parametsr, check them and set defaults if needed.
    int CheckInput(Box const&, int, double, double,
                   double, double, double,
                   double, double);
    /// Calculate sum q, sum q^2. 
    void CalculateCharges(Topology const&, AtomMask const&);
  private:

    double FindEwaldCoefficient(double, double);

    double ew_coeff_;                  ///< Ewald coefficient
    double lw_coeff_;                  ///< LJ Ewald coefficient
    double switch_width_; ///< Switching window size for LJ switch if active
    double cutoff_;       ///< Direct space cutoff
    double cut2_;         ///< Direct space cutoff squared.
    double cut2_0_;       ///< Direct space cutoff minus switch width, squared.
    double dsumTol_;      ///< Direct space sum tolerance.
    int debug_;

    ErfcFxn erfc_;              ///< Hold spline interpolation for erfc

    Darray Charge_;             ///< Hold charges for selected atoms
    Iarray TypeIndices_;        ///< Hold atom type indices for selected atoms
    NonbondParmType const* NB_; ///< Pointer to nonbonded parameters
    double sumq_;         ///< Sum of charges
    double sumq2_;        ///< Sum of charges squared
};
}
}
#endif
