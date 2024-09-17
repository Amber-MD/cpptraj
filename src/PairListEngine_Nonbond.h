#ifndef INC_PAIRLISTENGINE_NONBOND_H
#define INC_PAIRLISTENGINE_NONBOND_H
#include "Energy/Ene_LJ_6_12.h"
#include "ErfcFxn.h"
#include "Kernel_EwaldAdjust.h"
#include "Kernel_LJswitch.h"
#include "PairList.h"
#include "ParameterTypes.h"
#include <vector>
namespace Cpptraj {
/// Nonbond calculation for pairlist with Ewald and VDW LR correction
template <typename T>
class PairListEngine_Nonbond {
    typedef std::vector<int> Iarray;
  public:
    PairListEngine_Nonbond() : NB_(0) {}
    /// Call for atom 0 when looping over atoms of thisCell
    void SetupAtom0( PairList::AtmType const& atom0 ) {
      q0_ = Charge_[atom0.Idx()];
    }
    /// Call for atom 1 when looping over interaction atoms of this/other cell
    void SetupAtom1( PairList::AtmType const& atom1 ) {
      q1_ = Charge_[atom1.Idx()];
    }
    /// Call at the beginning of the frame calculation
    void FrameBeginCalc() { Evdw_ = 0; Eelec_ = 0; Eadjust_ = 0; }
    /// Call when cutoff is satisfied
    void CutoffSatisfied(T const& rij2,
                         PairList::AtmType const& atom0,
                         PairList::AtmType const& atom1)
    {
      double rij = sqrt( rij2 );
      double qiqj = q0_ * q1_;
      //double erfc = erfc_func(ew_coeff_ * rij);
      double erfc = erfc_.ERFC(ew_coeff_ * rij);
      double e_elec = qiqj * erfc / rij;
      Eelec_ += e_elec;

      int nbindex = NB_->GetLJindex(TypeIndices_[atom0.Idx()],
                                    TypeIndices_[atom1.Idx()]);
      if (nbindex > -1) {
        double vswitch = ljswitch_.switch_fn(rij2);
        NonbondType const& LJ = NB_->NBarray()[ nbindex ];
        T e_vdw = Cpptraj::Energy::Ene_LJ_6_12<T>(rij2, LJ.A(), LJ.B());
        Evdw_ += (e_vdw * vswitch);
      }
    }
    /// Call when cutoff is not satisfied
    void CutoffNotSatisfied(T const& rij2,
                            PairList::AtmType const& atom0,
                            PairList::AtmType const& atom1)
    {
      Eadjust_ += Kernel_EwaldAdjust<T>( q0_, q1_, sqrt(rij2), ew_coeff_, erfc_ );
    }
  private:
    
    T q0_;                  ///< Charge on atom 0
    T q1_;                  ///< Charge on atom 1
    T Evdw_;                ///< VDW sum for current frame
    T Eelec_;               ///< Coulomb sum for current frame
    T Eadjust_;             ///< Adjust energy sum for current frame

    T ew_coeff_;            ///< Ewald coefficient

    ErfcFxn erfc_; ///< Hold spline interpolation for erfc
    std::vector<T> Charge_; ///< Array of charges
    Iarray TypeIndices_;  ///< Hold atom type indices for selected atoms
    NonbondParmType const* NB_; ///< Pointer to nonbonded parameters
    Kernel_LJswitch<T> ljswitch_; ///< LJ switching function
};
} // END namespace Cpptraj
#endif
