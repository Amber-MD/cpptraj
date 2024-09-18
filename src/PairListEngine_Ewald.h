#ifndef INC_PAIRLISTENGINE_EWALD_H
#define INC_PAIRLISTENGINE_EWALD_H
#include "ErfcFxn.h"
#include "Kernel_LJswitch.h"
#include "PairList.h"
#include <vector>
class NonbondParmType;
namespace Cpptraj {
/// Base class for direct space nonbond calculation using pairlist with Ewald
class PairListEngine_Ewald : public PairListEngine {
    typedef std::vector<int> Iarray;
  public:
    PairListEngine_Ewald() : NB_(0) {}
    // -------------------------------------------
    /// Call at the beginning of the frame calculation
    void FrameBeginCalc() { Evdw_ = 0; Eelec_ = 0; Eadjust_ = 0; }
    /// Call for atom 0 when looping over atoms of thisCell
    void SetupAtom0( PairList::AtmType const& atom0 ) {
      q0_ = Charge_[atom0.Idx()];
    }
    /// Call for atom 1 when looping over interaction atoms of this/other cell
    void SetupAtom1( PairList::AtmType const& atom1 ) {
      q1_ = Charge_[atom1.Idx()];
    }
    // -------------------------------------------
    int CheckInput(Box const&, int, double, double,
                   double, double, double,
                   double, double);

  private:
    double q0_;                  ///< Charge on atom 0
    double q1_;                  ///< Charge on atom 1
    double Evdw_;                ///< VDW sum for current frame
    double Eelec_;               ///< Coulomb sum for current frame
    double Eadjust_;             ///< Adjust energy sum for current frame

    double ew_coeff_;                  ///< Ewald coefficient
    ErfcFxn erfc_;                     ///< Hold spline interpolation for erfc
    std::vector<double> Charge_;       ///< Array of charges
    Iarray TypeIndices_;               ///< Hold atom type indices for selected atoms
    NonbondParmType const* NB_;        ///< Pointer to nonbonded parameters
    Kernel_LJswitch<double> ljswitch_; ///< LJ switching function
};
} // END namespace Cpptraj
#endif
