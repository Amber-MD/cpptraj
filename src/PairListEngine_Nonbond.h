#ifndef INC_PAIRLISTENGINE_NONBOND_H
#define INC_PAIRLISTENGINE_NONBOND_H
#include "Kernel_LJswitch.h"
#include "PairList.h"
#include "ParameterTypes.h"
#include "SplineFxnTable.h"
#include <vector>
template <typename T>
class PairListEngine_Nonbond {
    typedef std::vector<int> Iarray;
  public:
    PairListEngine_Nonbond() {}
    /// Call for atom 0 when looping over atoms of thisCell
    void SetupAtom0( PairList::AtmType const& atom0 ) {
      q0_ = Charge_[atom0.Idx()];
    }
    /// Call for atom 1 when looping over interaction atoms of this/other cell
    void SetupAtom1( PairList::AtmType const& atom1 ) {
      q1_ = Charge_[atom1.Idx()];
    }
    /// Call when cutoff is satisfied
    void CutoffSatisfied(T const& rij2,
                         PairList::AtmType const& atom0,
                         PairList::AtmType const& atom1)
    {
      double rij = sqrt( rij2 );
      double qiqj = q0_ * q1_;
#     ifndef _OPENMP
      t_erfc_.Start();
#     endif
      //double erfc = erfc_func(ew_coeff_ * rij);
      double erfc = ERFC(ew_coeff_ * rij);
#     ifndef _OPENMP
      t_erfc_.Stop();
#     endif
      double e_elec = qiqj * erfc / rij;
      Eelec += e_elec;
      //mprintf("EELEC %4i%4i%12.5f%12.5f%12.5f%3.0f%3.0f%3.0f\n",
      //int ta0, ta1;
      //if (it0->Idx() < it1->Idx()) {
      //  ta0=it0->Idx(); ta1=it1->Idx();
      //} else {
      //  ta1=it0->Idx(); ta0=it1->Idx();
      //}
      //mprintf("PELEC %6i%6i%12.5f%12.5f%12.5f\n", ta0, ta1, rij, erfc, e_elec);
      int nbindex = NB_->GetLJindex(TypeIndices_[atom0.Idx()],
                                    TypeIndices_[atom1.Idx()]);
      if (nbindex > -1) {
        double vswitch = ljswitch.switch_fn(rij2);
        NonbondType const& LJ = NB_->NBarray()[ nbindex ];
        double r2    = 1.0 / rij2;
        double r6    = r2 * r2 * r2;
        double r12   = r6 * r6;
        double f12   = LJ.A() * r12;  // A/r^12
        double f6    = LJ.B() * r6;   // B/r^6
        double e_vdw = f12 - f6;      // (A/r^12)-(B/r^6)
        Evdw += (e_vdw * vswitch);
        //mprintf("PVDW %8i%8i%20.6f%20.6f\n", ta0+1, ta1+1, e_vdw, r2);
      }
    }
  private:
    double ERFC(double xIn) const {
      return table_.Yval( xIn);
    }

    T q0_;                  ///< Charge on atom 0
    T q1_;                  ///< Charge on atom 1

    SplineFxnTable table_; ///< Hold spline interpolation for erfc
    std::vector<T> Charge_; ///< Array of charges
    Iarray TypeIndices_;  ///< Hold atom type indices for selected atoms
    NonbondParmType const* NB_; ///< Pointer to nonbonded parameters
    Kernel_LJswitch<T> ljswitch_; ///< LJ switching function
};
#endif
