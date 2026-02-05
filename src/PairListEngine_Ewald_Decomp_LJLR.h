#ifndef INC_PAIRLISTENGINE_EWALD_DECOMP_LJLR_H
#define INC_PAIRLISTENGINE_EWALD_DECOMP_LJLR_H
#include "Energy/Ene_LJ_6_12.h"
#include "Energy/EwaldParams.h"
#include "Energy/Kernel_EwaldAdjust.h"
#include "PairList.h"
namespace Cpptraj {
/// Direct space nonbond calculation using pairlist with decomposable Ewald and VDW LR correction
template <typename T>
class PairListEngine_Ewald_Decomp_LJLR {
    typedef std::vector<T> Darray;
  public:
    PairListEngine_Ewald_Decomp_LJLR() {}
    // -------------------------------------------
    /// Call at the beginning of the frame calculation
    void FrameBeginCalc() { Evdw_ = 0; Eelec_ = 0; Eadjust_ = 0; 
                            atom_elec_.assign(EW_.Natom(), 0);
                            atom_evdw_.assign(EW_.Natom(), 0);
                            atom_eadj_.assign(EW_.Natom(), 0);
                          }
    /// Call for atom 0 when looping over atoms of thisCell
    void SetupAtom0( PairList::AtmType const& atom0 ) {
      q0_ = EW_.Charge(atom0.Idx());
    }
    /// Call for atom 1 when looping over interaction atoms of this/other cell
    void SetupAtom1( PairList::AtmType const& atom1 ) {
      q1_ = EW_.Charge(atom1.Idx());
    }
    /// Call when cutoff is satisfied
    void CutoffSatisfied(T const& rij2,
                         PairList::AtmType const& atom0,
                         PairList::AtmType const& atom1)
    {
      T rij = sqrt( rij2 );
      T qiqj = q0_ * q1_;
      //double erfc = erfc_func(ew_coeff_ * rij);
      T erfcval = EW_.ErfcEW( rij );
      T e_elec = qiqj * erfcval / rij;
      Eelec_ += e_elec;
      T e_half = e_elec * 0.5;
      atom_elec_[atom0.Idx()] += e_half;
      atom_elec_[atom1.Idx()] += e_half;

      int nbindex = EW_.NbIndex(atom0.Idx(), atom1.Idx());
      if (nbindex > -1) {
        double vswitch = EW_.Switch_Fn(rij2);
        NonbondType const& LJ = EW_.GetLJ( nbindex );
        T e_vdw = Cpptraj::Energy::Ene_LJ_6_12<T>(rij2, LJ.A(), LJ.B());
        e_vdw *= vswitch;
        Evdw_ += e_vdw;
        e_half = e_vdw * 0.5;
        atom_evdw_[atom0.Idx()] += e_half;
        atom_evdw_[atom1.Idx()] += e_half;
      }
    }
    /// Call when cutoff is not satisfied
    void AtomPairExcluded(T const& rij2,
                            PairList::AtmType const& atom0,
                            PairList::AtmType const& atom1)
    {
      T rij = sqrt(rij2);
      T erfcval = EW_.ErfcEW( rij );
      T e_adj = Cpptraj::Energy::Kernel_EwaldAdjust<T>( q0_, q1_, rij, erfcval );
      Eadjust_ += e_adj;
      T e_half = e_adj * 0.5;
      atom_eadj_[atom0.Idx()] += e_half;
      atom_eadj_[atom1.Idx()] += e_half;
    }
    // -------------------------------------------
    Cpptraj::Energy::EwaldParams& ModifyEwaldParams() { return EW_; }
    Cpptraj::Energy::EwaldParams const& EwaldParams() const { return EW_; }

    T Evdw() const { return Evdw_; }
    T Eelec() const { return Eelec_; }
    T Eadjust() const { return Eadjust_; }
    Darray const& Eatom_Elec() const { return atom_elec_; }
    Darray const& Eatom_EVDW() const { return atom_evdw_; }
    Darray const& Eatom_EAdjust() const { return atom_eadj_; }
#   ifdef _OPENMP
    static void sum_Darray(Darray& lhs, Darray const& rhs) {
      for (unsigned int idx = 0; idx != rhs.size(); ++idx)
        lhs[idx] += rhs[idx];
    }
    /// To allow reduction of the energy terms
    void operator+=(PairListEngine_Ewald_Decomp_LJLR const& rhs) {
      Evdw_ += rhs.Evdw_;
      Eelec_ += rhs.Eelec_;
      Eadjust_ += rhs.Eadjust_;
      sum_Darray(atom_elec_, rhs.atom_elec_);
      sum_Darray(atom_evdw_, rhs.atom_evdw_);
      sum_Darray(atom_eadj_, rhs.atom_eadj_);
    }
#   endif
  private:
    T q0_;                  ///< Charge on atom 0
    T q1_;                  ///< Charge on atom 1
    T Evdw_;                ///< VDW sum for current frame
    T Eelec_;               ///< Coulomb sum for current frame
    T Eadjust_;             ///< Adjust energy sum for current frame

    Darray atom_elec_;      ///< Sum of Coulomb E on each atom for current frame
    Darray atom_evdw_;      ///< Sum of VDW E on each atom for current frame
    Darray atom_eadj_;      ///< Sum of excluded atom adjust E on each atom for current frame

    Cpptraj::Energy::EwaldParams EW_;          ///< Hold Ewald parameters for PME
};
#ifdef _OPENMP
#pragma omp declare reduction( + : PairListEngine_Ewald_Decomp_LJLR<double> : omp_out += omp_in ) initializer( omp_priv = omp_orig )
#endif
} // END namespace Cpptraj
#endif
