#ifndef INC_GIST_PME_H
#define INC_GIST_PME_H
#ifdef LIBPME
#include <vector>
#include "AtomMask.h"
#include "Energy/EwaldParams.h"
#include "Energy/PME_RecipParams.h"
#include "Energy/VDW_LongRange_Correction.h"
#include "ExclusionArray.h"
#include "helpme_standalone.h"
#include "PairList.h"
#include "Timer.h"
class Box;
class Frame;
class Topology;
/// Class implementing the PME version of the nonbonded energy calc. for GIST
/** For more debug info, compile with:
  *  -DDEBUG_GIST_PME : Details on breakdown of PME calculation for the UV/VV grids.
  *  -DDEBUG_PAIRLIST : Details on the use of the pair list in the direct space calc.
  */
class GIST_PME : private Cpptraj::Energy::EwaldParams {
  public:
    GIST_PME();

    // Expose definitions/functions from Ewald_ParticleMesh
    using Cpptraj::Energy::EwaldParams::Darray;
    //using Ewald_ParticleMesh::Init;
    //using Ewald::Timing;

    typedef std::vector<float> Farray;

    int Init(Box const&, EwaldOptions const&, int);
    /// Setup PME calc. for top, all atoms. Allocate memory for internal arrays (# threads)
    int Setup_PME_GIST(Topology const&, unsigned int, double);
    /// Print timing data
    void Timing(double) const;
    /// Calculate nonbonded energy with PME for GIST
    int CalcNonbondEnergy_GIST(Frame const&, std::vector<int> const&,
                               std::vector<bool> const&,
                               std::vector<bool> const&,
                               // std::vector<Darray>&, std::vector<Darray>&,
                               std::vector<Darray>&, std::vector<Darray>&,
                               std::vector<Farray>&);
    /// \return Total energy of specified atom
    double E_of_atom(unsigned int idx) const {
      return (E_elec_self_[idx] + E_elec_direct_[0][idx] + E_elec_recip_[idx] +
              E_vdw_self_[idx]  + E_vdw_direct_[0][idx] +
              E_vdw_recip_[idx] + E_vdw_lr_cor_[idx]);
    }

    // Internal arrays
/*
    Darray const& E_Vdw_Direct()  const { return E_vdw_direct_[0]; }
    Darray const& E_Elec_Direct() const { return E_elec_direct_[0]; }

    Darray const& E_Vdw_Self()    const { return E_vdw_self_; }
    Darray const& E_Vdw_Recip()   const { return E_vdw_recip_; }
    Darray const& E_Vdw_LR_Corr() const { return E_vdw_lr_cor_; }

    Darray const& E_Elec_Self()   const { return E_elec_self_; }
    Darray const& E_Elec_Recip()  const { return E_elec_recip_; }
*/
  private:
    /// Atom to atom interaction type
    enum InteractionType { OTHER = 0,        ///< Interaction we do not care about
                           SOLUTE0_ONGRID1,  ///< Solute 0, on-grid solvent 1
                           SOLVENT0_ONGRID1, ///< Off-grid solvent 0, on-grid solvent 1
                           SOLUTE1_ONGRID0,  ///< Solute 1, on-grid solvent 0
                           SOLVENT1_ONGRID0, ///< Off-grid solvent 1, on-grid solvent 0
                           BOTH_ONGRID       ///< Both 0 and 1 on-grid solvent
                         };
    /// Strings corresponding to InteractionType
    static const char* InteractionTypeStr_[];

    typedef helpme::Matrix<double> MatType;

    /// Determine the interaction type between two atoms given grid voxels and solute status
    static inline InteractionType determineInteractionType(int, bool, int, bool);

    /// Nonbond energy kernel
    inline void Ekernel_NB(double&, double&, double, double, double, int, int, double*, double*,
                           InteractionType, double*, double*) const;
    // Adjust energy kernel
    inline void Ekernel_Adjust(double&, double, double, double, int, int, double*,
                               InteractionType, double*, double*);

    /// Electrostatic self energy, decomposed onto atoms.
    double Self_GIST(double, Darray&, std::vector<int> const&, std::vector<bool> const&, Darray&);
    /// Lennard-Jones self energy, decomposed onto atoms.
    //double Self6_GIST(Darray&);
    /// Reciprocal energy decomposed for every atom.
    double Recip_ParticleMesh_GIST(Box const&, std::vector<int> const&, std::vector<bool> const&,
                                   Darray&, Darray&);
    /// LJ reciprocal term, decomposed for every atom.
    //double LJ_Recip_ParticleMesh_GIST(Box const&, MatType&);
    /// VDW long range correction for GIST
    double Vdw_Correction_GIST(double, std::vector<bool> const&,
                               Darray&, Darray&);
    /// Calculate direct space energy for GIST, decomposed for every atom.
    double Direct_GIST(PairList const&, double&, std::vector<int> const&,
                       std::vector<bool> const&,
                       std::vector<bool> const&,
                       // std::vector<Darray>&, std::vector<Darray>&,
                       std::vector<Darray>&, std::vector<Darray>&,
                       std::vector<Farray>&);
    /// Calcualte direct space energy with long range VDW correction for GIST, decomposed for every atom.
    double Direct_VDW_LongRangeCorrection_GIST(PairList const&, double&, std::vector<int> const&,
                                               std::vector<bool> const&,
                                               std::vector<bool> const&,
                                               // std::vector<Darray>&, std::vector<Darray>&,
                                               std::vector<Darray>&, std::vector<Darray>&,
                                               std::vector<Farray>&);
    /// Calculate direct space energy with LJ PME for GIST, decomposed for every atom.
    //double Direct_VDW_LJPME_GIST(PairList const&, double&);

    // TODO could potentially make the calculation more efficient by skipping
    //      the per-atom arrays below and just passing in the GIST PME
    //      energy grids and summing into them directly.

    std::vector<Darray> E_vdw_direct_;  ///< VDW direct energy for each atom (sw, ww)
    std::vector<Darray> E_elec_direct_; ///< Elec. direct energy for each atom (sw, ww, ss)

    Darray E_vdw_self_;   ///< VDW self energy for each atom (LJ PME)
    Darray E_vdw_recip_;  ///< VDW recip energy for each atom (LJ PME)
    Darray E_vdw_lr_cor_; ///< VDW long range correction energy for each atom

    Darray E_elec_self_;  ///< Elec. self energy for each atom
    Darray E_elec_recip_; ///< Elec. recip energy for each atom

    MatType e_potentialD_; ///< Hold recip contributions for each atom

    AtomMask allAtoms_;    ///< Select all atoms

    double NeighborCut2_; ///< Cutoff for the O-O neighbor calculation
    Timer t_total_;
    Timer t_recip_;
    Timer t_direct_;
    Timer t_self_;

    PairList pairList_;
    Darray coordsD_;
    Cpptraj::Energy::PME_RecipParams recipParams_; ///< Hold parameters for recip part of PME
    PMEInstanceD pme_object_;
    Cpptraj::Energy::VDW_LongRange_Correction VDW_LR_; ///< For calculating the long range VDW correction
    ExclusionArray excluded_;
};
#endif /* LIBPME */
#endif
