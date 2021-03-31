#ifndef INC_GIST_PME_H
#define INC_GIST_PME_H
#ifdef LIBPME
#include <vector>
#include "Ewald_ParticleMesh.h"
class Frame;
class AtomMask;
class Box;
class Topology;
/// Class implementing the PME version of the nonbonded energy calc. for GIST
class GIST_PME : private Ewald_ParticleMesh {
  public:
    GIST_PME();

    // Expose definitions/functions from Ewald_ParticleMesh
    using Ewald::Darray;
    using Ewald_ParticleMesh::Init;
    using Ewald_ParticleMesh::Setup;
    using Ewald::Timing;

    /// Allocate memory for internal arrays (# atoms, # threads)
    int AllocateArrays(unsigned int, unsigned int);
    /// Calculate nonbonded energy with PME for GIST
    int CalcNonbondEnergy_GIST(Frame const&, AtomMask const&,
                                  double&, double&);
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
    typedef helpme::Matrix<double> MatType;

    /// Electrostatic self energy, decomposed onto atoms.
    double Self_GIST(double, Darray&);
    /// Lennard-Jones self energy, decomposed onto atoms.
    double Self6_GIST(Darray&);
    /// Reciprocal energy decomposed for every atom.
    double Recip_ParticleMesh_GIST(Box const&, MatType&);
    /// LJ reciprocal term, decomposed for every atom.
    double LJ_Recip_ParticleMesh_GIST(Box const&, MatType&);
    /// VDW long range correction for GIST
    double Vdw_Correction_GIST(double, Darray&);
    /// Calculate direct space energy for GIST, decomposed for every atom.
    double Direct_GIST(PairList const&, double&);
    /// Calcualte direct space energy with long range VDW correction for GIST, decomposed for every atom.
    double Direct_VDW_LongRangeCorrection_GIST(PairList const&, double&);
    /// Calculate direct space energy with LJ PME for GIST, decomposed for every atom.
    double Direct_VDW_LJPME_GIST(PairList const&, double&);

    std::vector<Darray> E_vdw_direct_;  ///< VDW direct energy for each atom (sw, ww)
    std::vector<Darray> E_elec_direct_; ///< Elec. direct energy for each atom (sw, ww, ss)

    Darray E_vdw_self_;   ///< VDW self energy for each atom (LJ PME)
    Darray E_vdw_recip_;  ///< VDW recip energy for each atom (LJ PME)
    Darray E_vdw_lr_cor_; ///< VDW long range correction energy for each atom

    Darray E_elec_self_;  ///< Elec. self energy for each atom
    Darray E_elec_recip_; ///< Elec. recip energy for each atom

    MatType e_potentialD_; ///< Hold recip contributions for each atom
};
#endif /* LIBPME */
#endif
