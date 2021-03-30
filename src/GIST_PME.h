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
    using Ewald::Iarray;
    using Ewald_ParticleMesh::Init;
    using Ewald_ParticleMesh::Setup;

    /// Calculate nonbonded energy with PME for GIST
    int CalcNonbondEnergy_GIST(Frame const&, AtomMask const&,
                                  double&, double&,
                                  Darray&, Darray&, Darray&, Darray&, Darray&, Darray&, Darray&,
                                  Iarray const&);
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
    double Direct_GIST(PairList const&, double&, Darray&, Darray&, Iarray const&);
    /// Calcualte direct space energy with long range VDW correction for GIST, decomposed for every atom.
    double Direct_VDW_LongRangeCorrection_GIST(PairList const&, double&, Iarray const&);
    /// Calculate direct space energy with LJ PME for GIST, decomposed for every atom.
    double Direct_VDW_LJPME_GIST(PairList const&, double&, Darray&, Darray&);

    std::vector<Darray> E_vdw_direct_;  ///< Only for interaction of water molecules involved in direct (sw, ww)
    std::vector<Darray> E_elec_direct_; ///< Count all interactions in direct (sw, ww, ss)
};
#endif /* LIBPME */
#endif
