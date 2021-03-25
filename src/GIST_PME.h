#ifndef INC_GIST_PME_H
#define INC_GIST_PME_H
#ifdef LIBPME
#include <vector>
#include "Ewald_ParticleMesh.h"
class Frame;
class AtomMask;
class Box;
/// Class implementing the PME version of the nonbonded energy calc. for GIST
class GIST_PME : private Ewald_ParticleMesh {
  public:
    GIST_PME();

    typedef std::vector<double> Darray;
    typedef std::vector<int> Iarray;

    /// Calculate nonbonded energy with PME for GIST
    int CalcNonbondEnergy_GIST(Frame const&, AtomMask const&,
                                  double&, double&,
                                  Darray&, Darray&, Darray&, Darray&, Darray&, Darray&, Darray&,
                                  Iarray const&);
  private:
    /// Electrostatic self energy, decomposed onto atoms.
    double Self_GIST(double, Darray&);
    /// Lennard-Jones self energy, decomposed onto atoms.
    double Self6_GIST(Darray&);
    /// Reciprocal energy decomposed for every atom.
    double Recip_ParticleMesh_GIST(Box const&, helpme::Matrix<double>&);
    /// LJ reciprocal term, decomposed for every atom.
    double LJ_Recip_ParticleMesh_GIST(Box const&, helpme::Matrix<double>&);
    /// VDW long range correction for GIST
    double Vdw_Correction_GIST(double, Darray&);

};
#endif /* LIBPME */
#endif
