#ifndef INC_GIST_PME_H
#define INC_GIST_PME_H
#ifdef LIBPME
#include <vector>
#include "Ewald_ParticleMesh.h"
class Frame;
class AtomMask;
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
    double Self_GIST(double, Darray&) const;
};
#endif /* LIBPME */
#endif
