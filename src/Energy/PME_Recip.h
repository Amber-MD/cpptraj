#ifndef INC_ENERGY_PME_RECIP_H
#define INC_ENERGY_PME_RECIP_H
#ifdef LIBPME
#include "../helpme_standalone.h"
#include "../Timer.h"
#include "PME_RecipParams.h"
class Box;
namespace Cpptraj {
namespace Energy {
/// Do the reciprocal part of a PME calculation
class PME_Recip {
    typedef std::vector<double> Darray;
  public:
    enum Type { COULOMB = 0, LJ };

    PME_Recip(Type);

    /// Initialize
    int InitRecip(EwaldOptions const& pmeOpts, int);

    double Recip_ParticleMesh(Darray&, Box const&, Darray&, double);
    double Recip_Decomp(Darray&, Darray&, Box const&, Darray&, double);

    Timer const& Timing_Total() const { return t_recip_; }
    //Timer const& Timing_Calc() const { return t_calc_; }
  private:
    /// Print options to stdout
    void PrintRecipOpts() const;

    static int set_lattice(PMEInstanceD::LatticeType&, Box const&);

    PMEInstanceD pme_object_;
    PME_RecipParams recipParams_; ///< Hold PME recip parameters
    Timer t_recip_;               ///< Recip calc timer
    //Timer t_calc_;
    int distKernelExponent_;      ///< Exponent of the distance kernel: 1 for Coulomb, 6 for LJ
    double scaleFac_;             ///< scale factor to be applied to all computed energies and derivatives thereof
};
}
}
#endif /* LIBPME */
#endif
