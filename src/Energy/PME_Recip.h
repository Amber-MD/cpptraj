#ifndef INC_ENERGY_PME_RECIP_H
#define INC_ENERGY_PME_RECIP_H
#include "../helpme_standalone.h"
#include "../Timer.h"
class Box;
namespace Cpptraj {
namespace Energy {
class PME_RecipParams;
/// Do the reciprocal part of a PME calculation
class PME_Recip {
    typedef std::vector<double> Darray;
  public:
    enum Type { COULOMB = 0, LJ };

    PME_Recip(Type);

    double Recip_ParticleMesh(PME_RecipParams const&, Darray&, Box const&, Darray&, double);
    double Recip_Decomp(PME_RecipParams const&, Darray&, Darray&, Box const&, Darray&, double);

    Timer const& Timing_Total() const { return t_recip_; }
    //Timer const& Timing_Calc() const { return t_calc_; }
  private:
    static int set_lattice(PMEInstanceD::LatticeType&, Box const&);

    PMEInstanceD pme_object_;
    Timer t_recip_; ///< Recip calc timer
    //Timer t_calc_;
    int distKernelExponent_; ///< Exponent of the distance kernel: 1 for Coulomb, 6 for LJ
    double scaleFac_; ///< scale factor to be applied to all computed energies and derivatives thereof
};
}
}
#endif
