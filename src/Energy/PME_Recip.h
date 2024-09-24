#ifndef INC_ENERGY_PME_RECIP_H
#define INC_ENERGY_PME_RECIP_H
#include "../helpme_standalone.h"
#include "../Timer.h"
class AtomMask;
class Box;
class Frame;
namespace Cpptraj {
namespace Energy {
/// Do the reciprocal part of a PME calculation
class PME_Recip {
    typedef std::vector<double> Darray;
  public:
    enum Type { COULOMB = 0, LJ };

    PME_Recip(Type);
    void SetDebug(int);
    double Recip_ParticleMesh(Darray&, Box const&, Darray&, const int*, double, int);

    Timer const& Timing_Total() const { return t_recip_; }
    //Timer const& Timing_Calc() const { return t_calc_; }
  private:
    static bool check_prime_factors(int);
    static int ComputeNFFT(double);
    int DetermineNfft(int&, int&, int&, Box const&) const;

    PMEInstanceD pme_object_;
    Timer t_recip_; ///< Recip calc timer
    //Timer t_calc_;
    int debug_;
    int distKernelExponent_; ///< Exponent of the distance kernel: 1 for Coulomb, 6 for LJ
    double scaleFac_; ///< scale factor to be applied to all computed energies and derivatives thereof
};
}
}
#endif
