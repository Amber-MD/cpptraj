#ifndef INC_ENERGY_ECALC_NONBOND_H
#define INC_ENERGY_ECALC_NONBOND_H
#include "../ExclusionArray.h"
#include "../PairList.h"
class CharMask;
class EwaldOptions;
class Topology;
namespace Cpptraj {
namespace Energy {
class EwaldCalc;
/// Calculate nonbonded energy for atoms
class Ecalc_Nonbond {
  public:
    typedef std::vector<double> Darray;

    enum CalcType { SIMPLE = 0, PME, LJPME, REGULAR_EWALD, UNSPECIFIED };
    /// CONSTRUCTOR
    Ecalc_Nonbond();
    /// DESTRUCTOR
    ~Ecalc_Nonbond();
    /// Init - Calc type, decomp?, box, options, debug
    int InitNonbondCalc(CalcType, bool, Box const&, EwaldOptions const&, int);
    /// Setup
    int SetupNonbondCalc(Topology const&, AtomMask const&);
    /// Calculate energy
    int NonbondEnergy(Frame const&, AtomMask const&, double&, double&);
    /// Calculate decomposed energy
    int DecomposedNonbondEnergy(Frame const&, CharMask const&, double&, double&,
                                Darray&, Darray&);
    /// Print timing to stdout
    void PrintTiming(double) const;
  private:
    EwaldCalc* calc_;
    Topology const* currentTop_; ///< Current Topology, set by SetupNonbondCalc()
    PairList pairList_;
    ExclusionArray Excluded_;
    Timer t_total_;
    CalcType type_;
    bool decompose_energy_;
};
}
}
#endif
