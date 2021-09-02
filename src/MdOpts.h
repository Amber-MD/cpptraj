#ifndef INC_MDOPTS_H
#define INC_MDOPTS_H
#include "Constraints.h"
class MdOpts {
  public:
    MdOpts();

    int GetOptsFromArgs(ArgList&);

    void PrintOpts() const;

    Constraints::ShakeType Shake() const { return shakeType_; }
    double ScaleEE()       const { return scaleEE_; }
    double ScaleNB()       const { return scaleNB_; }
    double CutEE()         const { return cutEE_; }
    double CutNB()         const { return cutNB_; }
    double CoulombFactor() const { return qfac_; }
    int N_Exclude()        const { return nExclude_; }
    /// Print options recognized by GetOptsFromArgs() to stdout.
    static void PrintHelp();
  private:
    Constraints::ShakeType shakeType_;
    double scaleEE_; ///< Global electrostatic 1-4 scaling factor
    double scaleNB_; ///< Global Lennard-Jones (VDW) 1-4 scaling factor
    double cutEE_;   ///< Nonbond electrostatic cutoff
    double cutNB_;   ///< Nonbond Lennard-Jones (VDW) cutoff
    double qfac_;    ///< Factor for getting kcal/mol units in Coulomb calculation (i.e. qfac*qi*qj)
    int nExclude_;   ///< Exclude nonbond interactions within this many bonded atoms
};
#endif
