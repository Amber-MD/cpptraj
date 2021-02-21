#ifndef INC_RANDOM_H
#define INC_RANDOM_H
namespace Cpptraj {
  class RNG;
}
/// Wrapper around RNG
class Random_Number {
  public:
    Random_Number();
    ~Random_Number();
    /// Different possible RNGs
    enum RngType { MARSAGLIAS=0, STDLIB, MERSENNE_TWISTER };

    static void SetDefaultRng(RngType);

    static void SetDefaultSeed(int);

    /// Allocate and initialize the random number generator with the given seed
    int rn_set(int);
    /// Initialize RN generator with 71277 (Amber default)
    //void rn_set();
    /// Generate a random number between 0.0 and 1.0
    double rn_gen() const;
    /// Generate a pseudo-random Gaussian sequence.
    double rn_gauss(double,double) const;
    /// \return true if RN generator has been set up.
    bool IsSet() const;
    /// \return Value of RN generator seed
    int Seed() const;
  private:
    /// Allocate the RNG
    void allocateRng();

    Cpptraj::RNG* rng_;          ///< Hold the random number generator implemenatation.

    static RngType defaultType_; ///< Default RNG type
    static int defaultSeed_;     ///< Default RNG seed
};
#endif
