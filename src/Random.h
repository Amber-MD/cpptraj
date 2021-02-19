#ifndef INC_RANDOM_H
#define INC_RANDOM_H
// Class: Random_Number
namespace Cpptraj {
class RNG;
}
/// Wrapper around RNG
class Random_Number {
  public:
    Random_Number();
    ~Random_Number();
    /// Different possible RNGs
    enum RngType { MARSAGLIAS=0, STDLIB };

    static void SetDefaultRng(RngType);

    /// Initialize the random number generator with the given seed
    void rn_set(int);
    /// Initialize RN generator with 71277 (Amber default)
    void rn_set();
    /// Generate a random number between 0.0 and 1.0
    double rn_gen();
    /// Generate a pseudo-random Gaussian sequence.
    double rn_gauss(double,double);
    /// \return true if RN generator has been set up.
    bool IsSet() const;
    /// \return Value of RN generator seed
    int Seed() const;
  private:
    void allocateRng();

    Cpptraj::RNG* rng_; ///< Hold the random number generator.

    static RngType default_; ///< Default RNG type
};
#endif
