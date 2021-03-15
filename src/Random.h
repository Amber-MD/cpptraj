#ifndef INC_RANDOM_H
#define INC_RANDOM_H
#include <vector>
namespace Cpptraj {
  class RNG;
}
/// Wrapper around RNG
class Random_Number {
  public:
    Random_Number();
    ~Random_Number();
    /// Different possible RNGs
    enum RngType { MARSAGLIA = 0, STDLIB, MERSENNE_TWISTER, PCG32, XOSHIRO128PP };

    static void SetDefaultRng(RngType);

    //static void SetDefaultSeed(int);

    static const char* CurrentDefaultRngStr();

    /// Allocate and initialize the random number generator with the given seed
    int rn_set(int);
    /// Initialize RN generator with 71277 (Amber default)
    //void rn_set();
    /// Generate a random number between 0.0 and 1.0
    double rn_gen() const;
    /// Generate a random integer
    unsigned int rn_num() const;
    /// Generate a random number on defined interval
    int rn_num_interval_signed(int, int) const;
    /// Generate a random number on defined interval
    unsigned int rn_num_interval(unsigned int, unsigned int) const;
    /// Generate random numbers in Gaussian distribution with given mean and SD.
    double GenerateGauss(double, double) const;
    /// Generate a pseudo-random Gaussian sequence. // TODO deprecate this version
    double rn_gauss(double,double) const;
    /// Shuffle given points in array
    void ShufflePoints(std::vector<int>&) const;
    /// \return true if RN generator has been set up.
    bool IsSet() const;
    /// \return Value of RN generator seed
    int Seed() const;
  private:
    /// Allocate the RNG
    void allocateRng();

    Cpptraj::RNG* rng_;          ///< Hold the random number generator implemenatation.

    static RngType defaultType_; ///< Default RNG type
    //static int defaultSeed_;     ///< Default RNG seed
};
#endif
