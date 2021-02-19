#ifndef INC_RNG_H
#define INC_RNG_H
namespace Cpptraj {
/// Interface to random number generators
class RNG {
  public:
    /// CONSTRUCTOR
    RNG();
    /// Initialize the random number generator with the given seed
    virtual void rn_set(int) = 0;
    /// Generate a random number between 0.0 and 1.0
    virtual double rn_gen() = 0;

    /// Initialize RN generator with 71277 (Amber default)
    void rn_set() { rn_set(71277); }

    /// Generate a pseudo-random Gaussian sequence.
    double rn_gauss(double,double);

    /// \return true if RN generator has been set up.
    bool IsSet() const { return (iseed_ != -1); }
    /// \return Value of RN generator seed
    int Seed() const { return iseed_; }
  private:
    int iseed_; ///< The random number generator seed. -1 means use time
};
}
#endif
