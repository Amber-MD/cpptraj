#ifndef INC_RNG_H
#define INC_RNG_H
namespace Cpptraj {
/// Interface to random number generators
class RNG {
  public:
    /// CONSTRUCTOR
    RNG();
    /// Generate a random number between 0.0 and 1.0
    virtual double Generate() = 0;

    /// Initialize the random number generator with the given seed
    void Set_Seed(int);
    /// Initialize RN generator with 71277 (Amber default)
    void Set_Seed() { Set_Seed(71277); }

    /// Generate a pseudo-random Gaussian sequence.
    double GenerateGauss(double,double);

    /// \return true if RN generator has been set up.
    bool IsSet() const { return (iseed_ != -1); }
    /// \return Value of RN generator seed
    int Seed() const { return iseed_; }
  protected:
    /// Any internal setup that needs to be done by the RNG
    virtual int SetupRng() = 0;
  private:
    int iseed_; ///< The random number generator seed. -1 means use time
};
}
#endif
