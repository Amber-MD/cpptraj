#ifndef INC_ENERGY_PME_RECIPPARAMS_H
#define INC_ENERGY_PME_RECIPPARAMS_H
class Box;
class EwaldOptions;
namespace Cpptraj {
namespace Energy {
/// Hold parameters for the reciprocal part of PME
/** NOTE: This is kept separate from PME_Recip so that classes that
  *       do the recip. calculation a little differently (like GIST_PME)
  *       can have access to the routines that assign NFFT etc.
  */
class PME_RecipParams {
  public:
    /// CONSTRUCTOR
    PME_RecipParams();
    /// Initilize recip options
    int InitRecip(EwaldOptions const&, int);
    /// Print recip options to stdout
    void PrintRecipOpts() const;
    /// Determine FFT grid points in X, Y, Z from box lengths
    int DetermineNfft(int&, int&, int&, Box const&) const;

    /// \return B-Spline order
    int Order() const { return order_; }
  private:
    static bool check_prime_factors(int);
    static int ComputeNFFT(double);

    int nfft_[3]; ///< Number of grid points for FFT in X, Y, Z
    int order_; ///< B-Spline order
    int debug_;
};
}
}
#endif
