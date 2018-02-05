#ifndef INC_EWALD_PARTICLEMESH_H
#define INC_EWALD_PARTICLEMESH_H
#ifdef LIBPME
#include "Ewald.h"
#include "libpme_standalone.h"
/// Class for calculating electrostatics with particle mesh Ewald.
class Ewald_ParticleMesh : public Ewald {
  public:
    Ewald_ParticleMesh();
    /// Box, cut, dsum tol, ew coeff, NB skin, erfc dx, order, debug, nfft
    int Init(Box const&, double, double, double, double, double, int, int, const int*);
    // ----- Inherited ---------------------------
    int Setup(Topology const&, AtomMask const&);
    double CalcEnergy(Frame const&, AtomMask const&, double&); // TODO const?
  private:
    /// Based on given length return number of grid points that is power of 2, 3, or 5
    static int ComputeNFFT(double);
    /// Determine grid points for FFT in each dimension
    int DetermineNfft(int&, int&, int&, Box const&) const;
    /// Particle mesh Ewald reciprocal energy
    double Recip_ParticleMesh(libpme::Mat<double> const&, libpme::Mat<double> const&, Box const&);

    libpme::Mat<double> coordsD_;   ///< Hold coordinates for selected atoms
    libpme::Mat<double> chargesD_;  ///< Will point to Ewald::Charge_

    int nfft_[3]; ///< Number of FFT grid points in each direction
    int order_;   ///< PME B spline order
};
#endif /* LIBPME */
#endif
