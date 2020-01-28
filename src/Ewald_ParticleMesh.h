#ifndef INC_EWALD_PARTICLEMESH_H
#define INC_EWALD_PARTICLEMESH_H
#ifdef LIBPME
#include "Ewald.h"
#include "helpme_standalone.h"
/// Class for calculating electrostatics with particle mesh Ewald.
class Ewald_ParticleMesh : public Ewald {
  public:
    Ewald_ParticleMesh();
    /// Box, cut, dsum tol, ew coeff, lj ew coeff, switch width, NB skin, erfc dx, order, dbg, nfft
    int Init(Box const&, double, double, double, double, double, double, double,
             int, int, const int*);
    // ----- Inherited ---------------------------
    int Setup(Topology const&, AtomMask const&);
    int CalcNonbondEnergy(Frame const&, AtomMask const&, double&, double&);
  private:
    typedef Ewald::Darray Darray;
    /// Based on given length return number of grid points that is power of 2, 3, or 5
    static int ComputeNFFT(double);
    /// Determine grid points for FFT in each dimension
    int DetermineNfft(int&, int&, int&, Box const&) const;
    /// Particle mesh Ewald reciprocal energy
    double Recip_ParticleMesh(Box const&);
    /// Particle mesh Ewald LJ recip energy
    double LJ_Recip_ParticleMesh(Box const&);

    Darray coordsD_;   ///< Hold coordinates for selected atoms

    int nfft_[3]; ///< Number of FFT grid points in each direction
    int order_;   ///< PME B spline order

    PMEInstanceD pme_object_;
    PMEInstanceD pme_vdw_;
};
#endif /* LIBPME */
#endif
