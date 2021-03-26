#ifndef INC_PMEOPTIONS_H
#define INC_PMEOPTIONS_H
class ArgList;
/// Class that will hold common PME options for Ewald_ParticleMesh class.
class PmeOptions {
  public:
    PmeOptions();

    static const char* Keywords1();
    static const char* Keywords2();
    static const char* Keywords3();

    /// Enable/disable LJ pme option
    void AllowLjPme(bool);
    /// Get options from arglist
    int GetOptions(ArgList&, const char*);
  private:
    bool allowLjPme_; ///< True if LJ pme is allowed.

    double cutoff_; ///< Direct space cutoff
    double dsumtol_; ///< Direct sum tolerance
    double ewcoeff_; ///< Ewald coefficient for electrostatics
    double lwcoeff_; ///< Ewald coefficient for LJ
    double ljswidth_; ///< LJ switch width
    double skinnb_;   ///< Nonbond "skin"
    double erfcDx_;   ///< Spline table spacing for ERFC function
    int npoints_;     ///< Spline order for PME grids
    int nfft1_;       ///< Number of grid points in X direction
    int nfft2_;       ///< Number of grid points in Y direction
    int nfft3_;       ///< Number of grid points in Z direction
};
#endif
