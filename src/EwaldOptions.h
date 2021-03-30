#ifndef INC_EWALDOPTIONS_H
#define INC_EWALDOPTIONS_H
class ArgList;
/// Class that will hold options for Ewald-derived classes.
class EwaldOptions {
  public:
    EwaldOptions();
    /// Type of options
    enum OptType { NOT_SET = 0, REG_EWALD, PME };

    static const char* KeywordsCommon1();
    static const char* KeywordsCommon2();
    static const char* KeywordsRegEwald();
    static const char* KeywordsPME();
    static const char* KeywordsLjpme();

    /// Enable/disable LJ pme option
    void AllowLjPme(bool);
    /// Get options from arglist
    int GetOptions(OptType, ArgList&, const char*);
    /// Print options to stdout
    void PrintOptions() const;

    OptType Type() const { return type_; }

    // Options common to all Ewald
    double Cutoff()     const { return cutoff_; }
    double DsumTol()    const { return dsumtol_; }
    double EwCoeff()    const { return ewcoeff_; }
    double ErfcDx()     const { return erfcDx_; }
    double SkinNB()     const { return skinnb_; }

    // Lennard-Jones (VDW) options 
    double LwCoeff()    const { return lwcoeff_; }
    double LJ_SwWidth() const { return ljswidth_; }

    // Options specific to regular Ewald
    double RsumTol()    const { return rsumtol_; }
    double MaxExp()     const { return maxexp_; }
    int Mlimits1()      const { return mlimits1_; }
    int Mlimits2()      const { return mlimits2_; }
    int Mlimits3()      const { return mlimits3_; }

    // Options specific to PME
    int SplineOrder()   const { return npoints_; }
    int Nfft1()         const { return nfft1_; }
    int Nfft2()         const { return nfft2_; }
    int Nfft3()         const { return nfft3_; }
  private:
    static int GetCommaSeparatedArgs(ArgList&, const char*, int&, int&, int&, int);

    bool allowLjPme_; ///< True if LJ pme is allowed.
    OptType type_;    ///< Options type (if any)
    // Common options
    double cutoff_;   ///< Direct space cutoff
    double dsumtol_;  ///< Direct sum tolerance
    double ewcoeff_;  ///< Ewald coefficient for electrostatics
    double erfcDx_;   ///< Spline table spacing for ERFC function
    double skinnb_;   ///< Nonbond "skin" (used for consitency with Amber)
    // LJ options
    double lwcoeff_;  ///< Ewald coefficient for LJ
    double ljswidth_; ///< LJ switch width
    // Ewald options
    double rsumtol_;  ///< Reciprocal sum tolerance
    double maxexp_;
    int mlimits1_;    ///< Number of reciprocal vectors in the X direction
    int mlimits2_;    ///< Number of reciprocal vectors in the Y direction
    int mlimits3_;    ///< Number of reciprocal vectors in the Z direction
    // PME options
    int npoints_;     ///< Spline order for PME grids
    int nfft1_;       ///< Number of grid points in the X direction
    int nfft2_;       ///< Number of grid points in the Y direction
    int nfft3_;       ///< Number of grid points in the Z direction
};
#endif
