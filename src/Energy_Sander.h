#ifndef INC_ESANDER_H
#define INC_ESANDER_H
#ifdef USE_SANDERLIB
#include "Topology.h"
#include "ArgList.h"
#include "sander.h"
/// Cpptraj interface to Sander API
class Energy_Sander {
  public:
    Energy_Sander();
    ~Energy_Sander();
    /// Set debug level.
    void SetDebug(int d) { debug_ = d; }
    /// Set input options from ArgList
    int SetInput(ArgList&);
    /// Initialize for given Topology and Frame
#   ifdef MPI
    int Initialize(Topology const&, Frame&, Parallel::Comm const&);
#   else
    int Initialize(Topology const&, Frame&);
#   endif
    /// Calculate energy for given Frame
    int CalcEnergy(Frame&);
    /// Calculate energy for given Frame, save forces
    int CalcEnergyForces(Frame&);
    /// Enumerate different types of energy calcd via libsander
    enum Etype { TOTAL = 0, VDW, ELEC, GB, BOND, ANGLE, DIHEDRAL, VDW14, ELEC14,
                 CONSTRAINT, POLAR, HBOND, SURF, CAVITY, SCF, DISP, DVDL, ANGLE_UB,
                 IMP, CMAP, EMAP, LES, NOE, PB, RISM, CT, AMD_BOOST, N_ENERGYTYPES };
    /// \return Value of specified energy term.
    double Energy(Etype) const;
    /// \return Pointer to specified energy term.
    const double* Eptr(Etype) const;
    /// \return DataSet aspect string for given energy term.
    static std::string Easpect(Etype);
    /// \return Label for given energy term.
    static const char* Elabel(Etype e) { return Estring_[e]; }
    /// \return True if corresponding energy term is active
    bool IsActive(Etype e) const { return isActive_[e]; }
    /// \return Temporary top file name
    FileName const& TopFilename() const { return top_filename_; }
  private:
    void SetDefaultInput();
    int WriteTop(Topology const&);
    int CommonInit(Topology const& topIn, Frame& fIn);

    static const char* Estring_[];
    sander_input input_;         ///< Sander input options
    pot_ene energy_;             ///< Sander energy terms
    FileName top_filename_;      ///< Current Topology file name
    std::vector<double> forces_; ///< Force array
    std::vector<bool> isActive_; ///< True if corresponding energy term is active.
    int debug_;                  ///< Debug level
    int top_pindex_;             ///< Current Topology internal index.
    bool specified_cut_;         ///< 'cut' was specified.
    bool specified_igb_;         ///< 'igb' was specified.
    bool specified_ntb_;         ///< 'ntb' was specified.
};
#endif /* USE_SANDERLIB */
#endif
