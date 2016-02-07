#ifndef INC_ESANDER_H
#define INC_ESANDER_H
#ifdef USE_SANDERLIB
#include "Topology.h"
#include "sander.h"
class Energy_Sander {
  public:
    Energy_Sander() : top_pindex_(-1) {}
    ~Energy_Sander();
    int Initialize(Topology const&, Frame&);
    int CalcEnergy(Frame&);
    /// Enumerate different types of energy calcd via libsander
    enum Etype { TOTAL = 0, VDW, ELEC, GB, BOND, ANGLE, DIHEDRAL, VDW14, ELEC14,
                 CONSTRAINT, POLAR, HBOND, SURF, SCF, DISP, DVDL, ANGLE_UB,
                 IMP, CMAP, EMAP, LES, NOE, PB, RISM, CT, AMD_BOOST, N_ENERGYTYPES };
    /// \return Value of specified energy term.
    double Energy(Etype) const;
    /// \return Pointer to specified energy term.
    const double* Eptr(Etype) const;
    /// \return DataSet aspect string for given energy term.
    static std::string Easpect(Etype);
  private:
    static const char* Estring_[];
    sander_input input_;         ///< Sander input options
    pot_ene energy_;             ///< Sander energy terms
    FileName top_filename_;      ///< Current Topology file name
    std::vector<double> forces_; ///< Force array TODO use Frame
    int top_pindex_;             ///< Current Topology internal index.
};
#endif /* USE_SANDERLIB */
#endif
