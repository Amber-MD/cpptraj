#ifndef INC_EXEC_PREPAREFORLEAP_H
#define INC_EXEC_PREPAREFORLEAP_H
#include "Exec.h"
class BondType;
namespace Cpptraj {
namespace Structure {
class SugarBuilder;
class ResStatArray;
class ResNameIndices;
}
}
/// Do common tasks to prepare a structure to be loaded into tleap 
class Exec_PrepareForLeap : public Exec {
  public:
    Exec_PrepareForLeap();
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_PrepareForLeap(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    typedef std::vector<int> Iarray;
    typedef std::set<NameType> SetType;
    typedef std::pair<NameType,Iarray> RpairType;
    typedef std::map<NameType,Iarray> RmapType;

    /// Set PDB residue names recognized by Amber FFs
    void SetPdbResNames();
    /// Load PDB residue names recognized by Amber FFs from dat file
    int LoadPdbResNames(std::string const&);
    /// \return true if residue name is a recognized PDB name
    bool IsRecognizedPdbRes(NameType const&, Cpptraj::Structure::SugarBuilder const&) const;
    /// \return Array of residue nums with unrecognized names
    Iarray GetUnrecognizedPdbResidues(Topology const&, Cpptraj::Structure::SugarBuilder const&) const;
    /// \return Array indices of isolated unrecognized residues
    Iarray GetIsolatedUnrecognizedResidues(Topology const&, Iarray const&) const;

    /// Try to determine where TER cards should be placed based on bonds
    int FindTerByBonds(Topology&, CharMask const&) const;

    /// Run leap to generate topology, perform any modifications
    int RunLeap(std::string const&, std::string const&) const;

    /// Download missing parameters
    int DownloadParameters(Cpptraj::Structure::ResStatArray&, RmapType const&,
                           Topology const&, CpptrajFile*, std::vector<BondType>&) const;

    // -----------------------
    SetType pdb_res_names_;      ///< PDB residue names recognized by Amber FFs

    std::string leapunitname_;   ///< Unit name to use when loading molecule
    bool errorsAreFatal_;        ///< If false, try to skip errors.
    bool downloadParams_;        ///< If true, try to download parameters for missing residues
    bool bondUnknownResidues_;   ///< If true, try to create bonds to unknown residues
    int debug_;                  ///< Debug level
    std::string solventResName_; ///< Solvent residue name
    std::string parameterURL_;   ///< URL to download parameters from
};
#endif
