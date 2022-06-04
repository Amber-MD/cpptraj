#ifndef INC_EXEC_PREPAREFORLEAP_H
#define INC_EXEC_PREPAREFORLEAP_H
#include "Exec.h"
#include <set>
#include <map>
#include <vector>
class AtomMask;
class CharMask;
class CpptrajFile;
/// Do common tasks to prepare a structure to be loaded into tleap 
class Exec_PrepareForLeap : public Exec {
  public:
    Exec_PrepareForLeap();
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_PrepareForLeap(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    typedef std::vector<int> Iarray;

    enum ResStatType { UNKNOWN = 0,
                       VALIDATED,
                       SUGAR_UNRECOGNIZED_LINK_RES,
                       SUGAR_UNRECOGNIZED_LINKAGE,
                       SUGAR_NO_LINKAGE,
                       SUGAR_NO_CHAIN_FOR_LINK,
                       SUGAR_NAME_MISMATCH,
                       //SUGAR_MISSING_C1X,
                       SUGAR_SETUP_FAILED };
    typedef std::vector<ResStatType> ResStatArray;

    
    
    /// Set PDB residue names recognized by Amber FFs
    void SetPdbResNames();
    /// Load PDB residue names recognized by Amber FFs from dat file
    int LoadPdbResNames(std::string const&);

    void LeapBond(int,int,Topology const&, CpptrajFile*) const;
//    int CalcStereocenterTorsion(double&, int, Topology const&, Frame const&) const;
    /// Try to find any missing bonds to C1 atoms
//    int FindSugarC1Linkages(Sugar const&, Topology&, Frame const&) const;

    int FindTerByBonds(Topology&, CharMask const&) const;
    int SearchForDisulfides(double, std::string const&, std::string const&, bool,
                            Topology&, Frame const&, CpptrajFile*);
    /// \return true if residue name is a recognized PDB name
    bool IsRecognizedPdbRes(NameType const&) const;
    /// \return Array of residue nums with unrecognized names
    Iarray GetUnrecognizedPdbResidues(Topology const&) const;
    /// \return Array indices of isolated unrecognized residues
    Iarray GetIsolatedUnrecognizedResidues(Topology const&, Iarray const&) const;
    /// Remove specified atoms
    int ModifyCoords(Topology&, Frame&, bool, std::string const&, std::string const&,
                     std::string const&, Iarray const&) const;
    int RemoveHydrogens(Topology&, Frame&) const;
    /** Try to determine protonation state of histidines from any hydrogens present. */
    int DetermineHisProt(Topology&,
                         NameType const&, NameType const&,
                         NameType const&, NameType const&, NameType const&, NameType const&) const;
    /// Run leap to generate topology, perform any modifications
    int RunLeap(std::string const&, std::string const&) const;
    /// Print a warning for residues that will need modification after leap
    static void LeapFxnGroupWarning(Topology const&, int);


    typedef std::set<NameType> SetType;
    SetType pdb_res_names_; ///< PDB residue names recognized by Amber FFs

    ResStatArray resStat_;  ///< Contain status of each residue

    static void PrintAtomNameMap(const char*, std::vector<NameMapType> const&);

    std::string leapunitname_;
    bool errorsAreFatal_;   ///< If false, try to skip errors.
    int debug_;             ///< Debug level
    std::string solventResName_; ///< Solvent residue name
};
#endif
