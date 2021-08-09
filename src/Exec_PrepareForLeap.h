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
    /// Set a reduced PDB res to glycam map when dat file not found.
    void SetGlycamPdbResMap();
    /// Load PDB res to glycam map from dat file
    int LoadGlycamPdbResMap(std::string const&);

    void LeapBond(int,int,Topology const&, CpptrajFile*) const;
    int CalcStereocenterTorsion(double&, int, Topology const&, Frame const&) const;
    int CalcAnomericRefTorsion(double&, int, int, Topology const&, Frame const&, std::vector<bool> const&) const;
    int CalcAnomericTorsion(double&, int, int, int, Topology const&, Frame const&, std::vector<bool> const&) const;
    int FindRemainingChainCarbons(std::vector<int>&, int, Topology const&, int,
                                  std::vector<bool> const&) const;
    int IdentifySugar(int, Topology&, Frame const&, CharMask const&, CpptrajFile*, std::set<BondType>&) const;
    int PrepareSugars(AtomMask&, Topology&, Frame const&, CpptrajFile*) const;
    int FindTerByBonds(Topology&, CharMask const&) const;
    int SearchForDisulfides(double, std::string const&, std::string const&, bool,
                            Topology&, Frame const&, CpptrajFile*) const;
    int ModifyCoords(Topology&, Frame&, bool, bool, char, std::string const&, std::string const&) const;

    std::string leapunitname_;
    typedef std::pair<NameType, char> PairType;
    typedef std::map<NameType, char> MapType;
    MapType pdb_to_glycam_; ///< Map PDB residue names to Glycam 1 char names
    bool errorsAreFatal_;   ///< If false, try to skip errors.
};
#endif
