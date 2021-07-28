#ifndef INC_EXEC_PREPAREFORLEAP_H
#define INC_EXEC_PREPAREFORLEAP_H
#include "Exec.h"
#include <set>
#include <map>
class CharMask;
class CpptrajFile;
/// Do common tasks to prepare a structure to be loaded into tleap 
class Exec_PrepareForLeap : public Exec {
  public:
    Exec_PrepareForLeap() : Exec(COORDS) { SetHidden(true); }
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_PrepareForLeap(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    /// Set a reduced PDB res to glycam map when dat file not found.
    void SetGlycamPdbResMap();
    /// Load PDB res to glycam map from dat file
    int LoadGlycamPdbResMap(std::string const&);

    void LeapBond(int,int,Topology const&, CpptrajFile*) const;
    int IdentifySugar(int, Topology*, Frame const&, CharMask const&, CpptrajFile*, std::set<BondType>&) const;
    int FindTerByBonds(Topology*, CharMask const&) const;

    std::string leapunitname_;
    typedef std::pair<NameType, char> PairType;
    typedef std::map<NameType, char> MapType;
    MapType pdb_to_glycam_; ///< Map PDB residue names to Glycam 1 char names
};
#endif
