#ifndef INC_EXEC_PREPAREFORLEAP_H
#define INC_EXEC_PREPAREFORLEAP_H
#include "Exec.h"
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
    void LeapBond(int,int,Topology const&, CpptrajFile*) const;
    int IdentifySugar(int, Topology*, Frame const&, CharMask const&, CpptrajFile*) const;
    int FindTerByBonds(Topology*, CharMask const&) const;

    std::string leapunitname_;
};
#endif
