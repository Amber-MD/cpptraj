#ifndef INC_EXEC_TOP_H
#define INC_EXEC_TOP_H
/*! \file Exec_Top.h
    \brief Exec classes for simple topology commands.
 */
#include "Exec.h"
/// Add topology to State.
class Exec_LoadParm : public Exec {
  public:
    Exec_LoadParm() : Exec(PARM) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_LoadParm(); }
    RetType Execute(CpptrajState& State, ArgList& argIn) {
      return (RetType)State.AddTopology( argIn.GetStringNext(), argIn );
    }
};
/// Print info for specified parm.
class Exec_ParmInfo : public Exec { 
  public:
    Exec_ParmInfo() : Exec(PARM) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_ParmInfo(); } 
    RetType Execute(CpptrajState&, ArgList&);
};
/// Print bond info for atoms in mask.
class Exec_BondInfo : public Exec {
  public:
    Exec_BondInfo() : Exec(PARM) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_BondInfo(); }
    RetType Execute(CpptrajState&, ArgList&);
};
/// Print angle info for atoms in mask.
class Exec_AngleInfo : public Exec {
  public:
    Exec_AngleInfo() : Exec(PARM) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_AngleInfo(); }
    RetType Execute(CpptrajState&, ArgList&);
};
/// Print dihedral info for atoms in mask.
class Exec_DihedralInfo : public Exec {
  public:
    Exec_DihedralInfo() : Exec(PARM) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_DihedralInfo(); }
    RetType Execute(CpptrajState&, ArgList&);
};
/// Print CHARMM improper info for atoms in mask.
class Exec_ImproperInfo : public Exec {
  public:
    Exec_ImproperInfo() : Exec(PARM) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_ImproperInfo(); }
    RetType Execute(CpptrajState&, ArgList&);
};
/// Print info for atoms in mask.
class Exec_AtomInfo : public Exec {
  public:
    Exec_AtomInfo() : Exec(PARM) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_AtomInfo(); }
    RetType Execute(CpptrajState&, ArgList&);
};
/// Print info for residues containing atoms in mask.
class Exec_ResInfo : public Exec {
  public:
    Exec_ResInfo() : Exec(PARM) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_ResInfo(); }
    RetType Execute(CpptrajState&, ArgList&);
};
/// Print info for molecules containing atoms in mask.
class Exec_MolInfo : public Exec {
  public:
    Exec_MolInfo() : Exec(PARM) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_MolInfo(); }
    RetType Execute(CpptrajState&, ArgList&);
};
/// Print the total charge of atoms in mask.
class Exec_ChargeInfo : public Exec {
  public:
    Exec_ChargeInfo() : Exec(PARM) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_ChargeInfo(); }
    RetType Execute(CpptrajState&, ArgList&);
};
/// Print the total mass of atoms in mask.
class Exec_MassInfo : public Exec {
  public:
    Exec_MassInfo() : Exec(PARM) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_MassInfo(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
