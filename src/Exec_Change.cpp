#include "Exec_Change.h"
#include "CpptrajStdio.h"

// Exec_Change::Help()
void Exec_Change::Help() const
{
  mprintf("\t[%s]\n"
          "\t{resname from <mask> to <value>}\n"
          "  Change specified parts of topology.\n", DataSetList::TopIdxArgs);
}

// Exec_Change::Execute()
Exec::RetType Exec_Change::Execute(CpptrajState& State, ArgList& argIn)
{
  // Change type
  enum ChangeType { UNKNOWN = 0, RESNAME, ATOMNAME };
  ChangeType type = UNKNOWN;
  if (argIn.hasKey("resname"))
    type = RESNAME;
  else if (argIn.hasKey("atomname"))
    type = ATOMNAME;
  if (type == UNKNOWN) {
    mprinterr("Error: No change type specified.\n");
    return CpptrajState::ERR;
  }
  Topology* parm = State.DSL().GetTopByIndex( argIn );
  if (parm == 0) return CpptrajState::ERR;
  int err = 0;
  switch (type) {
    case RESNAME : err = ChangeResidueName(*parm, argIn); break;
    case UNKNOWN : err = 1; // sanity check
  }
  if (err != 0) return CpptrajState::ERR;
  return CpptrajState::OK;
}

// Exec_Change::ChangeResidueName()
int Exec_Change::ChangeResidueName(Topology& topIn, ArgList& argIn)
const
{
  // Name to change to.
  std::string name = argIn.GetStringKey("to");
  if (name.empty()) {
    mprinterr("Error: Specify residue name to change to ('to <name>').\n");
    return 1;
  }
  NameType rname( name );
  // Residues to change
  std::string mexpr = argIn.GetStringKey("from");
  if (mexpr.empty()) {
    mprinterr("Error: Specify residue(s) to change names of ('from <mask>').\n");
    return 1;
  }
  CharMask mask(mexpr);
  if (topIn.SetupCharMask( mask )) return 1;
  mask.MaskInfo();
  if (mask.None()) {
    mprinterr("Error: No atoms selected by mask.\n");
    return 1;
  }
  for (int res = 0; res != topIn.Nres(); res++)
    if ( mask.AtomsInCharMask( topIn.Res(res).FirstAtom(), topIn.Res(res).LastAtom()-1 ) )
    {
      mprintf("\tChanging residue %s to %s\n", topIn.Res(res).c_str(), *rname);
      topIn.SetRes(res).SetName( rname );
    }
  return 0;
}
