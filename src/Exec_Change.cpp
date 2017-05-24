#include "Exec_Change.h"
#include "CpptrajStdio.h"

// Exec_Change::Help()
void Exec_Change::Help() const
{
  mprintf("\t[%s]\n"
          "\t{resname <value> <mask>}\n"
          "  Change specified parts of topology.\n", DataSetList::TopIdxArgs);
}

// Exec_Change::Execute()
Exec::RetType Exec_Change::Execute(CpptrajState& State, ArgList& argIn)
{
  // Change type
  enum ChangeType { UNKNOWN = 0, RESNAME };
  ChangeType type = UNKNOWN;
  std::string resname = argIn.GetStringKey("resname");
  if (!resname.empty()) type = RESNAME;
  if (type == UNKNOWN) {
    mprinterr("Error: No change type specified.\n");
    return CpptrajState::ERR;
  }
  Topology* parm = State.DSL().GetTopByIndex( argIn );
  if (parm == 0) return CpptrajState::ERR;
  int err = 0;
  switch (type) {
    case RESNAME : err = ChangeResidueName(resname, *parm, argIn); break;
    case UNKNOWN : err = 1; // sanity check
  }
  if (err != 0) return CpptrajState::ERR;
  return CpptrajState::OK;
}

int Exec_Change::ChangeResidueName(std::string const& name, Topology& topIn, ArgList& argIn)
const
{
  NameType rname( name );
  CharMask mask(argIn.GetMaskNext());
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
