#include "Exec_Graft.h"
#include "CpptrajStdio.h"
#include "DataSet_Coords.h"

// Exec_Graft::Help()
void Exec_Graft::Help() const
{
  mprintf("\tsrc <source COORDS> srcmask <srcmask>\n"
          "\ttgt <target COORDS>\n");
}

// Exec_Graft::Execute()
Exec::RetType Exec_Graft::Execute(CpptrajState& State, ArgList& argIn)
{
  // Get source coords
  std::string kw = argIn.GetStringKey("src");
  if (kw.empty()) {
    mprinterr("Error: Source COORDS must be specified with 'src'.\n");
    return CpptrajState::ERR;
  }
  DataSet_Coords* srcCoords = (DataSet_Coords*)State.DSL().FindSetOfGroup(kw, DataSet::COORDINATES);
  if (srcCoords == 0) {
    mprinterr("Error: Source COORDS %s not found.\n", kw.c_str());
    return CpptrajState::ERR;
  }
  // Get target coords
  kw = argIn.GetStringKey("tgt");
  if (kw.empty()) {
    mprinterr("Error: Target COORDS must be specified with 'tgt'.\n");
    return CpptrajState::ERR;
  }
  DataSet_Coords* tgtCoords = (DataSet_Coords*)State.DSL().FindSetOfGroup(kw, DataSet::COORDINATES);
  if (tgtCoords == 0) {
    mprinterr("Error: Target COORDS %s not found.\n", kw.c_str());
    return CpptrajState::ERR;
  }
  // Get other keywords
  AtomMask srcMask;
  srcMask.SetMaskString( argIn.GetStringKey("srcmask") );

  // Info
  mprintf("\tSource coords : %s\n", srcCoords->legend());
  mprintf("\tTarget coords : %s\n", tgtCoords->legend());
  mprintf("\tSource mask   : %s\n", srcMask.MaskString());

  return CpptrajState::OK;
}
