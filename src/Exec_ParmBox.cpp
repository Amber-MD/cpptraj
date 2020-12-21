#include "Exec_ParmBox.h"
#include "CpptrajStdio.h"
#include "BoxArgs.h"

void Exec_ParmBox::Help() const {
  mprintf("\t[%s]\n", DataSetList::TopIdxArgs);
  mprintf("\t{ nobox |\n"
          "\t  %s |\n"
          "\t  %s}\n", BoxArgs::Keywords_XyzAbg(), BoxArgs::Keywords_TruncOct());
  mprintf("  Set the box info for specified topology (currently only relevant for Amber\n"
          "  Topology/ASCII coords). If 'nobox' is specified, remove box info. If\n"
          "  'truncoct' specified, set truncated octahedron with lengths = <xval>.\n");
}

Exec::RetType Exec_ParmBox::Execute(CpptrajState& State, ArgList& argIn) {
  BoxArgs boxArgs;
  bool nobox = false;
  if ( argIn.hasKey("nobox") )
    nobox = true;
  else {
    if (boxArgs.SetBoxArgs( argIn )) return CpptrajState::ERR;
  }
  Topology* parm = State.DSL().GetTopByIndex( argIn );
  if (parm == 0) return CpptrajState::ERR;
  Box pbox;
  if (nobox)
    mprintf("\tRemoving box information from parm %i:%s\n", parm->Pindex(), parm->c_str());
  else {
    // Fill in missing parm box information from specified parm
    if (boxArgs.SetMissingInfo( parm->ParmBox() )) return CpptrajState::ERR;
    pbox.SetupFromXyzAbg( boxArgs.XyzAbg() );
  }
  parm->SetParmBox( pbox );
  parm->ParmBox().PrintInfo();
  return CpptrajState::OK;
}
