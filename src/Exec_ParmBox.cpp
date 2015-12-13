#include "Exec_ParmBox.h"
#include "CpptrajStdio.h"

void Exec_ParmBox::Help() const {
  mprintf("\t[%s] [nobox] [truncoct]\n", DataSetList::TopIdxArgs);
  mprintf("\t[x <xval>] [y <yval>] [z <zval>] [alpha <a>] [beta <b>] [gamma <g>]\n"
          "  Set the box info for specified topology (currently only relevant for Amber\n"
          "  Topology/ASCII coords). If 'nobox' is specified, remove box info. If\n"
          "  'truncoct' specified, set truncated octahedron with lengths = <xval>.\n");
}

Exec::RetType Exec_ParmBox::Execute(CpptrajState& State, ArgList& argIn) {
  Box pbox;
  bool nobox = false;
  if ( argIn.hasKey("nobox") )
    nobox = true;
  else {
    pbox.SetX( argIn.getKeyDouble("x",0) );
    pbox.SetY( argIn.getKeyDouble("y",0) );
    pbox.SetZ( argIn.getKeyDouble("z",0) );
    pbox.SetAlpha( argIn.getKeyDouble("alpha",0) );
    pbox.SetBeta(  argIn.getKeyDouble("beta",0)  );
    pbox.SetGamma( argIn.getKeyDouble("gamma",0) );
  }
  Topology* parm = State.DSL().GetTopByIndex( argIn );
  if (parm == 0) return CpptrajState::ERR;
  if (nobox)
    mprintf("\tRemoving box information from parm %i:%s\n", parm->Pindex(), parm->c_str());
  else
    // Fill in missing parm box information from specified parm
    pbox.SetMissingInfo( parm->ParmBox() );
  if (argIn.hasKey("truncoct")) pbox.SetTruncOct();
  parm->SetParmBox( pbox );
  parm->ParmBox().PrintInfo();
  return CpptrajState::OK;
}
