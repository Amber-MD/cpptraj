#include "Exec_ParmStrip.h"
#include "CpptrajStdio.h"

void Exec_ParmStrip::Help() const {
  mprintf("\t<mask> [%s]\n", DataSetList::TopIdxArgs);
  mprintf("  Strip atoms in mask from specified topology (first by default).\n");
}

Exec::RetType Exec_ParmStrip::Execute(CpptrajState& State, ArgList& argIn) {
    Topology* parm = State.DSL()->GetTopByIndex( argIn );
  if (parm == 0) return CpptrajState::ERR;
  // Check if this topology has already been used to set up an input
  // trajectory, as this will break the traj read.
  bool topology_in_use = false;
  const char* fname = 0;
  for (TrajinList::trajin_it tIn = State.InputTrajList().trajin_begin();
                             tIn != State.InputTrajList().trajin_end(); ++tIn)
    if ( (*tIn)->Traj().Parm() == parm ) {
      topology_in_use = true;
      fname = (*tIn)->Traj().Filename().full();
      break;
    }
  if (!topology_in_use) {
    for (TrajinList::ensemble_it eIn = State.InputTrajList().ensemble_begin();
                                 eIn != State.InputTrajList().ensemble_end(); ++eIn)
      if ( (*eIn)->Traj().Parm() == parm ) {
        topology_in_use = true;
        fname = (*eIn)->Traj().Filename().full();
        break;
      }
  }
  if (topology_in_use) {
    mprinterr("Error: Topology '%s' has already been used to set up trajectory '%s'.\n"
              "Error:   To strip this topology use the 'strip' action.\n",
              parm->c_str(), fname);
    return CpptrajState::ERR;
  }
  AtomMask tempMask( argIn.GetMaskNext() );
  // Since want to keep atoms outside mask, invert selection
  tempMask.InvertMaskExpression();
  if (parm->SetupIntegerMask( tempMask )) return CpptrajState::ERR;
  mprintf("\tStripping atoms in mask [%s] (%i) from %s\n",tempMask.MaskString(),
           parm->Natom() - tempMask.Nselected(), parm->c_str());
  Topology* tempParm = parm->modifyStateByMask(tempMask);
  if (tempParm==0) {
    mprinterr("Error: %s: Could not strip parm.\n", argIn.Command());
    return CpptrajState::ERR;
  } else {
    // Replace parm with stripped version
    *parm = *tempParm;
    parm->Brief("Stripped parm:");
    delete tempParm;
  }
  return CpptrajState::OK;
}
