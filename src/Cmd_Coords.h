#ifndef INC_CMD_COORDS_H
#define INC_CMD_COORDS_H
#include "Cmd.h"
void Help_CrdAction();
void Help_CrdOut();
void Help_LoadCrd();
void Help_LoadTraj();
void Help_CombineCoords();

/// Perform action on given COORDS dataset
/** NOTE: CrdAction is implemented differently because it requires access to
  *       Command::SearchTokenType()
  */
Cmd::RetType CrdAction(CpptrajState&, ArgList&, DataSet_Coords*, Action*, TrajFrameCounter const&);
/// Write out COORDS dataset
Cmd::RetType CrdOut(CpptrajState&, ArgList&, Cmd::AllocType);
/// Load single trajectory as DataSet_Coords
Cmd::RetType LoadCrd(CpptrajState&, ArgList&, Cmd::AllocType);
/// Load trajectory or convert input traj list to TRAJ data set
Cmd::RetType LoadTraj(CpptrajState&, ArgList&, Cmd::AllocType);
/// Combine two or more COORDS data sets
Cmd::RetType CombineCoords(CpptrajState&, ArgList&, Cmd::AllocType);
#endif
