#include "Exec_ParmWrite.h"
#include "CpptrajStdio.h"
#include "ParmFile.h"

void Exec_ParmWrite::Help() const {
  mprintf("\tout <filename> [{%s | crdset <setname>}] [<fmt>] [nochamber]\n",
          DataSetList::TopIdxArgs);
  mprintf("  Write specified topology or topology from COORDS set to <filename>.\n");
  ParmFile::WriteOptions();
}

Exec::RetType Exec_ParmWrite::Execute(CpptrajState& State, ArgList& argIn) {
  std::string outfilename = argIn.GetStringKey("out");
  if (outfilename.empty()) {
    mprinterr("Error: No output filename specified (use 'out <filename>').\n");
    return CpptrajState::ERR;
  }
  int err = 0;
  ParmFile pfile;
  // Check if a COORDS data set was specified.
  std::string crdset = argIn.GetStringKey("crdset");
  if (crdset.empty()) {
    Topology* parm = State.DSL()->GetTopByIndex( argIn );
    if (parm == 0) return CpptrajState::ERR;
    err = pfile.WriteTopology( *parm, outfilename, argIn, ParmFile::UNKNOWN_PARM, State.Debug() );
  } else {
    DataSet_Coords* ds = (DataSet_Coords*)State.DSL()->FindCoordsSet(crdset);
    if (ds == 0) return CpptrajState::ERR;
    mprintf("\tUsing topology from data set '%s'\n", ds->legend());
    err = pfile.WriteTopology(ds->Top(), outfilename, argIn, ParmFile::UNKNOWN_PARM, State.Debug());
  }
  if (err != 0)
    return CpptrajState::ERR;
  return CpptrajState::OK;
}
