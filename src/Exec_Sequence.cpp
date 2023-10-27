#include "Exec_Sequence.h"
#include "CpptrajStdio.h"

// Exec_Sequence::Help()
void Exec_Sequence::Help() const
{

}

// Exec_Sequence::Execute()
Exec::RetType Exec_Sequence::Execute(CpptrajState& State, ArgList& argIn)
{
  // Args
  Sarray LibSetNames;
  std::string libsetname = argIn.GetStringKey("libset");
  while (!libsetname.empty()) {
    LibSetNames.push_back( libsetname );
    libsetname = argIn.GetStringKey("libset");
  }
  std::string dsname = argIn.GetStringKey("name");
  if (dsname.empty()) {
    mprinterr("Error: No output set name specified with 'name'\n");
    return CpptrajState::ERR;
  }
  DataSet_Coords* OUT = (DataSet_Coords*)State.DSL().AddSet(DataSet::COORDS, MetaData(dsname));
  if (OUT == 0) {
    mprinterr("Error: Could not create output COORDS set named '%s'\n", dsname.c_str());
    return CpptrajState::ERR;
  }

  // Get the actual sequence from remaining args.
  Sarray main_sequence;
  ArgList remaining = argIn.RemainingArgs();
  std::string unit = remaining.GetStringNext();
  while (!unit.empty()) {
    main_sequence.push_back( unit );
    unit = remaining.GetStringNext();
  }
  if (main_sequence.empty()) {
    mprinterr("Error: No units specified.\n");
    return CpptrajState::ERR;
  }

  // Info
  if (!LibSetNames.empty()) {
    mprintf("\tLibrary set names:");
    for (Sarray::const_iterator it = LibSetNames.begin(); it != LibSetNames.end(); ++it)
      mprintf(" %s", it->c_str());
    mprintf("\n");
  }
  mprintf("\tMain sequence:");
  for (Sarray::const_iterator it = main_sequence.begin(); it != main_sequence.end(); ++it)
    mprintf(" %s", it->c_str());
  mprintf("\n");
  mprintf("\tOutput set name : %s\n", OUT->legend());

  return CpptrajState::OK;
}
