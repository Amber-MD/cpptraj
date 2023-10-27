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

  return CpptrajState::OK;
}
