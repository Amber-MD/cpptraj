#include "Exec_HmassRepartition.h"
#include "CharMask.h"
#include "CpptrajStdio.h"

// Exec_HmassRepartition::Help()
void Exec_HmassRepartition::Help() const
{
  mprintf("\t[%s] [<mask>] [hmass <hydrogen mass>]\n", DataSetList::TopIdxArgs);
  mprintf("  Perform hydrogen mass repartitioning for atoms selected by <mask>\n"
          "  (all solute atoms by default).\n");
}

// Exec_HmassRepartition::Execute()
Exec::RetType Exec_HmassRepartition::Execute(CpptrajState& State, ArgList& argIn)
{
  double hmass = argIn.getKeyDouble("hmass", 3.024);

  std::string maskArg = argIn.GetMaskNext();
  CharMask mask;
  if (mask.SetMaskString(maskArg)) {
    mprinterr("Error: Could not process mask string.\n");
    return CpptrajState::ERR;
  }
  // Get Topology
  Topology* topIn = State.DSL().GetTopByIndex( argIn );
  if (topIn == 0) return CpptrajState::ERR;

  mprintf("\tPerforming hydrogen mask repartitioning for atoms selected by mask '%s'\n",
          mask.MaskString());
  mprintf("\tTopology is '%s'\n", topIn->c_str());
  mprintf("\tHydrogen masses will be changed to %f amu\n", hmass);

  return CpptrajState::OK;
}
