#include "Exec_Sequence.h"
#include "CpptrajStdio.h"

// Exec_Sequence::Help()
void Exec_Sequence::Help() const
{

}

/** Generate and build the specified sequence. */
int Exec_Sequence::generate_sequence(DataSet_Coords* OUT,
                                     DataSetList const& DSL,
                                     Sarray const& main_sequence,
                                     Sarray const& LibSetNames)
const
{
  // First, get all units in order.
  typedef std::vector<DataSet_Coords*> Uarray;
  Uarray Units;
  Units.reserve( main_sequence.size() );

  for (Sarray::const_iterator it = main_sequence.begin(); it != main_sequence.end(); ++it)
  {
    DataSet_Coords* unit = 0;
    if (LibSetNames.empty()) {
      // Look for name
      unit = (DataSet_Coords*)DSL.FindSetOfGroup( *it, DataSet::COORDINATES );
    } else {
      // Look for name[aspect]
      DataSet* ds = 0;
      for (Sarray::const_iterator name = LibSetNames.begin(); name != LibSetNames.end(); ++name)
      {
        MetaData meta(*name, *it);
        ds = DSL.CheckForSet( meta );
        if (ds != 0) break;
      }
      if (ds != 0) {
        if (ds->Group() != DataSet::COORDINATES) {
          mprinterr("Error: Set '%s' is not of type Coordinates.\n", ds->legend());
          return 1;
        }
        unit = (DataSet_Coords*)ds;
      }
    }
    if (unit == 0) {
      mprinterr("Error: Unit '%s' not found.\n", it->c_str());
      return 1;
    }
    // Needs to have connect associated data
    AssociatedData* ad = unit->GetAssociatedData(AssociatedData::CONNECT);
    if (ad == 0) {
      mprinterr("Error: Unit '%s' does not have CONNECT data.\n");
      return 1;
    }
    Units.push_back( unit );
  } // END loop over sequence
  mprintf("\tFound %zu units.\n", Units.size());
  return 0;
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
