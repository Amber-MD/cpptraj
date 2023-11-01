#include "Exec_Sequence.h"
#include "CpptrajStdio.h"
#include "AssociatedData_Connect.h"
#include "Structure/Builder.h"

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
  typedef std::vector<int> Iarray;
  Iarray connectAt0, connectAt1;
  connectAt0.reserve( Units.size() );
  connectAt1.reserve( Units.size() );
  int total_natom = 0;

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
    if (unit->Size() < 1) {
      mprinterr("Error: Unit '%s' is empty.\n", unit->legend());
      return 1;
    }
    if (unit->Size() > 1) {
      mprintf("Warning: Unit '%s' has more than 1 frame. Only using first frame.\n", unit->legend());
    }
    // Needs to have connect associated data
    AssociatedData* ad = unit->GetAssociatedData(AssociatedData::CONNECT);
    if (ad == 0) {
      mprinterr("Error: Unit '%s' does not have CONNECT data.\n", unit->legend());
      return 1;
    }
    AssociatedData_Connect const& C = static_cast<AssociatedData_Connect const&>( *ad );
    if (C.NconnectAtoms() < 2) {
      mprinterr("Error: Not enough connect atoms in unit '%s'\n", unit->legend());
      return 1;
    }
    // Update connect atom 1 indices based on their position in the sequence.
    // Do not update connect atom 0 indices since they will not yet be connected.
    connectAt0.push_back( C.Connect()[0] );
    connectAt1.push_back( C.Connect()[1] + total_natom );
    Units.push_back( unit );
    total_natom += unit->Top().Natom();
  } // END loop over sequence
  mprintf("\tFound %zu units.\n", Units.size());
  if (Units.empty()) {
    mprinterr("Error: No units.\n");
    return 1;
  }
  for (unsigned int idx = 0; idx < Units.size(); idx++)
    mprintf("\tUnit %s HEAD %i TAIL %i\n", Units[idx]->legend(), connectAt0[idx]+1, connectAt1[idx]+1);

  Topology combinedTop;
  combinedTop.SetDebug( debug_ );
  combinedTop.SetParmName( OUT->Meta().Name(), FileName() );
  combinedTop.AppendTop( Units.front()->Top() );
  //combinedTop.SetParmBox( Units // TODO
  combinedTop.Brief("Sequence topology:");

  Frame CombinedFrame = Units.front()->AllocateFrame();
  Units.front()->GetFrame(0, CombinedFrame);


  using namespace Cpptraj::Structure;
  Builder builder;
  for (unsigned int idx = 1; idx < Units.size(); idx++) {
    mprintf("\tConnect %s atom %i to %s atom %i\n",
            Units[idx-1]->legend(), connectAt1[idx-1]+1,
            Units[idx]->legend(),   connectAt0[idx]  +1);
    Frame mol1frm = Units[idx]->AllocateFrame();
    Units[idx]->GetFrame(0, mol1frm);
    int bondat0 = connectAt1[idx-1];
    int bondat1 = connectAt0[idx];
    if (bondat0 < 0 || bondat1 < 0) {
      mprinterr("Error: Invalid connect atom(s) between %s atom %i to %s atom %i\n",
                Units[idx-1]->legend(), bondat0+1, Units[idx]->legend(), bondat1+1);
      return 1;
    }
    if (builder.Combine( combinedTop, CombinedFrame, Units[idx]->Top(), mol1frm,
                         connectAt1[idx-1], connectAt0[idx] )) {
      mprinterr("Error: Sequence combine between units %u %s and %u %s failed.\n",
                idx, Units[idx-1]->legend(), idx+1, Units[idx]->legend());
      return 1;
    }
  }

  OUT->CoordsSetup(combinedTop, CombinedFrame.CoordsInfo());
  OUT->AddFrame( CombinedFrame );

  return 0;
}

// Exec_Sequence::Help()
void Exec_Sequence::Help() const
{
  mprintf("\tname <output set name> <unit0> <unit1> ...\n"
          "\t[{libset <libsetname>} ...]\n"
          "  Create a molecule from a sequence of units.\n");
}


// Exec_Sequence::Execute()
Exec::RetType Exec_Sequence::Execute(CpptrajState& State, ArgList& argIn)
{
  debug_ = State.Debug();
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

  // Execute
  if (generate_sequence(OUT, State.DSL(), main_sequence, LibSetNames)) {
    mprinterr("Error: Could not generate sequence.\n");
    return CpptrajState::ERR;
  }

  return CpptrajState::OK;
}
