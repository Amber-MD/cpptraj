#include "Action_Remap.h"
#include "CpptrajStdio.h"
#include "DataSet_1D.h"

void Action_Remap::Help() const {
  mprintf("\tdata <setname>\n"
          "  Re-map atoms according to the given reference data set which is of the format:\n"
          "    Reference[Target]\n"
          "  with atom numbering starting from 1. E.g. Reference[1] = 10 would mean remap\n"
          "  atom 10 in target to position 1.\n");
}

Action::RetType Action_Remap::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Get Keywords
  std::string mapsetname = actionArgs.GetStringKey("data");
  if (mapsetname.empty()) {
    mprinterr("Error: Atom map data set name not specified.\n");
    return Action::ERR;
  }
  // Get dataset
  DataSet* mapset = 0;
  if (!mapsetname.empty()) {
    mapset = init.DSL().GetDataSet( mapsetname );
    if (mapset == 0) {
      mprinterr("Error: Atom map set '%s' not found.\n", mapsetname.c_str());
      return Action::ERR;
    }
    if (mapset->Group() != DataSet::SCALAR_1D) {
      mprinterr("Error: Atom map set '%s' is not a 1D scalar set.\n", mapset->legend());
      return Action::ERR;
    }
    if (mapset->Size() < 1) {
      mprinterr("Error: Atom map set '%s' contains no data.\n", mapset->legend());
      return Action::ERR;
    }
    DataSet_1D const& ds = static_cast<DataSet_1D const&>( *mapset );
    Map_.reserve( ds.Size() );
    for (unsigned int i = 0; i != ds.Size(); i++)
      Map_.push_back( (int)ds.Dval(i) );
  }

  mprintf("    REMAP:");
  if (mapset != 0) {
    mprintf("Remapping atoms according to positions specified by data set '%s' (%zu atoms).\n",
            mapset->legend(), Map_.size());
  } else
    return Action::ERR; // Sanity check
  return Action::OK;
}

Action::RetType Action_Remap::Setup(ActionSetup& setup) {
  return Action::ERR;
}

Action::RetType Action_Remap::DoAction(int frameNum, ActionFrame& frm) {
  return Action::ERR;
}

