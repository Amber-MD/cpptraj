#include "Action_Remap.h"
#include "CpptrajStdio.h"
#include "DataSet_1D.h"

// CONSTRUCTOR
Action_Remap::Action_Remap() : newParm_(0) {}

// DESTRUCTOR
Action_Remap::~Action_Remap() {
  if (newParm_ != 0) delete newParm_;
}

// Action_Remap::Help()
void Action_Remap::Help() const {
  mprintf("\tdata <setname>\n");
  mprintf("%s", ActionTopWriter::Keywords());
  mprintf("  Re-map atoms according to the given reference data set which is of the format:\n"
          "    Reference[Target]\n"
          "  with atom numbering starting from 1. E.g. Reference[1] = 10 would mean remap\n"
          "  atom 10 in target to position 1.\n");
  mprintf("%s", ActionTopWriter::Options());
}

// Action_Remap::Init()
Action::RetType Action_Remap::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Get Keywords
  std::string mapsetname = actionArgs.GetStringKey("data");
  if (mapsetname.empty()) {
    mprinterr("Error: Atom map data set name not specified.\n");
    return Action::ERR;
  }
  topWriter_.InitTopWriter(actionArgs, "re-mapped", debugIn);
  // Get dataset
  DataSet* mapset = init.DSL().GetDataSet( mapsetname );
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
  // User atom #s start from 1
  for (unsigned int i = 0; i != ds.Size(); i++)
    Map_.push_back( (int)ds.Dval(i) - 1 );

  mprintf("    REMAP:\n"
          "\tRemapping atoms according to positions specified by data set '%s' (%zu atoms).\n",
          mapset->legend(), Map_.size());
  topWriter_.PrintOptions();
  return Action::OK;
}

// Action_Remap::Setup()
Action::RetType Action_Remap::Setup(ActionSetup& setup) {
  if (setup.Top().Natom() != (int)Map_.size()) {
    mprintf("Warning: Topology '%s' size (%i) does not match map size (%zu). Skipping.\n",
            setup.Top().c_str(), setup.Top().Natom(), Map_.size());
    return Action::SKIP;
  }
  // Attempt to create remapped topology
  if (newParm_ != 0) delete newParm_;
  newParm_ = setup.Top().ModifyByMap( Map_ );
  if (newParm_ == 0) {
    mprinterr("Error: Could not create re-mapped topology.\n");
    return Action::ERR;
  }
  setup.SetTopology( newParm_ );
  newParm_->Brief("Re-mapped topology:");
  // Allocate space for new frame
  newFrame_.SetupFrameV(setup.Top().Atoms(), setup.CoordInfo());
  // Write output topology if specified
  topWriter_.WriteTops( setup.Top() );

  return Action::MODIFY_TOPOLOGY;
}

// Action_Remap::DoAction()
Action::RetType Action_Remap::DoAction(int frameNum, ActionFrame& frm) {
  newFrame_.SetCoordinatesByMap( frm.Frm(), Map_ );
  frm.SetFrame( &newFrame_ );
  return Action::MODIFY_COORDS;
}
