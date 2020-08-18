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
  mprintf("\t{data <setname> | mask <mask> toindex <#>\n");
  mprintf("%s", ActionTopWriter::Keywords());
  mprintf("  Re-map atoms according to the given reference data set which is of the format:\n"
          "    Reference[Target]\n"
          "  with atom numbering starting from 1. E.g. Reference[1] = 10 would mean remap\n"
          "  atom 10 in target to position 1.\n"
          "  If 'mask' is specified, re-map the atoms in <mask> to positions starting\n"
          "  from <#>.\n");
  mprintf("%s", ActionTopWriter::Options());
}

// Action_Remap::Init()
Action::RetType Action_Remap::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Get Keywords
  mode_ = DATASET;
  std::string mapsetname = actionArgs.GetStringKey("data");
  std::string maskExpr = actionArgs.GetStringKey("mask");
  if (mapsetname.empty() && maskExpr.empty()) {
    mprinterr("Error: 'data' or 'mask' must be specified.\n");
    return Action::ERR;
  } else if (!mapsetname.empty() && !maskExpr.empty()) {
    mprinterr("Error: Specify either 'data' or 'mask'.\n");
    return Action::ERR;
  }
  if (!maskExpr.empty()) {
    mode_ = BY_MASK;
    startIdx_ = actionArgs.getKeyInt("toindex", 0);
  }
  topWriter_.InitTopWriter(actionArgs, "re-mapped", debugIn);
  // Get dataset
  DataSet* mapset = 0;
  if (mode_ == DATASET) {
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
    // User atom #s start from 1
    for (unsigned int i = 0; i != ds.Size(); i++)
      Map_.push_back( (int)ds.Dval(i) - 1 );
  } else if (mode_ == BY_MASK) {
    if (mask_.SetMaskString(maskExpr)) return Action::ERR;
  } else {
    mprinterr("Internal Error: Unhandled remap mode.\n");
    return Action::ERR;
  }

  mprintf("    REMAP:\n");
  if (mode_ == DATASET) {
    mprintf("\tRemapping atoms according to positions specified by data set '%s' (%zu atoms).\n",
            mapset->legend(), Map_.size());
  } else if (mode_ == BY_MASK) {
    mprintf("\tRemapping atoms in mask '%s' to starting position %i\n",
            mask_.MaskString(), startIdx_+1);
  }
  topWriter_.PrintOptions();
  return Action::OK;
}

// Action_Remap::Setup()
Action::RetType Action_Remap::Setup(ActionSetup& setup) {
  // Set up mask if needed
  if (mode_ == BY_MASK) {
    if (setup.Top().SetupIntegerMask(mask_)) return Action::ERR;
    mask_.MaskInfo();
    if (mask_.None()) {
      mprintf("Warning: No atoms selected, skipping.\n");
      return Action::SKIP;
    }
    Map_.assign(setup.Top().Natom(), -1);
    std::vector<bool> used(setup.Top().Natom(), false);
    // Put atoms in mask starting from position.
    int idx = startIdx_;
    for (AtomMask::const_iterator at = mask_.begin(); at != mask_.end(); ++at, ++idx)
    {
      Map_[idx] = *at;
      used[*at] = true;
    }
    // Find the first unused index
    idx = 0;
    for (; idx != setup.Top().Natom(); idx++)
      if (!used[idx]) break;
    // Fill in the rest of the map
    for (unsigned int n = 0; n != Map_.size(); ++n)
    {
      if (Map_[n] == -1) {
        Map_[n] = idx;
        used[idx] = true;
        for (; idx != setup.Top().Natom(); idx++)
          if (!used[idx]) break;
      }
    }
  }

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
