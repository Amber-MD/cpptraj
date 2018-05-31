#include "Action_CreateCrd.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // ByteString

Action_CreateCrd::Action_CreateCrd() : pindex_(-1), check_(true) {}

void Action_CreateCrd::Help() const {
  mprintf("\t[<name>] [ parm <name> | parmindex <#> ] [nocheck]\n"
          "  Create a COORDS data set named <name> for frames associated with the\n"
          "  specified topology.\n");
}

Action::RetType Action_CreateCrd::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Keywords
  Topology* parm = init.DSL().GetTopology( actionArgs );
  if (parm == 0) {
    mprinterr("Error: createcrd: No parm files loaded.\n");
    return Action::ERR;
  }
  pindex_ = parm->Pindex();
  check_ = !actionArgs.hasKey("nocheck");
  // DataSet
  std::string setname = actionArgs.GetStringNext();
  bool append = false;
  coords_ = 0;
  if (setname == "_DEFAULTCRD_") {
    // Special case: Creation of COORDS DataSet has been requested by an
    //               analysis and should already be present.
    coords_ = (DataSet_Coords_CRD*)init.DSL().FindSetOfType(setname, DataSet::COORDS);
  } else {
    if (!setname.empty()) {
      DataSet* ds = init.DSL().FindSetOfType( setname, DataSet::COORDS );
      if (ds != 0) {
        append = true;
        coords_ = (DataSet_Coords_CRD*)ds;
        pindex_ = coords_->Top().Pindex();
      }
    }
    if (coords_ == 0)
      coords_ = (DataSet_Coords_CRD*)init.DSL().AddSet(DataSet::COORDS, setname, "CRD");
    if (coords_ == 0) return Action::ERR;
  }
  // Do not set topology here since it may be modified later.
  if (append)
    mprintf("    CREATECRD: Appending coordinates to \"%s\"\n", coords_->legend());
  else
    mprintf("    CREATECRD: Saving coordinates from Top %s to \"%s\"\n",
            parm->c_str(), coords_->legend());
  if (!check_)
    mprintf("\tNot strictly enforcing that all frames have same # atoms.\n");
# ifdef MPI
  if (init.TrajComm().Size() > 1)
    mprintf("Warning: Synchronization of COORDS data sets over multiple threads is\n"
            "Warning:   experimental and may be slower than reading in via a single\n"
            "Warning:   thread. Users are encouraged to run benchmarks before\n"
            "Warning:   extensive usage.\n");
# endif
  return Action::OK;
}

Action::RetType Action_CreateCrd::Setup(ActionSetup& setup) {
  // Set COORDS topology now if not already set.
  if (setup.Top().Pindex() == pindex_ && coords_->Top().Natom() == 0) {
    coords_->CoordsSetup( setup.Top(), setup.CoordInfo() );
    // Estimate memory usage
    mprintf("\tEstimated memory usage (%i frames): %s\n",
            setup.Nframes(),
            ByteString(coords_->SizeInBytes(setup.Nframes()), BYTE_DECIMAL).c_str());
  }
  // If # atoms in currentParm does not match coords, warn user.
  if (setup.Top().Natom() != coords_->Top().Natom()) {
    if (check_) {
      mprinterr("Error: # atoms in current topology (%i) != # atoms in coords set \"%s\" (%i)\n",
                setup.Top().Natom(), coords_->legend(), coords_->Top().Natom());
      return Action::ERR;
    } else {
      mprintf("Warning: # atoms in current topology (%i) != # atoms in coords set \"%s\" (%i)\n"
              "Warning:   The resulting COORDS data set may have problems.\n",
              setup.Top().Natom(), coords_->legend(), coords_->Top().Natom());
    }
  }
  return Action::OK;
}

Action::RetType Action_CreateCrd::DoAction(int frameNum, ActionFrame& frm) 
{
  coords_->AddFrame( frm.Frm() );
  return Action::OK;
}
