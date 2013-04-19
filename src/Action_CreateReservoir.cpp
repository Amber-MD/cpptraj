#include "Action_CreateReservoir.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_CreateReservoir::Action_CreateReservoir() :
  original_trajparm_(0),
  ene_(0),
  bin_(0),
  trajIsOpen_(false),
  nframes_(0)
{}

void Action_CreateReservoir::Help() {
  mprintf("\t<filename> ene <energy data set> [bin <cluster bin data set>]\n");
  mprintf("\t[parm <parmfile> | parmindex <#>]\n");
}

// Action_CreateReservoir::Init()
Action::RetType Action_CreateReservoir::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get keywords
  filename_ = actionArgs.GetStringNext();
  if (filename_.empty()) {
    mprinterr("Error: createreservoir: No filename specified.\n");
    return Action::ERR;
  }
  // Get parm for reservoir traj
  original_trajparm_ = PFL->GetParm( actionArgs );
  if (original_trajparm_ == 0) {
    mprinterr("Error: createreservoir: no topology.\n");
    return Action::ERR;
  }
  // Get energy data set
  std::string ene_dsname = actionArgs.GetStringKey("ene");
  ene_ = DSL->GetDataSet( ene_dsname );
  if (ene_ == 0) {
    mprinterr("Error: could not get energy data set %s\n", ene_dsname.c_str());
    return Action::ERR;
  }
  if (ene_->Type() != DataSet::FLOAT && ene_->Type() != DataSet::DOUBLE) {
    mprinterr("Error: energy data set %s must be type FLOAT or DOUBLE.\n", ene_dsname.c_str());
    return Action::ERR;
  }
  // Get bin data set
  std::string bin_dsname = actionArgs.GetStringKey("bin");
  if (!bin_dsname.empty()) {
    bin_ = DSL->GetDataSet( bin_dsname );
    if (bin_ == 0) {
      mprinterr("Error: could not get bin data set %s\n", bin_dsname.c_str());
      return Action::ERR;
    }
    if (bin_->Type() != DataSet::INT) {
      mprinterr("Error: bin data set %s must be type INTEGER.\n", bin_dsname.c_str());
      return Action::ERR;
    }
  }
  trajIsOpen_ = false;
  nframes_ = 0;
  // Setup output reservoir file
  reservoir_.SetDebug( debugIn );
  // Set title
  reservoir_.SetTitle( actionArgs.GetStringKey("title") );
  // Process additional netcdf traj args
  reservoir_.processWriteArgs( actionArgs );

  mprintf("    CREATERESERVOIR: %s, energy data %s", filename_.c_str(),
          ene_->Legend().c_str());
  if (!bin_dsname.empty())
    mprintf(", bin data %s", bin_->Legend().c_str());
  mprintf("\n\tTopology: %s\n", original_trajparm_->c_str());
  return Action::OK;
}

// Action_CreateReservoir::Setup()
Action::RetType Action_CreateReservoir::Setup(Topology* currentParm, Topology** parmAddress)
{
  // Check that input parm matches current parm
  if (original_trajparm_->Pindex() != currentParm->Pindex()) {
    mprintf("Info: createreservoir was set up for topology %s\n", original_trajparm_->c_str());
    mprintf("Info: skipping topology %s\n", currentParm->c_str());
    return Action::ERR;
  }
  if (!trajIsOpen_) {
    mprintf("\tCreating reservoir file %s\n", filename_.c_str());
    // Use parm to set up box info for the reservoir. FIXME: Necessary?
    reservoir_.SetBox( currentParm->ParmBox() );
    // Set up write and open - no append.
    if (reservoir_.setupTrajout( filename_, currentParm, currentParm->Nframes(), false))
      return Action::ERR;
    trajIsOpen_ = true;
    nframes_ = 0;
  }
  return Action::OK;
}

// Action_CreateReservoir::DoAction()
Action::RetType Action_CreateReservoir::DoAction(int frameNum, Frame* currentFrame, 
                                                 Frame** frameAddress) 
{
  if (reservoir_.writeFrame(nframes_++, currentFrame->xAddress(), currentFrame->vAddress(),
                                        currentFrame->bAddress(), currentFrame->Temperature()))
    return Action::ERR;
  return Action::OK;
}

// Action_CreateReservoir::Print()
void Action_CreateReservoir::Print() {
  mprintf("\tReservoir %s: %u frames.\n", filename_.c_str(), nframes_);
  reservoir_.closeTraj();
  trajIsOpen_=false;
}
