// Action_Average
#include "Action_Average.h"
#include "CpptrajStdio.h"
#include "Trajout.h"

// CONSTRUCTOR
Action_Average::Action_Average() :
  AvgFrame_(NULL),
  AvgParm_(NULL),
  parmStripped_(false),
  Natom_(0),
  Nframes_(0),
  start_(0),
  stop_(0),
  offset_(1),
  targetFrame_(0)
{ } 

void Action_Average::Help() {
  mprintf("average <filename> [mask] [start <start>] [stop <stop>] [offset <offset>]\n");
  mprintf("        [TRAJOUT ARGS]\n");
  mprintf("\tCalculate the average structure over input frames.\n");
}

// DESTRUCTOR
Action_Average::~Action_Average() {
  //fprintf(stderr,"Average Destructor.\n");
  if (AvgFrame_!=NULL) delete AvgFrame_;
  if (parmStripped_ && AvgParm_!=NULL) delete AvgParm_;
}

// Action_Average::init()
Action::RetType Action_Average::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get Keywords
  avgfilename_ = actionArgs.GetStringNext();
  if (avgfilename_.empty()) {
    mprinterr("Error: average: No filename given.\n");
    return Action::ERR;
  }
  // TODO: Create frame tracker class for actions
  // User start/stop args are +1
  start_ = actionArgs.getKeyInt("start",1);
  --start_;
  targetFrame_ = start_;
  stop_ = actionArgs.getKeyInt("stop",-1);
  if (stop_!=-1) --stop_;
  offset_ = actionArgs.getKeyInt("offset",1);

  // Get Masks
  Mask1_.SetMaskString( actionArgs.getNextMask() );

  // Save all remaining arguments for setting up the trajectory at the end.
  ArgList::ConstArg arg;
  while ( (arg = actionArgs.getNextString() ) != NULL )
    trajArgs_.AddArg( arg );

  mprintf("    AVERAGE: Averaging over coordinates in mask [%s]",Mask1_.MaskString());
  if (stop_==-1) 
    mprintf(", starting from frame %i",start_+1);
  else
    mprintf(", frames %i-%i",start_+1,stop_+1);
  if (offset_!=1)
    mprintf(", offset %i",offset_);
  mprintf(".\n             Writing averaged coords to [%s]\n",avgfilename_.c_str());

  Nframes_ = 0;

  return Action::OK;
}

// Action_Average::setup()
/** On first call, set up Frame according to first parmtop. This will be
  * used for coordinate output.
  * On subsequent calls, determine if the number of atoms is greater than or
  * less than the original # atoms. Never calculate more than the original
  * # atoms.
  */
Action::RetType Action_Average::Setup(Topology* currentParm, Topology** parmAddress) {

  if ( currentParm->SetupIntegerMask( Mask1_ ) ) return Action::ERR;

  if (Mask1_.None()) {
    mprintf("    Error: Average::setup: No Atoms in mask.\n");
    return Action::ERR;
  }

  mprintf("    AVERAGE:");

  if (AvgFrame_==NULL) {
    mprintf(" Averaging over %i atoms.\n",Mask1_.Nselected());
    AvgFrame_ = new Frame(Mask1_.Nselected());
    AvgFrame_->ZeroCoords();
    Natom_ = AvgFrame_->Natom(); // Initiially equal to Mask1.Nselected
    // AvgParm will be used for coordinate output
    // If the number of selected atoms is less than the current parm, strip
    // the parm for output purposes.
    if (Mask1_.Nselected()<currentParm->Natom()) {
      mprintf("             Atom selection < natom, stripping parm for averaging only:\n");
      AvgParm_ = currentParm->modifyStateByMask(Mask1_);
      parmStripped_=true;
      AvgParm_->Summary();
    } else 
      AvgParm_ = currentParm; 
  } else {
    // If the frame is already set up, check to see if the current number
    // of atoms is bigger or smaller. If bigger, only average Natom coords.
    // If smaller, only average P->natom coords.
    if (Mask1_.Nselected() > AvgFrame_->Natom()) {
      Natom_ = AvgFrame_->Natom();
      mprintf("Warning: Average [%s]: Parm %s selected # atoms (%i) > original parm %s\n",
              avgfilename_.c_str(), currentParm->c_str(),
              Mask1_.Nselected(), AvgParm_->c_str());
      mprintf("         selected# atoms (%i).\n",AvgFrame_->Natom());
    } else if (Mask1_.Nselected() < AvgFrame_->Natom()) {
      Natom_ = Mask1_.Nselected();
      mprintf("Warning: Average[%s]: Parm %s selected # atoms (%i) < original parm %s\n",
              avgfilename_.c_str(), currentParm->c_str(), 
              Mask1_.Nselected(), AvgParm_->c_str());
      mprintf("         selected # atoms (%i).\n",AvgFrame_->Natom());
    } else {
      Natom_ = AvgFrame_->Natom();
    }
    mprintf("    AVERAGE: %i atoms will be averaged for this parm.\n",Natom_);
  }
        
  return Action::OK;  
}

// Action_Average::action()
Action::RetType Action_Average::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  if (frameNum != targetFrame_) return Action::OK;

  AvgFrame_->AddByMask(*currentFrame, Mask1_);
  ++Nframes_; 

  targetFrame_ += offset_;
  // Since frameNum will never be -1 this effectively disables the routine
  if (targetFrame_>stop_ && stop_!=-1) targetFrame_=-1;
  return Action::OK;
} 

// Action_Average::print()
void Action_Average::Print() {
  Trajout outfile;
  double d_Nframes;

  if (Nframes_ < 1) return;
  d_Nframes = (double) Nframes_;
  AvgFrame_->Divide(d_Nframes);

  mprintf("    AVERAGE: [%s %s]\n",avgfilename_.c_str(), trajArgs_.ArgLine());

  if (outfile.SetupTrajWrite(avgfilename_, &trajArgs_, AvgParm_, TrajectoryFile::UNKNOWN_TRAJ)) 
  {
    mprinterr("Error: AVERAGE: Could not set up %s for write.\n",avgfilename_.c_str());
    return;
  }

  outfile.PrintInfo(0);

  outfile.WriteFrame(0, AvgParm_, *AvgFrame_);

  outfile.EndTraj();
}
