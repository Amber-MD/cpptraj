// Action_Average
#include "Action_Average.h"
#include "CpptrajStdio.h"
#include "Trajout.h"

// CONSTRUCTOR
Action_Average::Action_Average() :
  debug_(0),
  AvgFrame_(0),
  AvgParm_(0),
  parmStripped_(false),
  Natom_(0),
  Nframes_(0)
{ } 

void Action_Average::Help() {
  mprintf("\t<filename> [<mask>] %s\n", ActionFrameCounter::HelpText);
  mprintf("\t[TRAJOUT ARGS]\n");
  mprintf("\tCalculate the average structure of atoms in <mask> over specified input frames.\n");
}

// DESTRUCTOR
Action_Average::~Action_Average() {
  //fprintf(stderr,"Average Destructor.\n");
  if (AvgFrame_!=0) delete AvgFrame_;
  if (parmStripped_ && AvgParm_!=0) delete AvgParm_;
}

// Action_Average::init()
Action::RetType Action_Average::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  debug_ = debugIn;
  // Get Keywords
  avgfilename_ = actionArgs.GetStringNext();
  if (avgfilename_.empty()) {
    mprinterr("Error: average: No filename given.\n");
    return Action::ERR;
  }
  // Get start/stop/offset args
  if (InitFrameCounter(actionArgs)) return Action::ERR;

  // Get Masks
  Mask1_.SetMaskString( actionArgs.GetMaskNext() );

  // Save all remaining arguments for setting up the trajectory at the end.
  trajArgs_ = actionArgs.RemainingArgs();

  mprintf("    AVERAGE: Averaging over coordinates in mask [%s]\n",Mask1_.MaskString());
  FrameCounterInfo();
  mprintf("\tWriting averaged coords to [%s]\n",avgfilename_.c_str());

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

  if (AvgFrame_==0) {
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
      if (debug_ > 0)
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
Action::RetType Action_Average::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) 
{
  if ( CheckFrameCounter( frameNum ) ) return Action::OK;

  AvgFrame_->AddByMask(*currentFrame, Mask1_);
  ++Nframes_; 

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
