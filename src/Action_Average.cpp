#include "Action_Average.h"
#include "CpptrajStdio.h"
#include "Trajout_Single.h"
#include "DataSet_Coords_REF.h"

// CONSTRUCTOR
Action_Average::Action_Average() :
  ensembleNum_(-1),
  debug_(0),
  AvgFrame_(0),
  Natom_(0),
  Nframes_(0),
  crdset_(0)
{ } 

void Action_Average::Help() const {
  mprintf("\t{crdset <set name> | <filename>} [<mask>]\n\t%s\n\t[TRAJOUT ARGS]\n"
          "  Calculate the average structure of atoms in <mask> over specified input frames.\n"
          "  If 'crdset' is specified a reference COORDS data set will be created with name\n"
          "  <set name>, otherwise the averaged coords will be written to <filename>.\n",
          ActionFrameCounter::HelpText);
}

// DESTRUCTOR
Action_Average::~Action_Average() {
  //fprintf(stderr,"Average Destructor.\n");
  if (AvgFrame_!=0) delete AvgFrame_;
}

// Action_Average::Init()
Action::RetType Action_Average::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  ensembleNum_ = init.DSL().EnsembleNum();
  debug_ = debugIn;
  // Get Keywords
  std::string crdName = actionArgs.GetStringKey("crdset");
  if (crdName.empty()) {
    crdset_ = 0;
    avgfilename_ = actionArgs.GetStringNext();
    if (avgfilename_.empty()) {
      mprinterr("Error: average: No filename given.\n");
      return Action::ERR;
    }
  } else {
    // Create REF_FRAME data set.
    crdset_ = init.DSL().AddSet(DataSet::REF_FRAME, crdName, "AVGSTRUCT");
    if (crdset_ == 0) {
      mprinterr("Error: Could not allocate average coordinate data set '%s'\n", crdName.c_str());
      return Action::ERR;
    }
#   ifdef MPI
    // crdset_ is not written to until Print(), no sync needed.
    crdset_->SetNeedsSync( false );
#   endif
  }
  // Get start/stop/offset args
  if (InitFrameCounter(actionArgs)) return Action::ERR;

  // Get Masks
  Mask1_.SetMaskString( actionArgs.GetMaskNext() );

  // Save all remaining arguments for setting up the trajectory at the end.
  if (crdset_ == 0)
    trajArgs_ = actionArgs.RemainingArgs();

  mprintf("    AVERAGE: Averaging over coordinates in mask [%s]\n",Mask1_.MaskString());
  FrameCounterInfo();
  if (crdset_ == 0)
    mprintf("\tWriting averaged coords to file '%s'\n",avgfilename_.c_str());
  else
    mprintf("\tSaving averaged coords to set '%s'\n", crdset_->legend());

  Nframes_ = 0;

  return Action::OK;
}

// Action_Average::Setup()
/** On first call, set up Frame according to first parmtop. This will be
  * used for coordinate output.
  * On subsequent calls, determine if the number of atoms is greater than or
  * less than the original # atoms. Never calculate more than the original
  * # atoms.
  */
Action::RetType Action_Average::Setup(ActionSetup& setup) {

  if ( setup.Top().SetupIntegerMask( Mask1_ ) ) return Action::ERR;

  if (Mask1_.None()) {
    mprinterr("Warning: Cannot create average: No Atoms in mask.\n");
    return Action::SKIP;
  }

  if (AvgFrame_==0) {
    mprintf("\tAveraging over %i atoms.\n",Mask1_.Nselected());
    AvgFrame_ = new Frame(Mask1_.Nselected());
    AvgFrame_->ZeroCoords();
    Natom_ = AvgFrame_->Natom(); // Initiially equal to Mask1.Nselected
    // AvgParm will be used for coordinate output
    // If the number of selected atoms is less than the current parm, strip
    // the parm for output purposes.
    if (Mask1_.Nselected()<setup.Top().Natom()) {
      mprintf("Warning: Atom selection < natom, stripping parm for averaging only:\n");
      Topology* aparm = setup.Top().modifyStateByMask(Mask1_);
      if (aparm == 0) return Action::ERR;
      AvgParm_ = *aparm;
      delete aparm;
      if (debug_ > 0)
        AvgParm_.Summary();
    } else 
      AvgParm_ = setup.Top();
  } else {
    // If the frame is already set up, check to see if the current number
    // of atoms is bigger or smaller. If bigger, only average Natom coords.
    // If smaller, only average P->natom coords.
    if (Mask1_.Nselected() > AvgFrame_->Natom()) {
      Natom_ = AvgFrame_->Natom();
      mprintf("Warning: Parm '%s' selected # atoms (%i) > original parm '%s'\n",
              setup.Top().c_str(), Mask1_.Nselected(), AvgParm_.c_str());
      mprintf("Warning:   selected# atoms (%i).\n",AvgFrame_->Natom());
    } else if (Mask1_.Nselected() < AvgFrame_->Natom()) {
      Natom_ = Mask1_.Nselected();
      mprintf("Warning: Parm '%s' selected # atoms (%i) < original parm '%s'\n",
              setup.Top().c_str(), Mask1_.Nselected(), AvgParm_.c_str());
      mprintf("Warning:   selected # atoms (%i).\n",AvgFrame_->Natom());
    } else {
      Natom_ = AvgFrame_->Natom();
    }
    mprintf("\t%i atoms will be averaged for '%s'.\n",Natom_, setup.Top().c_str());
  }
        
  return Action::OK;  
}

// Action_Average::DoAction()
Action::RetType Action_Average::DoAction(int frameNum, ActionFrame& frm) {
  if ( CheckFrameCounter( frm.TrajoutNum() ) ) return Action::OK;

  if (AvgFrame_->AddByMask(frm.Frm(), Mask1_)) return Action::ERR;
  ++Nframes_; 

  return Action::OK;
}

#ifdef MPI
int Action_Average::SyncAction(Parallel::Comm const& commIn) {
  int total_frames = 0;
  commIn.Reduce( &total_frames, &Nframes_, 1, MPI_INT, MPI_SUM );
  AvgFrame_->SumToMaster(commIn);
  if (commIn.Master())
    Nframes_ = total_frames;
  return 0;
}
#endif


// Action_Average::Print()
void Action_Average::Print() {
  if (Nframes_ < 1) return;
  double d_Nframes = (double) Nframes_;
  AvgFrame_->Divide(d_Nframes);
  // NOTE: AvgFrame_ has coordinates only, blank CoordinateInfo is fine.
  mprintf("    AVERAGE: %i frames,", Nframes_);
  if (crdset_ == 0) {
    Trajout_Single outfile;
    mprintf(" [%s %s]\n",avgfilename_.c_str(), trajArgs_.ArgLine());
    if (outfile.PrepareEnsembleTrajWrite(avgfilename_, trajArgs_, &AvgParm_,
                                         CoordinateInfo(), 1, 
                                         TrajectoryFile::UNKNOWN_TRAJ, ensembleNum_)) 
    {
      mprinterr("Error: AVERAGE: Could not set up %s for write.\n",avgfilename_.c_str());
      return;
    }
    outfile.PrintInfo(0);
    outfile.WriteSingle(0, *AvgFrame_);
    outfile.EndTraj();
  } else {
    mprintf(" COORDS set '%s'\n", crdset_->legend());
    DataSet_Coords_REF& ref = static_cast<DataSet_Coords_REF&>( *crdset_ );
    ref.CoordsSetup( AvgParm_, CoordinateInfo() ); // Coords Only
    ref.AddFrame( *AvgFrame_ );
  }
}
