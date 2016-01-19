#include "Action_Average.h"
#include "CpptrajStdio.h"
#include "Trajout_Single.h"
#include "DataSet_Coords_REF.h"

// CONSTRUCTOR
Action_Average::Action_Average() :
  debug_(0),
  Nframes_(0),
  crdset_(0)
{} 

void Action_Average::Help() const {
  mprintf("\t{crdset <set name> | <filename>} [<mask>]\n\t%s\n\t[TRAJOUT ARGS]\n"
          "  Calculate the average structure of atoms in <mask> over specified input frames.\n"
          "  If 'crdset' is specified a reference COORDS data set will be created with name\n"
          "  <set name>, otherwise the averaged coords will be written to <filename>.\n",
          ActionFrameCounter::HelpText);
}

// Action_Average::Init()
Action::RetType Action_Average::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
  // Get Keywords
  std::string crdName = actionArgs.GetStringKey("crdset");
  std::string avgfilename;
  if (crdName.empty()) {
    // Create output trajectory
    crdset_ = 0;
    avgfilename = actionArgs.GetStringNext();
    if (avgfilename.empty()) {
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

  // Initialize output trajectory from remaining arguments.
  if (crdset_ == 0) {
    outtraj_.SetDebug( debug_ );
    if (outtraj_.InitEnsembleTrajWrite(avgfilename, actionArgs.RemainingArgs(),
                                       TrajectoryFile::UNKNOWN_TRAJ, init.DSL().EnsembleNum()))
      return Action::ERR;
  }

  mprintf("    AVERAGE: Averaging over coordinates in mask [%s]\n",Mask1_.MaskString());
  FrameCounterInfo();
  if (crdset_ == 0)
    mprintf("\tWriting averaged coords to file '%s'\n", outtraj_.Traj().Filename().full());
  else
    mprintf("\tSaving averaged coords to set '%s'\n", crdset_->legend());

  Nframes_ = 0;

  return Action::OK;
}

// Action_Average::Setup()
/** On first call, set up Frame according to first topology. This will be
  * used for coordinate output.
  * On subsequent calls, determine if the number of atoms is greater than or
  * less than the original # atoms and reallocate if necessary.
  */
Action::RetType Action_Average::Setup(ActionSetup& setup) {
  if ( setup.Top().SetupIntegerMask( Mask1_ ) ) return Action::ERR;
  if (Mask1_.None()) {
    mprinterr("Warning: Cannot create average: No Atoms in mask.\n");
    return Action::SKIP;
  }
  Mask1_.MaskInfo();
  // Frame allocation/reallocation
  if (AvgFrame_.empty()) {
    mprintf("\tAveraging over %i atoms.\n", Mask1_.Nselected());
    AvgFrame_.SetupFrame( Mask1_.Nselected() );
    AvgFrame_.ZeroCoords();
  } else {
    // If the frame is already set up, check to see if the current number
    // of atoms is bigger or smaller. 
    if (Mask1_.Nselected() > AvgFrame_.Natom()) {
      // More atoms; reallocate the frame and copy existing coords.
      mprintf("Warning: Topology '%s' selected # atoms (%i) > original topology '%s' (%i)\n",
              setup.Top().c_str(), Mask1_.Nselected(), AvgParm_.c_str(), AvgFrame_.Natom());
      Frame tmp = AvgFrame_;
      AvgFrame_.SetupFrame( Mask1_.Nselected() );
      AvgFrame_.ZeroCoords();
      std::copy( tmp.xAddress(), tmp.xAddress() + tmp.size(), AvgFrame_.xAddress() );
    } else if (Mask1_.Nselected() < AvgFrame_.Natom()) {
      // Fewer atoms; just print a warning.
      mprintf("Warning: Topology '%s' selected # atoms (%i) < original topology '%s' (%i)\n",
              setup.Top().c_str(), Mask1_.Nselected(), AvgParm_.c_str(), AvgFrame_.Natom());
    }
    mprintf("\t%i atoms will be averaged for '%s'.\n", AvgFrame_.Natom(), setup.Top().c_str());
  }
  // AvgParm will be used for coordinate output.
  if (AvgParm_.Natom() < AvgFrame_.Natom()) {
    // If the number of selected atoms is less than the current parm, strip
    // the parm for output purposes.
    if (Mask1_.Nselected() < setup.Top().Natom()) {
      mprintf("Warning: Atom selection < total # atoms, stripping parm for averaging only:\n");
      Topology* aparm = setup.Top().modifyStateByMask( Mask1_ );
      if (aparm == 0) return Action::ERR;
      AvgParm_ = *aparm;
      delete aparm;
      if (debug_ > 0)
        AvgParm_.Summary();
    } else 
      AvgParm_ = setup.Top();
  }
  return Action::OK;  
}

// Action_Average::DoAction()
Action::RetType Action_Average::DoAction(int frameNum, ActionFrame& frm) {
  if ( CheckFrameCounter( frm.TrajoutNum() ) ) return Action::OK;

  if (AvgFrame_.AddByMask(frm.Frm(), Mask1_)) return Action::ERR;
  ++Nframes_; 

  return Action::OK;
}

#ifdef MPI
int Action_Average::SyncAction(Parallel::Comm const& commIn) {
  int total_frames = 0;
  commIn.Reduce( &total_frames, &Nframes_, 1, MPI_INT, MPI_SUM );
  AvgFrame_.SumToMaster(commIn);
  if (commIn.Master())
    Nframes_ = total_frames;
  return 0;
}
#endif

// Action_Average::Print()
void Action_Average::Print() {
  if (Nframes_ < 1) return;
  AvgFrame_.Divide( (double)Nframes_ );
  // NOTE: AvgFrame_ has coordinates only, blank CoordinateInfo is fine.
  mprintf("    AVERAGE: %i frames,", Nframes_);
  if (crdset_ == 0) {
    if (outtraj_.SetupTrajWrite(&AvgParm_, CoordinateInfo(), 1)) {
      mprinterr("Error: AVERAGE: Could not set up %s for write.\n",
                 outtraj_.Traj().Filename().full());
      return;
    }
    outtraj_.PrintInfo(0);
    outtraj_.WriteSingle(0, AvgFrame_);
    outtraj_.EndTraj();
  } else {
    mprintf(" COORDS set '%s'\n", crdset_->legend());
    DataSet_Coords_REF& ref = static_cast<DataSet_Coords_REF&>( *crdset_ );
    ref.CoordsSetup( AvgParm_, CoordinateInfo() ); // Coords Only
    ref.AddFrame( AvgFrame_ );
  }
}
