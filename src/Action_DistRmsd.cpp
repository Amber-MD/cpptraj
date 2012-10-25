// DISTRMSD
#include "Action_DistRmsd.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_DistRmsd::Action_DistRmsd() :
  drmsd_(NULL),
  refmode_(UNKNOWN_REF)
{}

void Action_DistRmsd::Help() {
  mprintf("drmsd <name> <mask> [<refmask>] [out filename]\n");
  mprintf("      [ first | ref <filename> | refindex <#> |\n");
  mprintf("      [ reftraj <filename> [parm <parmname> | parmindex <#>] ]\n"); 
}

/** Setup RefMask based on given Topology. Allocate space for selected
  * reference atoms. 
  */
int Action_DistRmsd::SetRefMask( Topology* RefParm ) {
  if (RefParm == NULL) return 1;
  if (RefParm->SetupIntegerMask( RefMask_ )) return 1;
  if (RefMask_.None()) {
    mprinterr("Error: distrmsd: No reference atoms selected for parm %s, [%s]\n",
              RefParm->c_str(), RefMask_.MaskString());
    return 1;
  }
  SelectedRef_.SetupFrameFromMask( RefMask_, RefParm->Atoms() );
  return 0;
}

/** Setup selected reference coordinates based on given frame
  * and RefMask. 
  */
void Action_DistRmsd::SetRefStructure( Frame& frameIn ) {
  RefFrame_ = frameIn;
  SelectedRef_.SetCoordinates( RefFrame_, RefMask_ );
}

// Action_DistRmsd::init()
/** Called once before traj processing. Set up reference info. */
Action::RetType Action_DistRmsd::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  std::string refname, reftrajname;
  int refindex = -1;
  Topology* RefParm = NULL;

  // Check for keywords
  ArgList::ConstArg rmsdFile = actionArgs.getKeyString("out");
  // Reference keywords
  refmode_ = UNKNOWN_REF;
  if ( actionArgs.hasKey("first") ) {
    refmode_ = FIRST;
  } else {
    refindex = actionArgs.getKeyInt("refindex", -1);
    if (actionArgs.hasKey("reference")) refindex = 0;
    refname = actionArgs.GetStringKey("ref");
    if (refindex==-1 && refname.empty()) {
      reftrajname = actionArgs.GetStringKey("reftraj");
      if (!reftrajname.empty()) {
        RefParm = PFL->GetParm( actionArgs );
        refmode_ = REFTRAJ;
      }
    } else
      refmode_ = REF;
  }

  // Get the RMS mask string for target 
  ArgList::ConstArg mask0 = actionArgs.getNextMask();
  TgtMask_.SetMaskString(mask0);
  // Get the RMS mask string for reference
  ArgList::ConstArg mask1 = actionArgs.getNextMask();
  if (mask1==NULL)
    mask1 = mask0;
  RefMask_.SetMaskString( mask1 );

  // Check that a reference keyword was specified. If not, default to first.
  if (refmode_==UNKNOWN_REF) {
    mprintf("Warning: distrmsd: No reference structure given. Defaulting to first.\n");
    refmode_ = FIRST;
  }
  // Initialize reference if not 'first'.
  if (refmode_ != FIRST) {
    if ( !reftrajname.empty() ) {
      // Reference trajectory
      if (RefParm == NULL) {
        mprinterr("Error: distrmsd: Could not get parm for reftraj %s\n", reftrajname.c_str());
        return Action::ERR;
      }
      if (SetRefMask( RefParm )!=0) return Action::ERR;
      // Attempt to open reference traj.
      if (RefTraj_.SetupTrajRead( reftrajname, NULL, RefParm)) {
        mprinterr("Error: distrmsd: Could not set up reftraj %s\n", reftrajname.c_str());
        return Action::ERR;
      }
      RefFrame_.SetupFrameV(RefParm->Atoms(), RefTraj_.HasVelocity());
      if (RefTraj_.BeginTraj(false)) {
        mprinterr("Error: distrmsd: Could not open reftraj %s\n", reftrajname.c_str());
        return Action::ERR;
      }
    } else {
      // Reference by name/tag
      if (!refname.empty())
        refindex = FL->FindName( refname );
      // Get reference by index
      Frame* TempFrame = FL->GetFrame( refindex );
      // Get parm for reference
      RefParm = FL->GetFrameParm( refindex );
      if (RefParm == NULL) {
        mprinterr("Error: distrmsd: Could not get parm for frame %s\n", FL->FrameName(refindex));
        return Action::ERR;
      }
      if (SetRefMask( RefParm )!=0) return Action::ERR;
      SetRefStructure( *TempFrame );
    }
  }

  // Set up the RMSD data set
  drmsd_ = DSL->Add(DataSet::DOUBLE, actionArgs.getNextString(),"DRMSD");
  if (drmsd_==NULL) return Action::ERR;
  // Add dataset to data file list
  DFL->Add(rmsdFile,drmsd_);

  mprintf("    DISTRMSD: (%s), reference is",TgtMask_.MaskString());
  if (refmode_ == FIRST)
    mprintf(" first frame");
  else if (refmode_==REFTRAJ)
    mprintf(" trajectory %s",RefTraj_.FullTrajStr());
  else {
    mprintf(" reference frame");
    if (!refname.empty())
      mprintf(" %s",refname.c_str());
    else
      mprintf(" index %i", refindex);
  }
  mprintf(" (%s)\n",RefMask_.MaskString());

  return Action::OK;
}

// Action_DistRmsd::setup()
/** Called every time the trajectory changes. Set up TgtMask for the new 
  * parmtop and allocate space for selected atoms from the Frame.
  */
Action::RetType Action_DistRmsd::Setup(Topology* currentParm, Topology** parmAddress) {

  if ( currentParm->SetupIntegerMask(TgtMask_) ) return Action::ERR;
  if ( TgtMask_.None() ) {
    mprintf("    Error: DistRmsd::setup: No atoms in mask.\n");
    return Action::ERR;
  }
  // Allocate space for selected atoms in the frame. This will also put the
  // correct masses in based on the mask.
  SelectedTgt_.SetupFrameFromMask(TgtMask_, currentParm->Atoms());

  // Reference setup if 'first'
  if (refmode_ == FIRST) {
    if ( SetRefMask( currentParm )!=0 ) return Action::ERR;
  }

  // Check that num atoms in frame mask from this parm match ref parm mask
  if ( RefMask_.Nselected() != TgtMask_.Nselected() ) {
    mprintf( "    Error: Number of atoms in RMS mask (%i) does not \n",TgtMask_.Nselected());
    mprintf( "           equal number of atoms in Ref mask (%i).\n",RefMask_.Nselected());
    return Action::ERR;
  }

  return Action::OK;
}

// Action_DistRmsd::action()
/** Called every time a frame is read in. Calc distance RMSD.
  * If first is true, set the first frame read in as reference.
  */
Action::RetType Action_DistRmsd::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  // Perform any needed reference actions
  if (refmode_ == FIRST) {
    SetRefStructure( *currentFrame );
    refmode_ = REF;
  } else if (refmode_ == REFTRAJ) {
    RefTraj_.GetNextFrame( RefFrame_ );
    SelectedRef_.SetCoordinates(RefFrame_, RefMask_);
  }

  // Set selected frame atoms. Masses have already been set.
  SelectedTgt_.SetCoordinates(*currentFrame, TgtMask_);

  // DEBUG
/*  mprintf("  DEBUG: RMSD: First atom coord in SelectedTgt is : "); 
  SelectedTgt->printAtomCoord(0);
  mprintf("  DEBUG: RMSD: First atom coord in SelectedRef is : ");
  SelectedRef->printAtomCoord(0);
*/

  double DR = SelectedTgt_.DISTRMSD( SelectedRef_ );

  drmsd_->Add(frameNum, &DR);

  return Action::OK;
}

