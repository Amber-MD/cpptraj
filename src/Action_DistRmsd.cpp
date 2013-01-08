// DISTRMSD
#include "Action_DistRmsd.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_DistRmsd::Action_DistRmsd() :
  drmsd_(0),
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
  if (RefParm == 0) return 1;
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
  std::string reftrajname;
  Topology* RefParm = 0;
  ReferenceFrame REF;
  // Check for keywords
  DataFile* outfile = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  // Reference keywords
  refmode_ = UNKNOWN_REF;
  if ( actionArgs.hasKey("first") ) {
    refmode_ = FIRST;
  } else {
    REF = FL->GetFrame( actionArgs );
    if (REF.error()) return Action::ERR; 
    if (REF.empty()) {
      reftrajname = actionArgs.GetStringKey("reftraj");
      if (!reftrajname.empty()) {
        RefParm = PFL->GetParm( actionArgs );
        refmode_ = REFTRAJ;
      } else {
        // No reference keywords specified. Default to first.
        mprintf("Warning: distrmsd: No reference structure given. Defaulting to first.\n");
        refmode_ = FIRST;
      }
    } else
      refmode_ = REFFRAME;
  }
  // Get the RMS mask string for target 
  std::string mask0 = actionArgs.GetMaskNext();
  TgtMask_.SetMaskString(mask0);
  // Get the RMS mask string for reference
  std::string mask1 = actionArgs.GetMaskNext();
  if (mask1.empty())
    mask1 = mask0;
  RefMask_.SetMaskString( mask1 );

  // Initialize reference if not 'first'.
  if (refmode_ != FIRST) {
    if ( !reftrajname.empty() ) {
      // Reference trajectory
      if (RefParm == 0) {
        mprinterr("Error: distrmsd: Could not get parm for reftraj %s\n", reftrajname.c_str());
        return Action::ERR;
      }
      if (SetRefMask( RefParm )!=0) return Action::ERR;
      // Attempt to open reference traj.
      if (RefTraj_.SetupTrajRead( reftrajname, &actionArgs, RefParm)) {
        mprinterr("Error: distrmsd: Could not set up reftraj %s\n", reftrajname.c_str());
        return Action::ERR;
      }
      RefFrame_.SetupFrameV(RefParm->Atoms(), RefTraj_.HasVelocity());
      if (RefTraj_.BeginTraj(false)) {
        mprinterr("Error: distrmsd: Could not open reftraj %s\n", reftrajname.c_str());
        return Action::ERR;
      }
    } else {
      // Reference Frame
      if (SetRefMask( REF.Parm() ) != 0) return Action::ERR;
      SetRefStructure( *(REF.Coord()) );
    }
  }

  // Set up the RMSD data set
  drmsd_ = DSL->AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(),"DRMSD");
  if (drmsd_==0) return Action::ERR;
  // Add dataset to data file list
  if (outfile != 0) outfile->AddSet( drmsd_ );

  mprintf("    DISTRMSD: (%s), reference is",TgtMask_.MaskString());
  if (refmode_ == FIRST)
    mprintf(" first frame");
  else if (refmode_==REFTRAJ)
    mprintf(" trajectory %s",RefTraj_.FullTrajStr());
  else // REFFRAME
    mprintf(" reference frame %s", REF.FrameName());
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
    refmode_ = REFFRAME;
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

