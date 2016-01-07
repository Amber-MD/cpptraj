#include "ReferenceAction.h"
#include "CpptrajStdio.h"
#include "DataSet_Coords_TRJ.h"
#include "ReferenceFrame.h"

ReferenceAction::ReferenceAction() :
  refMode_( FIRST ),
  refCrd_( 0 ),
  traj_( 0 ),
  previous_( false ),
  needsSetup_( false),
  fitRef_( false ),
  useMass_( false )
{}

ReferenceAction::~ReferenceAction() { if (traj_ != 0) delete traj_; }

// ReferenceAction::RefModeString()
std::string ReferenceAction::RefModeString() const {
  std::string modeString;
  if (previous_)
    modeString = "previous frame";
  else if (refMode_ == FIRST)
    modeString = "first frame";
  else if (refMode_ == TRAJ)
    modeString = "trajectory " + refCrd_->Meta().Legend();
  else // FRAME
    modeString = "\"" + refCrd_->Meta().Legend() + "\"";
  modeString += " (" + refMask_.MaskExpression() + ")";
  return modeString;
}

std::string ReferenceAction::help_ =
  "\t[ first | " + std::string(DataSetList::RefArgs) + " | previous |\n" +
  "\t  reftraj <filename> [" + std::string(DataSetList::TopArgs) + "] ]\n";

// ReferenceAction::InitRef()
int ReferenceAction::InitRef(ArgList& argIn, DataSetList const& DSLin,
                             bool fitRefIn, bool useMassIn)
{
  fitRef_ = fitRefIn;
  useMass_ = useMassIn;
  previous_ = argIn.hasKey("previous");
  // Attempt to determine reference mode.
  refMode_ = FIRST; // Default
  if (!argIn.hasKey("first")) {
    // Attempt to set refMode_ and refCrd_
    if (argIn.Contains("reftraj")) {
      // Reference trajectory. First try to find COORDS or TRAJ.
      std::string reftraj = argIn.GetStringKey("reftraj");
      refCrd_ = (DataSet_Coords*)DSLin.FindSetOfType( reftraj, DataSet::COORDS );
      if (refCrd_ == 0)
        refCrd_ = (DataSet_Coords*)DSLin.FindSetOfType( reftraj, DataSet::TRAJ );
      if (refCrd_ == 0) {
        mprintf("\tLoading reference trajectory '%s'\n", reftraj.c_str());
        // No existing COORDS set in DSL. Try to load as a TRAJ set.
        DataSet_Coords_TRJ* trj = new DataSet_Coords_TRJ();
        if (trj == 0) {
          mprinterr("Internal Error: Could not allocate memory for reftraj.\n");
          return 1;
        }
        Topology* RefParm = DSLin.GetTopology( argIn );
        if (RefParm == 0) {
          mprinterr("Error: No topology found for reftraj %s. Ensure topologies are loaded.\n",
                  reftraj.c_str());
          return 1;
        }
        if (trj->AddSingleTrajin( reftraj, argIn, RefParm )) return 1;
        trj->SetMeta( reftraj );
        traj_ = (DataSet_Coords*)trj;
        refCrd_ = traj_;
      } else
        mprintf("\tUsing set '%s' as reference trajectory.\n", refCrd_->legend());
      refMode_ = TRAJ;
    } else {
      // Look for existing reference structure in DSL.
      ReferenceFrame REF = DSLin.GetReferenceFrame( argIn );
      if (REF.error()) return 1;
      if (!REF.empty()) {
        refCrd_ = (DataSet_Coords*)REF.RefPtr();
        refMode_ = FRAME;
      }
    }
    /// Allocate space for reading Frame from COORDS if necessary
    if (refCrd_ != 0) {
      refFrame_ = refCrd_->AllocateFrame();
      needsSetup_ = true;
    }
  }
  return 0;
}

// ReferenceAction::SetupRefMask()
int ReferenceAction::SetupRefMask(Topology const& topIn) {
  if (topIn.SetupIntegerMask( refMask_ )) return 1;
  mprintf("\tReference mask:");
  refMask_.BriefMaskInfo();
  mprintf("\n");
  if (refMask_.None()) {
    mprinterr("Error: No reference atoms selected for parm %s, [%s]\n",
              topIn.c_str(), refMask_.MaskString());
    return 1;
  }
  selectedRef_.SetupFrameFromMask( refMask_, topIn.Atoms() );
  return 0;
}

// ReferenceAction::SelectRefAtoms()
void ReferenceAction::SelectRefAtoms(Frame const& frameIn) {
  refFrame_ = frameIn;
# ifdef MPI
  // Ensure all threads are using the same reference
  if (trajComm_.Master()) { // TODO MasterBcast
    rprintf("DEBUG: Sending reference frame to children.\n");
    for (int rank = 1; rank != trajComm_.Size(); rank++)
      refFrame_.SendFrame(rank, trajComm_);
  } else {
    rprintf("DEBUG: Receiving reference frame from master.\n");
    refFrame_.RecvFrame(0, trajComm_);
  }
  trajComm_.Barrier();
# endif
  selectedRef_.SetCoordinates( refFrame_, refMask_ );
  if (fitRef_)
    refTrans_ = selectedRef_.CenterOnOrigin( useMass_ );
}

// ReferenceAction::SetupRef()
int ReferenceAction::SetupRef(Topology const& topIn, int Ntgt) {
  if (refMode_ == FIRST) {
    if ( SetupRefMask( topIn ) ) return 1;
  } else if (previous_) {
    mprintf("Warning: %s: 'previous' may not work properly for changing topologies.\n");
    if ( SetupRefMask( topIn ) ) return 1;
  } else if (needsSetup_) { // FRAME, TRAJ; only set up once.
    if ( SetupRefMask( refCrd_->Top() ) ) return 1;
    if (refMode_ == FRAME)
      SelectRefAtoms( ((DataSet_Coords_REF*)refCrd_)->RefFrame() );
    needsSetup_ = false;
  }
  // Check that num atoms in target mask from this parm match ref parm mask
  if ( refMask_.Nselected() != Ntgt ) {
    mprintf("Warning: Number of atoms in target mask (%i) does not equal\n"
            "Warning:   number of atoms in reference mask (%i).\n",
            Ntgt, refMask_.Nselected());
    return 1;
  }
  return 0;
}
