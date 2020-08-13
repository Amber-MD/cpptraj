#include "ReferenceAction.h"
#include "DataSetList.h"
#include "CpptrajStdio.h"
#include "DataSet_Coords_TRJ.h"
#include "ReferenceFrame.h"
#include "ArgList.h"

ReferenceAction::ReferenceAction() :
  refMode_( FIRST ),
  refCrd_( 0 ),
  traj_( 0 ),
  needsSetup_( false),
  fitRef_( false ),
  useMass_( false )
{}

ReferenceAction::~ReferenceAction() { if (traj_ != 0) delete traj_; }

// ReferenceAction::RefModeString()
std::string ReferenceAction::RefModeString() const {
  std::string modeString;
  if (refMode_ == PREVIOUS)
    modeString = "previous frame";
  else if (refMode_ == FIRST)
    modeString = "first frame";
  else if (refMode_ == TRAJ)
    modeString = "trajectory " + refCrd_->Meta().Legend();
  else // FRAME
    modeString = "\"" + refCrd_->Meta().Legend() + "\"";
  if (refMask_.MaskStringSet())
    modeString += " (" + refMask_.MaskExpression() + ")";
  return modeString;
}

std::string ReferenceAction::help_ =
  "\t[ first | " + std::string(DataSetList::RefArgs) + " | previous |\n" +
  "\t  reftraj <name> [" + std::string(DataSetList::TopArgs) + "] ]\n";

#ifdef MPI
/** Should be called after InitRef(). */
int ReferenceAction::SetTrajComm( Parallel::Comm const& commIn ) {
  // Does not work for previous
  if (refMode_ == PREVIOUS && commIn.Size() > 1) {
    mprinterr("Error: 'previous' reference does not work in parallel.\n");
    return 1;
  }
  trajComm_ = commIn;
  // Sanity check. Ensure refMode_ matches on all threads.
  std::vector<int> all_modes( trajComm_.Size() );
  trajComm_.AllGather( &refMode_, 1, MPI_INT, &all_modes[0] );
  for (int rank = 1; rank < trajComm_.Size(); rank++)
    if (all_modes[rank] != all_modes[0]) { 
      mprinterr("Error: Reference mode on rank does not match master mode.\n");
      return 1;
    }
  if (refMode_ == FRAME) {
    // Ensure natom matches on all threads.
    int natom = ((DataSet_Coords_REF*)refCrd_)->RefFrame().Natom();
    trajComm_.AllGather( &natom, 1, MPI_INT, &all_modes[0] );
    for (int rank = 1; rank < trajComm_.Size(); rank++)
      if (all_modes[rank] != all_modes[0]) {
        mprinterr("Error: Reference mode # atoms on rank (%i) != # atoms on master (%i)\n",
                  all_modes[rank], all_modes[0]);
        return 1;
      }
  }
  return 0;
}
#endif

// ReferenceAction::InitRef()
int ReferenceAction::InitRef(ArgList& argIn, DataSetList const& DSLin,
                             bool fitRefIn, bool useMassIn)
{
  fitRef_ = fitRefIn;
  useMass_ = useMassIn;
  // Attempt to determine reference mode.
  refMode_ = FIRST; // Default
  if (argIn.hasKey("previous")) {
    refMode_ = PREVIOUS;
  } else if (!argIn.hasKey("first")) {
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
          delete trj;
          return 1;
        }
        if (trj->AddSingleTrajin( reftraj, argIn, RefParm )) {
          delete trj;
          return 1;
        }
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
    }
  }
  needsSetup_ = true;
  return 0;
}

// ReferenceAction::SetupRefMask()
int ReferenceAction::SetupRefMask(Topology const& topIn) {
  mprintf("\tReference topology: %s\n", topIn.c_str());
  if (refMask_.MaskStringSet()) {
    if (topIn.SetupIntegerMask( refMask_ )) return 1;
    mprintf("\tReference mask:");
    refMask_.BriefMaskInfo();
    mprintf("\n");
  } else {
    refMask_.ResetMask();
    refMask_.SetNatoms( topIn.Natom() );
    refMask_.AddAtomRange(0, topIn.Natom() );
  }
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
    for (int rank = 1; rank < trajComm_.Size(); rank++)
      refFrame_.SendFrame(rank, trajComm_);
  } else
    refFrame_.RecvFrame(0, trajComm_);
  trajComm_.Barrier();
# endif
  selectedRef_.SetCoordinates( refFrame_, refMask_ );
  if (fitRef_)
    refTrans_ = selectedRef_.CenterOnOrigin( useMass_ );
}

// ReferenceAction::SetupRef()
int ReferenceAction::SetupRef(Topology const& topIn, int Ntgt) {
  if (needsSetup_) {
    if (refMode_ == FIRST || refMode_ == PREVIOUS) {
      // First frame or previous frame. Topology will be topIn.
      if ( SetupRefMask( topIn ) ) return 1;
    } else {
      // Specified frame or trajectory. Accessed via refCrd_.
      if ( SetupRefMask( refCrd_->Top() ) ) return 1;
      if (refMode_ == FRAME)
        SelectRefAtoms( ((DataSet_Coords_REF*)refCrd_)->RefFrame() );
    }
    needsSetup_ = false;
  } else {
    if (refMode_ == PREVIOUS)
      mprintf("Warning: 'previous' may not work properly for changing topologies.\n");
  }
  // Check that num atoms in target mask from this parm match ref parm mask
  if ( Ntgt != -1 && refMask_.Nselected() != Ntgt ) {
    mprintf("Warning: Number of atoms in target mask (%i) does not equal\n"
            "Warning:   number of atoms in reference mask (%i).\n",
            Ntgt, refMask_.Nselected());
    return 1;
  }
  return 0;
}
