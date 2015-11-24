#include "ReferenceAction.h"
#include "CpptrajStdio.h"
#ifdef MPI
# include "Parallel.h"
#endif

// ReferenceAction::SetRefMask()
int ReferenceAction::SetRefMask(Topology const& topIn, const char* call) {
  if (topIn.SetupIntegerMask( refMask_ )) return 1;
  mprintf("\tReference mask:");
  refMask_.BriefMaskInfo();
  mprintf("\n");
  if (refMask_.None()) {
    mprinterr("Error: %s: No reference atoms selected for parm %s, [%s]\n",
              call, topIn.c_str(), refMask_.MaskString());
    return 1;
  }
  selectedRef_.SetupFrameFromMask( refMask_, topIn.Atoms() );
  return 0;
}

// ReferenceAction::InitRef()
int ReferenceAction::InitRef(bool previousIn, bool firstIn, bool massIn, bool fitIn,
                             std::string const& reftrajname, ReferenceFrame const& REF, 
                             Topology* RefParm, std::string const& refmaskIn, 
                             ArgList& actionArgs, const char* call)
{
  refmode_ = UNKNOWN_REF;
  previous_ = previousIn;
  if (firstIn || previous_)
    refmode_ = FIRST;
  else {
    if (REF.error()) return 1; 
    if (REF.empty()) {
      if (!reftrajname.empty()) {
        if (RefParm == 0) {
          mprinterr("Error: %s: No parm found for reftraj %s. Make sure parm has been loaded.\n",
                    call, reftrajname.c_str());
          return 1;
        }
        refmode_ = REFTRAJ;
      } else {
        // No reference keywords specified. Default to first.
        mprintf("Warning: %s: No reference structure given. Defaulting to first.\n",call);
        refmode_ = FIRST;
      }
    } else
      refmode_ = REFFRAME;
  }
  // Set the reference mask expression
  refMask_.SetMaskString(refmaskIn);
  // Initialize reference if not 'first'
  if (refmode_ != FIRST) {
    if ( !reftrajname.empty() ) {
      // Reference trajectory
      if (SetRefMask( *RefParm, call )!=0) return 1;
      // Attempt to open reference traj.
      if (refTraj_.SetupTrajRead( reftrajname, actionArgs, RefParm)) {
        mprinterr("Error: %s: Could not set up reftraj %s\n", call, reftrajname.c_str());
        return 1;
      }
      refFrame_.SetupFrameV(RefParm->Atoms(), refTraj_.TrajCoordInfo());
      if (refTraj_.BeginTraj()) {
        mprinterr("Error: %s: Could not open reftraj %s\n", call, reftrajname.c_str());
        return 1;
      }
    } else {
      // Reference Frame
      if (SetRefMask( REF.Parm(), call ) != 0) return 1;
      SetRefStructure( REF.Coord(), fitIn, massIn );
    }
  }
  // Set reference mode string
  if (previous_)
    modeString_ = "previous frame";
  else if (refmode_ == FIRST)
    modeString_ = "first frame";
  else if (refmode_==REFTRAJ)
    modeString_ = "trajectory " + refTraj_.Traj().Filename().Full();
  else // REFFRAME
    modeString_ = "\"" + REF.RefName() + "\"";
  modeString_ += " (" + refMask_.MaskExpression() + ")";

  return 0;
}

// ReferenceAction::SetupRef()
int ReferenceAction::SetupRef(Topology const& topIn, int Ntgt, const char* call) {
  if (refmode_ == FIRST) {
    if ( SetRefMask( topIn, call )!=0 ) return 1;
  } else if (previous_) {
    mprintf("Warning: %s: 'previous' may not work properly for changing topologies.\n",call);
    if ( SetRefMask( topIn, call )!=0 ) return 1;
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

// ReferenceAction::SetRefStructure()
void ReferenceAction::SetRefStructure(Frame const& frameIn, bool fitIn, bool useMassIn)
{
  refFrame_ = frameIn;
# ifdef MPI
  // Ensure all threads are using the same reference
  if (Parallel::World().Master()) { // TODO MasterBcast
    rprintf("DEBUG: Sending reference frame to children.\n");
    for (int rank = 1; rank != Parallel::World().Size(); rank++)
      refFrame_.SendFrame(rank);
  } else {
    rprintf("DEBUG: Receiving reference frame from master.\n");
    refFrame_.RecvFrame(0);
  }
  Parallel::World().Barrier();
# endif
  selectedRef_.SetCoordinates( refFrame_, refMask_ );
  if (fitIn)
    refTrans_ = selectedRef_.CenterOnOrigin( useMassIn );
}
