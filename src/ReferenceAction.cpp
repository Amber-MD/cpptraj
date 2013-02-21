#include "ReferenceAction.h"
#include "CpptrajStdio.h"

// ReferenceAction::SetRefMask()
int ReferenceAction::SetRefMask(Topology const& topIn, const char* call) {
  if (topIn.SetupIntegerMask( RefMask_ )) return 1;
  if (RefMask_.None()) {
    mprinterr("Error: %s: No reference atoms selected for parm %s, [%s]\n",
              call, topIn.c_str(), RefMask_.MaskString());
    return 1;
  }
  SelectedRef_.SetupFrameFromMask( RefMask_, topIn.Atoms() );
  return 0;
}

// ReferenceAction::InitRef()
int ReferenceAction::InitRef(bool previousIn, bool firstIn, bool massIn, bool nofitIn,
                             std::string const& reftrajname, ReferenceFrame& REF, 
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
  RefMask_.SetMaskString(refmaskIn);
  // Initialize reference if not 'first'
  if (refmode_ != FIRST) {
    if ( !reftrajname.empty() ) {
      // Reference trajectory
      if (SetRefMask( *RefParm, call )!=0) return 1;
      // Attempt to open reference traj.
      if (RefTraj_.SetupTrajRead( reftrajname, &actionArgs, RefParm)) {
        mprinterr("Error: %s: Could not set up reftraj %s\n", call, reftrajname.c_str());
        return 1;
      }
      RefFrame_.SetupFrameV(RefParm->Atoms(), RefTraj_.HasVelocity());
      if (RefTraj_.BeginTraj(false)) {
        mprinterr("Error: %s: Could not open reftraj %s\n", call, reftrajname.c_str());
        return 1;
      }
    } else {
      // Reference Frame
      if (SetRefMask( *(REF.Parm()), call ) != 0) return 1;
      SetRefStructure( *(REF.Coord()), nofitIn, massIn );
    }
  }

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
  if ( RefMask_.Nselected() != Ntgt ) {
    mprintf( "Warning: Number of atoms in target mask (%i) does not \n", Ntgt);
    mprintf( "         equal number of atoms in Ref mask (%i).\n",RefMask_.Nselected());
    return 1;
  }
  return 0;
}
