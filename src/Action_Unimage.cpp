#include "Action_Unimage.h"
#include "CpptrajStdio.h"
#include "ImageRoutines.h" // CreatePairList

// Action_Unimage::Help()
void Action_Unimage::Help() const {

}

// Action_Unimage::Init()
Action::RetType Action_Unimage::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Always image
  image_.InitImaging( true );
  useCenter_ = actionArgs.hasKey("center");
  if (actionArgs.hasKey("bymol"))
    imageMode_ = Image::BYMOL;
  else if (actionArgs.hasKey("byres"))
    imageMode_ = Image::BYRES;
  else if (actionArgs.hasKey("byatom")) {
    imageMode_ = Image::BYATOM;
    // Unwrapping to center by atom makes no sense
    if (useCenter_) useCenter_ = false;
  } else
    imageMode_ = Image::BYATOM;
  // Get reference
  ReferenceFrame REF = init.DSL().GetReferenceFrame( actionArgs );
  if (REF.error()) return Action::ERR;
  // Get mask
  maskExp_ = actionArgs.GetMaskNext();
  // Set up previous coords
  if (!REF.empty()) {
    AtomMask refMask( maskExp_ );
    if (REF.Parm().SetupIntegerMask( refMask, REF.Coord() )) return Action::ERR;
    refMask.MaskInfo();
    if (refMask.None()) {
      mprinterr("Error: No atoms selected by '%s' for '%s'\n",
                refMask.MaskString(), REF.refName());
      return Action::ERR;
    }
    previous_.SetupFrame( refMask.Nselected() );
    unsigned int idx = 0;
    for (AtomMask::const_iterator at = refMask.begin(); at != refMask.end(); ++at, idx += 3)
    {
      const double* XYZ = REF.Coord().XYZ( *at );
      previous_[idx  ] = XYZ[0];
      previous_[idx+1] = XYZ[1];
      previous_[idx+2] = XYZ[2];
    }
  }

  mprintf("    UNIMAGE: By %s", Image::ModeString(imageMode_));
  if (maskExp_.empty())
    mprintf(" using all atoms");
  else
    mprintf(" using mask '%s'", maskExp_.c_str());
  if (imageMode_ != Image::BYATOM) {
    if (useCenter_)
      mprintf(" based on center of mass.");
    else
      mprintf(" based on first atom position.");
  }
  mprintf("\n");
  if ( !REF.empty())
    mprintf("\tReference is %s", REF.refName());
  else
    mprintf("\tReference is first frame.");
  mprintf("\n");
  return Action::OK;
}

// Action_Unimage::Setup()
Action::RetType Action_Unimage::Setup(ActionSetup& setup)
{
  // Need box info
  if (setup.CoordInfo().TrajBox().Type()==Box::NOBOX) {
    mprintf("Warning: Topology '%s' does not contain box information; required for imaging.\n",
            setup.Top().c_str());
    return Action::SKIP;
  }
  // Set up atom pairs to be unwrapped 
  atomPairs_ = Image::CreatePairList(setup.Top(), imageMode_, maskExp_);
  if (atomPairs_.empty()) {
    mprintf("Warning: No atoms selected by '%s'\n", maskExp_.c_str());
    return Action::SKIP;
  }
  natoms_ = atomPairs_.size() / 2;
  mprintf("\tNumber of %ss to be unimaged is %i\n",
          Image::ModeString(imageMode_), natoms_);
  // If previous frame is set up ensure same number of atoms
  if (!previous_.empty()) {
    if (previous_.Natom() != natoms_) {
      mprinterr("Error: # reference atoms %i != # selected atoms.\n", natoms_);
      return Action::ERR;
    }
  }
  return Action::OK;
}

// Action_Unimage::DoAction()
Action::RetType Action_Unimage::DoAction(int frameNum, ActionFrame& frm)
{
  if (previous_.empty()) {
    previous_.SetupFrame( natoms_ );
    // Load this frame as reference
    unsigned int idx = 0;
    for (Iarray::const_iterator ap = atomPairs_.begin(); ap != atomPairs_.end(); ap += 2)
    {
      for (int at = atomPairs_[*ap]; at != atomPairs_[*(ap+1)]; at++, idx += 3)
      {
        const double* XYZ = frm.Frm().XYZ( at );
        previous_[idx  ] = XYZ[0];
        previous_[idx+1] = XYZ[1];
        previous_[idx+2] = XYZ[2];
      }
    }
    return Action::OK;
  } else {
    Vec3 currXYZ, prevXYZ;
    boxCenter_ = frm.Frm().BoxCrd().Center();
    // Calculate matrices if necessary
    frm.Frm().BoxCrd().ToRecip(ucell_, recip_);
    // Perform unwrapping
    unsigned int pidx = 0;
    for (Iarray::const_iterator ap = atomPairs_.begin(); ap != atomPairs_.end(); ap += 2)
    {
      int eltsize = *(ap+1) - *ap;
      if (useCenter_) {
        // Use c.o.m. of elements
        currXYZ = frm.Frm().VGeometricCenter(*ap, *(ap+1));
        prevXYZ = previous_.VGeometricCenter(pidx, pidx + eltsize);
      } else {
        // Use first atom position only
        currXYZ = frm.Frm().XYZ( *ap );
        prevXYZ = previous_.XYZ( pidx );
      }
      // Calculate distance to previous frames elements coordinates
      double delx = currXYZ[0] - prevXYZ[0];
      double dely = currXYZ[1] - prevXYZ[1];
      double delz = currXYZ[2] - prevXYZ[2];
      if ( image_.ImageType() == ORTHO ) {
        // ----- Orthorhombic imaging ------------
        // If the element moved more than half the box, assume it was imaged
        // and adjust the current position.
        if      (delx >  boxCenter_[0]) boxTrans_[0] = -frm.Frm().BoxCrd().BoxX();
        else if (delx < -boxCenter_[0]) boxTrans_[0] =  frm.Frm().BoxCrd().BoxX();
        if      (dely >  boxCenter_[1]) boxTrans_[1] = -frm.Frm().BoxCrd().BoxY();
        else if (dely < -boxCenter_[1]) boxTrans_[1] =  frm.Frm().BoxCrd().BoxY();
        if      (delz >  boxCenter_[2]) boxTrans_[2] = -frm.Frm().BoxCrd().BoxZ();
        else if (delz < -boxCenter_[2]) boxTrans_[2] =  frm.Frm().BoxCrd().BoxZ();
      } else {
        // ----- Non-orthorhombic imaging --------
        // If the element moved more than half the box, assume it was imaged.
        if (delx > boxCenter_[0] || delx < -boxCenter_[0] ||
            dely > boxCenter_[1] || dely < -boxCenter_[1] ||
            delz > boxCenter_[2] || delz < -boxCenter_[2])
        {
          // Previous position in Cartesian space is in prevXYZ
          // Current position in fractional coords
          Vec3 cFrac = recip_ * currXYZ;
          // Look for imaged distance closer to previous than current position
          double minDist2 = frm.Frm().BoxCrd().BoxX() *
                            frm.Frm().BoxCrd().BoxY() *
                            frm.Frm().BoxCrd().BoxZ();
          Vec3 minCurr(0.0);
          for (int ix = -1; ix < 2; ix++) {
            for (int iy = -1; iy < 2; iy++) {
              for (int iz = -1; iz < 2; iz++) {
                if (ix != 0 || iy != 0 || iz != 0) { // Ignore current position
                  Vec3 ixyz(ix, iy, iz);
                  // Current position shifted and back in Cartesian space
                  Vec3 IMG = ucell_.TransposeMult(cFrac + ixyz);
                  // Distance from previous position to imaged current position
                  Vec3 dxyz = IMG - prevXYZ;
                  double dist2 = dxyz.Magnitude2();
                  if (dist2 < minDist2) {
                    minDist2 = dist2;
                    minCurr = IMG;
                  }
                }
              }
            }
          }
          boxTrans_[0] = minCurr[0] - currXYZ[0];
          boxTrans_[1] = minCurr[1] - currXYZ[1];
          boxTrans_[2] = minCurr[2] - currXYZ[2];
        }
      }
      // Move each atom of the element
      frm.ModifyFrm().Translate(boxTrans_, *ap, *(ap+1));
      // Save the new positions
      int i3 = *ap * 3;
      std::copy( frm.Frm().xAddress()+i3, frm.Frm().xAddress()+(*(ap+1)*3),
                 previous_.xAddress()+pidx*3 );

      pidx += eltsize;
    } // END loop over atom pairs
  }
 
  return Action::MODIFY_COORDS;
}
