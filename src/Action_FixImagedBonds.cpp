#include <stack>
#include "Action_FixImagedBonds.h"
#include "CpptrajStdio.h"

// Action_FixImagedBonds::Help()
void Action_FixImagedBonds::Help() const {

}

// Action_FixImagedBonds::Init()
Action::RetType Action_FixImagedBonds::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Always image
  image_.InitImaging( true );

  mask_.SetMaskString( actionArgs.GetMaskNext() );

  mprintf("\tChecking all bonds selected by mask '%s'\n", mask_.MaskString());

  return Action::OK;
}

// Action_FixImagedBonds::Setup()
Action::RetType Action_FixImagedBonds::Setup(ActionSetup& setup)
{
  // Need box info
  if (setup.CoordInfo().TrajBox().Type()==Box::NOBOX) {
    mprintf("Warning: Topology '%s' does not contain box information; required for imaging.\n",
            setup.Top().c_str());
    return Action::SKIP;
  }
  // Setup mask
  if ( setup.Top().SetupCharMask( mask_ ) ) return Action::ERR;
  mask_.MaskInfo();
  if ( mask_.None() ) {
    mprintf("Warning: No atoms selected.\n");
    return Action::SKIP;
  }
  // Figure out the first and last selected atoms.
  firstSelected_ = -1;
  lastSelected_ = -1;
  for (int at = 0; at != setup.Top().Natom(); at++) {
    if (mask_.AtomInCharMask( at )) {
      if (firstSelected_ == -1)
        firstSelected_ = at;
      lastSelected_ = at;
    }
  }
  mprintf("\tFirst selected atom %i, last selected atom %i\n", firstSelected_+1, lastSelected_+1);
  // Set up imaging info for this parm
  image_.SetupImaging( setup.CoordInfo().TrajBox().Type() );

  CurrentParm_ = setup.TopAddress();

  atomVisited_.assign( lastSelected_+1, false );

  return Action::OK;
}

// Action_FixImagedBonds::DoAction()
Action::RetType Action_FixImagedBonds::DoAction(int frameNum, ActionFrame& frm)
{
  Box const& box = frm.Frm().BoxCrd();
  // Calculate box info needed for imaging based on cell type
  if (image_.ImageType() == NONORTHO)
    box.ToRecip(ucell_, recip_);
  else
    boxCenter_ = box.Center();
  // Starting with the first atom, check every atom bonded to that atom
  // pseudo-recursively. Ensure that no bond length is longer than half
  // the box size. If it is adjust the position of the bonded atom to
  // correct this. This has to be done 1 atom at a time because if an atom
  // is moved then the atoms it is bonded to must move as well.
  std::stack<unsigned int> nextAtomToSearch;
  Topology const& Top = *CurrentParm_;
  unsigned int Natoms = (unsigned int)lastSelected_ + 1;
  // Set any unselected atoms to true so they will be skipped, all others
  // to false.
  for (unsigned int i = 0; i != Natoms; i++)
    if ( mask_.AtomInCharMask(i) )
      atomVisited_[i] = false;
    else
      atomVisited_[i] = true;
  
  bool uncheckedAtomsRemain = true;
  unsigned int currentAtom = (unsigned int)firstSelected_;
  unsigned int lowestUnassignedAtom = currentAtom;
  // BEGIN main loop
  while (uncheckedAtomsRemain) {
    atomVisited_[currentAtom] = true;
    Atom const& AT = Top[currentAtom];
    Vec3 currXYZ( frm.Frm().XYZ( currentAtom ) );
    mprintf("ANCHOR: %i\n", currentAtom);
    // All atoms bonded to this one are in molecule.
    for (Atom::bond_iterator batom = AT.bondbegin(); batom != AT.bondend(); ++batom)
    {
      if (!atomVisited_[*batom]) {
        if ( Top[*batom].Nbonds() > 1 )
          // Bonded atom has more than 1 bond; needs to be searched.
          nextAtomToSearch.push( *batom );
        else
          // Bonded atom only bonded to current atom. No more search needed.
          atomVisited_[*batom] = true;
        // Check if bonded atom is too far away
        Vec3 bondXYZ( frm.Frm().XYZ( *batom ) );
        Vec3 delta = bondXYZ - currXYZ;
        Vec3 boxTrans(0.0);
        if ( image_.ImageType() == ORTHO ) {
          // ----- Orthorhombic imaging ------------
          // If the distance between current and bonded atom is more than half the box,
          // adjust the position of the bonded atom.
          while (delta[0] >  boxCenter_[0]) { delta[0] -= box.BoxX(); boxTrans[0] -= box.BoxX(); }
          while (delta[0] < -boxCenter_[0]) { delta[0] += box.BoxX(); boxTrans[0] += box.BoxX(); }
          while (delta[1] >  boxCenter_[1]) { delta[1] -= box.BoxY(); boxTrans[1] -= box.BoxY(); }
          while (delta[1] < -boxCenter_[1]) { delta[1] += box.BoxY(); boxTrans[1] += box.BoxY(); }
          while (delta[2] >  boxCenter_[2]) { delta[2] -= box.BoxZ(); boxTrans[2] -= box.BoxZ(); }
          while (delta[2] < -boxCenter_[2]) { delta[2] += box.BoxZ(); boxTrans[2] += box.BoxZ(); }
        } else {
          // ----- Non-orthorhombic imaging --------
          Vec3 fdelta = recip_ * delta;
          // DEBUG
//          Vec3 dbgdelta = (recip_ * bondXYZ) - (recip_ * currXYZ);
//          fdelta.Print("fdelta");
//          dbgdelta.Print("dbgdelta");
          // If the distance between current and bonded atom is more than half the cell,
          // adjust position of the bonded atom.
          while (fdelta[0] >  0.5) { fdelta[0] -= 1.0; boxTrans[0] -= 1.0; }
          while (fdelta[0] < -0.5) { fdelta[0] += 1.0; boxTrans[0] += 1.0; }
          while (fdelta[1] >  0.5) { fdelta[1] -= 1.0; boxTrans[1] -= 1.0; }
          while (fdelta[1] < -0.5) { fdelta[1] += 1.0; boxTrans[1] += 1.0; }
          while (fdelta[2] >  0.5) { fdelta[2] -= 1.0; boxTrans[2] -= 1.0; }
          while (fdelta[2] < -0.5) { fdelta[2] += 1.0; boxTrans[2] += 1.0; }
          boxTrans = ucell_.TransposeMult( boxTrans );
          delta = ucell_.TransposeMult( fdelta ); // DEBUG
        }
        // Translate the atom
        if (boxTrans[0] != 0.0 || boxTrans[1] != 0.0 || boxTrans[2] != 0.0)
          mprintf("\tMove atom %6i d={%8.3g %8.3g %8.3g} t={%8.3g %8.3g %8.3g}\n",
                  *batom, delta[0], delta[1], delta[2], boxTrans[0], boxTrans[1], boxTrans[2]);
        frm.ModifyFrm().Translate(boxTrans, *batom);
      } // END if not visited
    } // END loop over bonded atoms
    if (nextAtomToSearch.empty()) {
      // No more atoms to search. Find next unmarked atom.
      unsigned int idx = lowestUnassignedAtom;
      for (; idx != Natoms; idx++)
        if (!atomVisited_[idx]) break;
      if (idx == Natoms)
        uncheckedAtomsRemain = false;
      else {
        currentAtom = idx;
        lowestUnassignedAtom = idx + 1;
      }
    } else {
      currentAtom = nextAtomToSearch.top();
      nextAtomToSearch.pop();
    }
  } // END main loop
  return Action::MODIFY_COORDS;
}
