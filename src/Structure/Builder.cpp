#include "Builder.h"
#include "../CpptrajStdio.h"
#include "../Topology.h"
#include "../Frame.h"
#include <algorithm> // std::copy

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
Builder::Builder() :
  debug_(0)
{}

/// Wrapper for DetermineChirality that checks number of bonds TODO bake into DetermineChirality?
static inline
Cpptraj::Chirality::ChiralType determineChiral(int at, Topology const& topIn, Frame const& frmIn, int dbg)
{
  if (topIn[at].Nbonds() > 2)
    return Cpptraj::Chirality::DetermineChirality(at, topIn, frmIn, dbg);
  else
    return Cpptraj::Chirality::IS_UNKNOWN_CHIRALITY;
}

/** \return Heavy atom count */
int Builder::heavy_atom_count(Topology const& topIn) {
  int hac = 0;
  for (int i = 0; i != topIn.Natom(); i++) {
    if (topIn[i].Element() != Atom::HYDROGEN &&
        topIn[i].Element() != Atom::EXTRAPT)
      hac++;
  }
  return hac;
}

/** Combine two units. The smaller fragment will be merged into 
  * the larger fragment.
  */
int Builder::Combine(Topology& frag0Top, Frame& frag0frm,
                     Topology& frag1Top, Frame& frag1frm,
                     int bondAt0, int bondAt1)
{
  if (heavy_atom_count(frag0Top) < heavy_atom_count(frag1Top)) {
    return combine01(frag1Top, frag1frm, frag0Top, frag0frm, bondAt1, bondAt0);
  } else {
    return combine01(frag0Top, frag0frm, frag1Top, frag1frm, bondAt0, bondAt1);
  }
}

/** Combine two units. Fragment 1 will be merged into Fragment 0. */
int Builder::combine01(Topology& frag0Top, Frame& frag0frm,
                     Topology const& frag1Top, Frame const& frag1frm,
                     int bondAt0, int bondAt1)
{
  using namespace Cpptraj::Chirality;
  int newNatom = frag0Top.Natom() + frag1Top.Natom();

  // The positions of atoms in fragment 0 will be "known"
  Barray posKnown( newNatom, false );
  for (int at = 0; at != frag0Top.Natom(); at++)
    posKnown[at] = true;

  // Determine what bondAt1 will be in the combined topology.
  int bondAt1_new = bondAt1 + frag0Top.Natom();

  // Reserve space in atom chirality and priority arrays.
  atomChirality_.assign( newNatom, IS_UNKNOWN_CHIRALITY );
  atomPriority_.resize( newNatom );

  // Get the chirality of each atom before the bond is added.
  // This is because when the bond is added the geometry will likely
  // not be correct, so while priority can be determined, the
  // actual chirality can not.
  // TODO store priorities as well?
  // TODO only store for fragment 1?
  for (int at = 0; at != frag0Top.Natom(); ++at)
    atomChirality_[at] = determineChiral(at, frag0Top, frag0frm, debug_);
  int at1 = 0;
  for (int at = frag0Top.Natom(); at != newNatom; at++, at1++)
    atomChirality_[at1] = determineChiral(at, frag1Top, frag1frm, debug_);

  // Merge fragment 1 topology into fragment 0
  frag0Top.AppendTop( frag1Top );
  // Merge fragment 1 coords into fragment 0
  double* tmpcrd = new double[ frag0frm.size() ];
  unsigned int tmpsize = frag0frm.size();
  std::copy( frag0frm.xAddress(), frag0frm.xAddress()+frag0frm.size(), tmpcrd );
  frag0frm.SetupFrameV( frag0Top.Atoms(), CoordinateInfo(frag0frm.BoxCrd(), false, false, false) );
  std::copy( tmpcrd,              tmpcrd+tmpsize,                      frag0frm.xAddress() );
  std::copy( frag1frm.xAddress(), frag1frm.xAddress()+frag1frm.size(), frag0frm.xAddress()+tmpsize );

  // Create the bond
  mprintf("DEBUG: Bond %i %s to %i %s\n",
          bondAt0+1,     frag0Top.AtomMaskName(bondAt0).c_str(),
          bondAt1_new+1, frag0Top.AtomMaskName(bondAt1_new).c_str());
  frag0Top.AddBond( bondAt0, bondAt1_new ); // TODO create pseudo-parameter?
  // Regenerate the molecule info TODO should AddBond do this automatically?
  if (frag0Top.DetermineMolecules()) {
    mprinterr("Error: Could not determine molecule info after combining fragments.\n");
    return 1;
  }

  return 0;
}
