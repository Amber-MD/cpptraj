#include "Builder.h"
#include "BuildAtom.h"
#include "Chirality.h"
#include "Zmatrix.h"
#include "../CpptrajStdio.h"
#include "../Frame.h"
#include "../Topology.h"
#include <algorithm> // std::copy

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
Builder::Builder() :
  debug_(0)
{}

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
  int natom0 = frag0Top.Natom();
  int newNatom = natom0 + frag1Top.Natom();

  // Reserve space in build atom array (chirality, priority, position status).
  // The positions of atoms in fragment 0 will be "known"
  atoms_.assign( newNatom, BuildAtom() );
  Barray posKnown( newNatom, false );
  for (int at = 0; at != natom0; at++)
    posKnown[at] = true;

  // Determine what bondAt1 will be in the combined topology.
  int bondAt1_new = bondAt1 + natom0;

  // Get the chirality of each atom before the bond is added.
  // This is because when the bond is added the geometry between the
  // bonded atoms will likely not be correct, so while priority can be
  // determined, the actual chirality can not.
  // TODO store priorities as well?
  // TODO only store for fragment 1?
  for (int at = 0; at != natom0; ++at)
    if (frag0Top[at].Nbonds() > 2)
      atoms_[at].SetChirality( DetermineChirality(at, frag0Top, frag0frm, debug_) );
  int at1 = 0;
  for (int at = natom0; at != newNatom; at++, at1++)
    if (frag1Top[at1].Nbonds() > 2)
      atoms_[at].SetChirality( DetermineChirality(at1, frag1Top, frag1frm, debug_) );

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

  // Store priorities around each atom.
  // TODO only fragment 1?
  double tors = 0;
  for (int at = 0; at != natom0; ++at) {
    if (frag0Top[at].Nbonds() > 2) {
      atoms_[at].SetNbonds(frag0Top[at].Nbonds());
      DetermineChirality(tors, atoms_[at].PriorityPtr(), at, frag0Top, frag0frm, debug_);
    }
  }
  for (int at = natom0; at != newNatom; at++) {
    if (frag0Top[at].Nbonds() > 2) {
      atoms_[at].SetNbonds(frag0Top[at].Nbonds());
      DetermineChirality(tors, atoms_[at].PriorityPtr(), at, frag0Top, frag0frm, debug_);
    }
  }

  mprintf("DEBUG: Atoms:\n");
  for (int at = 0; at != frag0Top.Natom(); ++at) {
    mprintf("DEBUG:\t\t%s (%i) %s", frag0Top.AtomMaskName(at).c_str(), (int)posKnown[at],
            chiralStr(atoms_[at].Chirality()));
    if (!atoms_[at].Priority().empty()) {
      mprintf(" priority:");
      for (std::vector<int>::const_iterator it = atoms_[at].Priority().begin();
                                            it != atoms_[at].Priority().end(); ++it)
        mprintf(" %s", frag0Top.AtomMaskName(*it).c_str());
      mprintf("\n");
    }
  }

  // Generate a zmatrix for the smaller fragment
  Zmatrix frag1Zmatrix;
  frag1Zmatrix.SetDebug( 2 ); // FIXME
  if (frag1Zmatrix.SetupICsAroundBond( bondAt0, bondAt1_new, frag0frm, frag0Top, posKnown, atoms_ ))
  {
    mprinterr("Error: Zmatrix setup for ICs around %s - %s failed.\n",
              frag0Top.AtomMaskName(bondAt0).c_str(),
              frag0Top.AtomMaskName(bondAt1).c_str());
    return 1;
  }
  frag1Zmatrix.print( &frag0Top );
  if (frag1Zmatrix.SetToFrame( frag0frm, posKnown )) {
    mprinterr("Error: Conversion from fragment 1 Zmatrix to Cartesian coords failed.\n");
    return 1;
  }

  return 0;
}
