#include "BondSearch.h"
#include "DistRoutines.h"
#include "CpptrajStdio.h"
#ifdef TIMER
# include "Timer.h"
#endif

/** Search for bonds between atoms in residues and atoms in adjacent residues
  * using distance-based criterion that depends on atomic elements.
  * \param top Topology to add bonds to.
  * \param frameIn Frame containing atomic coordinates.
  * \param offset Offset to add when determining if a bond is present.
  * \param debug If > 0 print extra info.
  */
int BondSearch( Topology& top, Frame const& frameIn, double offset, int debug) {
  mprintf("\t%s: determining bond info from distances.\n", top.c_str());
  if (frameIn.empty()) {
    mprinterr("Internal Error: No coordinates set; cannot search for bonds.\n");
    return 1;
  }
# ifdef TIMER
  Timer time_total, time_within, time_between;
  time_total.Start();
  time_within.Start();
# endif
  // ----- STEP 1: Determine bonds within residues
  for (Topology::res_iterator res = top.ResStart(); res != top.ResEnd(); ++res)
  {
    int stopatom = res->LastAtom();
    // Check for bonds between each atom in the residue.
    for (int atom1 = res->FirstAtom(); atom1 != stopatom; ++atom1) {
      Atom::AtomicElementType a1Elt = top[atom1].Element();
      // If this is a hydrogen and it already has a bond, move on.
      if (a1Elt==Atom::HYDROGEN && top[atom1].Nbonds() > 0 )
        continue;
      for (int atom2 = atom1 + 1; atom2 != stopatom; ++atom2) {
        Atom::AtomicElementType a2Elt = top[atom2].Element();
        double D2 = DIST2_NoImage(frameIn.XYZ(atom1), frameIn.XYZ(atom2) );
        double cutoff2 = Atom::GetBondLength(a1Elt, a2Elt) + offset;
        cutoff2 *= cutoff2;
        if (D2 < cutoff2) {
          top.AddBond(atom1, atom2);
          // Once a bond has been made to hydrogen move on.
          if (a1Elt==Atom::HYDROGEN) break;
        }
      }
    }
  }
# ifdef TIMER
  time_within.Stop();
  time_between.Start();
# endif
  // ----- STEP 2: Determine bonds between adjacent residues
  Topology::mol_iterator nextmol = top.MolStart();
  if (top.Nmol() > 0)
    ++nextmol;
  for (Topology::res_iterator res = top.ResStart() + 1; res != top.ResEnd(); ++res)
  {
    // If molecule information is already present, check if first atom of 
    // this residue >= first atom of next molecule, which indicates this
    // residue and the previous residue are in different molecules.
    if ( (nextmol != top.MolEnd()) &&
         (res->FirstAtom() >= nextmol->BeginAtom()) )
    {
      ++nextmol;
      continue;
    }
    // If this residue is recognized as solvent, no need to check previous or
    // next residue
    if ( res->NameIsSolvent() ) {
      ++res;
      if (res == top.ResEnd()) break;
      continue;
    }
    // Get previous residue
    Topology::res_iterator previous_res = res - 1;
    // If previous residue is recognized as solvent, no need to check previous.
    if ( previous_res->NameIsSolvent() ) continue;
    // Get previous residue start atom
    int startatom = previous_res->FirstAtom();
    // Previous residue stop atom, this residue start atom
    int midatom = res->FirstAtom();
    // This residue stop atom
    int stopatom = res->LastAtom();
    // Check for bonds between adjacent residues
    for (int atom1 = startatom; atom1 != midatom; atom1++) {
      Atom::AtomicElementType a1Elt = top[atom1].Element();
      if (a1Elt==Atom::HYDROGEN) continue;
      for (int atom2 = midatom; atom2 != stopatom; atom2++) {
        Atom::AtomicElementType a2Elt = top[atom2].Element();
        if (a2Elt==Atom::HYDROGEN) continue;
        double D2 = DIST2_NoImage(frameIn.XYZ(atom1), frameIn.XYZ(atom2) );
        double cutoff2 = Atom::GetBondLength(a1Elt, a2Elt) + offset;
        cutoff2 *= cutoff2;
        if (D2 < cutoff2)
          top.AddBond(atom1, atom2);
      }
    }
  }
# ifdef TIMER
  time_between.Stop();
  time_total.Stop();
  time_within.WriteTiming(2, "Distances within residues", time_total.Total());
  time_between.WriteTiming(2, "Distances between residues", time_total.Total());
  time_total.WriteTiming(1, "Total for determining bonds via distances");
# endif
  if (debug > 0)
    mprintf("\t%s: %zu bonds to hydrogen, %zu other bonds.\n", top.c_str(),
            top.BondsH().size(), top.Bonds().size());
  return 0;
}
