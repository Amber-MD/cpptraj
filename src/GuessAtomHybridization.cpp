#include "GuessAtomHybridization.h"
#include "Atom.h"

AtomType::HybridizationType Cpptraj::GuessAtomHybridization(Atom const& AJ, std::vector<Atom> const& atoms)
{
  AtomType::HybridizationType hybrid = AtomType::UNKNOWN_HYBRIDIZATION;
  // Handle specific elements
  switch (AJ.Element()) {
    case Atom::CARBON :
      switch (AJ.Nbonds()) {
        case 2 : hybrid = AtomType::SP; break;
        case 3 : hybrid = AtomType::SP2; break;
        case 4 : hybrid = AtomType::SP3; break;
      }
      break;
    case Atom::NITROGEN :
      switch (AJ.Nbonds()) {
        case 2 : hybrid = AtomType::SP2; break;
        case 3 :
          // Check for potential SP2. If only 1 of the bonded atoms is
          // hydrogen, assume SP2. TODO actually check for aromaticity.
          int n_hydrogens = 0;
          for (Atom::bond_iterator bat = AJ.bondbegin(); bat != AJ.bondend(); ++bat)
            if (atoms[*bat].Element() == Atom::HYDROGEN)
              n_hydrogens++;
          if (n_hydrogens == 1)
            hybrid = AtomType::SP2;
          else
            hybrid = AtomType::SP3;
          break;
      }
      break;
    case Atom::OXYGEN :
      switch (AJ.Nbonds()) {
        case 2 : hybrid = AtomType::SP3; break;
      }
      break;
    case Atom::SULFUR :
      switch (AJ.Nbonds()) {
        case 2 : hybrid = AtomType::SP3; break;
      }
      break;
    default: hybrid = AtomType::UNKNOWN_HYBRIDIZATION; break;
  }
/*
if (hybrid == AtomType::UNKNOWN_HYBRIDIZATION) {
    // Assign a theta based on number of bonds 
    switch (AJ.Nbonds()) {
      case 4 : hybrid = AtomType::SP3; break;
      case 3 : hybrid = AtomType::SP2; break;
      case 2 : hybrid = AtomType::SP; break;
      default : mprinterr("Internal Error: GuessAtomHybridization(): Unhandled # bonds for %s (%i)\n", topIn.AtomMaskName(aj).c_str(), AJ.Nbonds()); return 1;
    }
*/

  return hybrid;
}
