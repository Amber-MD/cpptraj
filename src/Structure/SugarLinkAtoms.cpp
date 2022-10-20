#include "SugarLinkAtoms.h"
#include "../CpptrajStdio.h"
#include "../Topology.h"
#include "../StringRoutines.h"

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
SugarLinkAtoms::SugarLinkAtoms(int debugIn) : debug_(debugIn) {}

/// \return Glycam linkage code for given linked atoms
std::string SugarLinkAtoms::GlycamLinkageCode(Topology const& topIn)
const
{
  std::string linkcode;

  // Try to create a link string based on link atom element and position.
  // Check for any unknown positions.
  std::string linkstr;
  for (std::set<LinkAtom>::const_iterator it = linkages_.begin(); it != linkages_.end(); ++it) {
    if (it->Position() < 1) {
      mprinterr("Error: Linkage for atom '%s' has undetermined position in sugar.\n",
                topIn.AtomMaskName(it->Idx()).c_str());
      return linkcode;
    }
    // Carbon is terminal
    if (topIn[it->Idx()].Element() == Atom::CARBON)
      linkstr.append("T");
    else
      linkstr.append( std::string(topIn[it->Idx()].ElementName()) +
                      integerToString(it->Position()) );
  }

  if (debug_ > 0)
    mprintf("DEBUG:\t  linkstr= '%s'\n", linkstr.c_str());
  if      (linkstr == "T") linkcode = "0";
  else if (linkstr == "O1") linkcode = "1";
  else if (linkstr == "TO2") linkcode = "2";
  else if (linkstr == "O2")  linkcode = "2"; // Furanose C2-O2-X
  else if (linkstr == "TO3") linkcode = "3";
  else if (linkstr == "TO4") linkcode = "4";
  else if (linkstr == "TO5") linkcode = "5";
  else if (linkstr == "TO6") linkcode = "6";
  else if (linkstr == "TO2O3") linkcode = "Z";
  else if (linkstr == "TO2O4") linkcode = "Y";
  else if (linkstr == "TO2O6") linkcode = "X";
  else if (linkstr == "TO3O4") linkcode = "W";
  else if (linkstr == "TO3O6") linkcode = "V";
  else if (linkstr == "TO4O6") linkcode = "U";
  else if (linkstr == "TO2O3O4") linkcode = "T";
  else if (linkstr == "TO2O3O6") linkcode = "S";
  else if (linkstr == "TO2O4O6") linkcode = "R";
  else if (linkstr == "TO3O4O6") linkcode = "Q";
  else if (linkstr == "TO2O3O4O6") linkcode = "P";
  if (linkcode.empty())
    mprintf("Warning: Could not determine link code for link atoms '%s'.\n", linkstr.c_str());
  return linkcode;
}
