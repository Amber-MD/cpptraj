#include "Model.h"
#include "../CpptrajStdio.h"
#include "../Topology.h"
#include "../Chirality.h"
#include "../TorsionRoutines.h"
#include "../Constants.h"

/** Attempt to assign a reasonable value for phi internal coordinate for atom i
  * given that atoms j k and l have known positions.
  *   j - k
  *  /     \
  * i       l
  */
int Cpptraj::Structure::Model::AssignPhi(double& phi, int ai, int aj, int ak, int al, Topology const& topIn, Frame const& frameIn, std::vector<bool> const& atomPositionKnown)
{
  using namespace Cpptraj::Chirality;
  // Figure out hybridization and chirality of atom j.
  mprintf("DEBUG: AssignPhi for atom j : %s\n", topIn.AtomMaskName(aj).c_str());

  Atom const& AJ = topIn[aj];
  mprintf("DEBUG:\t\tNbonds: %i\n", AJ.Nbonds());

  // FIXME aj ak and al should be known
  // TODO check that atom i actually ends up on the list?
  std::vector<int> Priority( AJ.Nbonds() + 1 );
  int* priority = static_cast<int*>( &Priority[0] );
  double tors = 0;
  // DetermineChirality will always put aj as the last index.
  ChiralType chirality = DetermineChirality(tors, priority, aj, topIn, frameIn, 0); // FIXME debug
  if (chirality == ERR) {
    mprinterr("Error: Problem determining chirality around atom %s\n", topIn.AtomMaskName(aj).c_str());
    return 1;
  } else if (chirality == IS_UNKNOWN_CHIRALITY) {
    mprintf("DEBUG:\t\tChirality is unknown\n");
  } else if (chirality == IS_S) {
    mprintf("DEBUG:\t\tChirality is S\n");
  } else {
    mprintf("DEBUG:\t\tChirality is R\n");
  }
  mprintf("DEBUG:\t\tPriority around J:");
  for (int idx = 0; idx < AJ.Nbonds(); idx++)
    mprintf(" %s(%i)", topIn.AtomMaskName(priority[idx]).c_str(), (int)atomPositionKnown[priority[idx]]);
  mprintf("\n");

  // Fill in what values we can for known atoms
  std::vector<double> knownPhi( AJ.Nbonds() );
  for (int idx = 0; idx < AJ.Nbonds(); idx++) {
    int atnum = priority[idx];
    if (atnum != ak && atomPositionKnown[atnum]) {
      knownPhi[idx] = Torsion(frameIn.XYZ(atnum),
                              frameIn.XYZ(aj),
                              frameIn.XYZ(ak),
                              frameIn.XYZ(al));
      mprintf("DEBUG:\t\tKnown phi for %s = %g\n", topIn.AtomMaskName(atnum).c_str(), knownPhi[idx]*Constants::RADDEG);
    }
  }

  if (AJ.Nbonds() > 1) {
    // The interval will be 360 / (number of bonds - 1)
    double interval = Constants::TWOPI / (AJ.Nbonds() - 1);
    mprintf("DEBUG:\t\tInterval is %g degrees\n", interval * Constants::RADDEG);
  }

  return 0;
}