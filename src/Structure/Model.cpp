#include "Model.h"
#include "../CpptrajStdio.h"
#include "../Topology.h"
#include "../Chirality.h"

using namespace Cpptraj::Structure::Model;

/** Attempt to assign a reasonable value for phi internal coordinate for atom i
  * given that atoms j k and l have known positions.
  *   j - k
  *  /     \
  * i       l
  */
int AssignPhi(double& phi, int ai, int aj, int ak, int al, Topology const& topIn, Frame const& frameIn)
{
  using namespace Cpptraj::Chirality;
  // Figure out hybridization and chirality of atom j.
  mprintf("DEBUG: AssignPhi for atom j : %s\n", topIn.AtomMaskName(aj).c_str());

  Atom const& AJ = topIn[aj];
  mprintf("DEBUG:\t\tNbonds: %i\n", AJ.Nbonds());

  int priority[4];
  double tors = 0;
  ChiralType chirality = DetermineChirality(tors, priority, aj, topIn, frameIn, 0); // FIXME debug
  if (chirality == ERR) {
    mprinterr("Error: Could not determine chirality around atom %s\n", topIn.AtomMaskName(aj).c_str());
    return 1;
  } else if (chirality == IS_S) {
    mprintf("DEBUG:\t\tChirality is S\n");
  } else {
    mprintf("DEBUG:\t\tChirality is R\n");
  }
  mprintf("DEBUG:\t\tPriority around J:");
  for (int idx = 0; idx < 3; idx++)
    mprintf(" %s", topIn.AtomMaskName(priority[idx]).c_str());
  mprintf("\n");

  return 0;
}
