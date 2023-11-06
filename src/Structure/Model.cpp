#include "Model.h"
#include "BuildAtom.h"
#include "../CpptrajStdio.h"
#include "../Topology.h"
#include "../TorsionRoutines.h"
#include "../Constants.h"

/** Attempt to assign a reasonable value for theta internal coordinate for
  * atom i given that atoms j and k have known positions.
  */
int Cpptraj::Structure::Model::AssignTheta(double& theta, int ai, int aj, int ak, Topology const& topIn, Frame const& frameIn, std::vector<bool> const& atomPositionKnown, int debug)
{
  // Figure out hybridization and chirality of atom j.
  if (debug > 0)
    mprintf("DEBUG: AssignTheta for atom j : %s\n", topIn.AtomMaskName(aj).c_str());

  enum HybridizationType { SP = 0, SP2, SP3, UNKNOWN_HYBRIDIZATION };

  Atom const& AJ = topIn[aj];
  if (debug > 0) {
    mprintf("DEBUG:\t\tI %s Nbonds: %i\n", topIn[ai].ElementName(), topIn[ai].Nbonds());
    mprintf("DEBUG:\t\tJ %s Nbonds: %i\n", AJ.ElementName(), AJ.Nbonds());
    mprintf("DEBUG:\t\tK %s Nbonds: %i\n", topIn[ak].ElementName(), topIn[ak].Nbonds());
  }
  // Sanity check
  if (AJ.Nbonds() < 2) {
    mprinterr("Internal Error: AssignTheta() called for atom J %s with fewer than 2 bonds.\n", topIn.AtomMaskName(aj).c_str());
    return 1;
  }

  HybridizationType hybrid = UNKNOWN_HYBRIDIZATION;
  // Handle specific elements
  switch (AJ.Element()) {
    case Atom::CARBON :
      switch (AJ.Nbonds()) {
        case 2 : hybrid = SP; break;
        case 3 : hybrid = SP2; break;
        case 4 : hybrid = SP3; break;
      }
      break;
    case Atom::NITROGEN :
      switch (AJ.Nbonds()) {
        case 2 : hybrid = SP2; break;
        case 3 :
          // Check for potential SP2. If only 1 of the bonded atoms is
          // hydrogen, assume SP2. TODO actually check for aromaticity.
          int n_hydrogens = 0;
          for (Atom::bond_iterator bat = AJ.bondbegin(); bat != AJ.bondend(); ++bat)
            if (topIn[*bat].Element() == Atom::HYDROGEN)
              n_hydrogens++;
          if (n_hydrogens == 1)
            hybrid = SP2;
          else
            hybrid = SP3;
          break;
      }
      break;
    case Atom::OXYGEN :
      switch (AJ.Nbonds()) {
        case 2 : hybrid = SP3; break;
      }
      break;
    case Atom::SULFUR :
      switch (AJ.Nbonds()) {
        case 2 : hybrid = SP3; break;
      }
      break;
    default: hybrid = UNKNOWN_HYBRIDIZATION; break;
  }
  // Fill in what values we can for known atoms
/*  std::vector<double> knownTheta( AJ.Nbonds() );
  int knownIdx = -1;
  for (int idx = 0; idx < AJ.Nbonds(); idx++) {
    int atnum = AJ.Bond(idx);
    if (atnum != ak && atomPositionKnown[atnum]) {
      knownTheta[idx] = CalcAngle(frameIn.XYZ(atnum),
                                  frameIn.XYZ(aj),
                                  frameIn.XYZ(ak));
      mprintf("DEBUG:\t\tKnown theta for %s = %g\n", topIn.AtomMaskName(atnum).c_str(), knownTheta[idx]*Constants::RADDEG);
      if (knownIdx == -1) knownIdx = idx; // FIXME handle more than 1 known
    }
  }
  if (knownIdx == -1) {*/
    //mprintf("DEBUG:\t\tNo known theta.\n");
  if (hybrid == UNKNOWN_HYBRIDIZATION) {
    // Assign a theta based on number of bonds 
    switch (AJ.Nbonds()) {
      case 4 : hybrid = SP3; break;
      case 3 : hybrid = SP2; break;
      case 2 : hybrid = SP; break;
      default : mprinterr("Internal Error: AssignTheta(): Unhandled # bonds for %s (%i)\n", topIn.AtomMaskName(aj).c_str(), AJ.Nbonds()); return 1;
    }
  }
  // Assign a theta based on hybridization
  switch (hybrid) {
    case SP3 : theta = 109.5 * Constants::DEGRAD; break;
    case SP2 : theta = 120.0 * Constants::DEGRAD; break;
    case SP  : theta = 180.0 * Constants::DEGRAD; break;
    default : mprinterr("Internal Error: AssignTheta(): Unhandled hybridization for %s (%i)\n", topIn.AtomMaskName(aj).c_str(), AJ.Nbonds()); return 1;
  }/*
  } else {
    theta = knownTheta[knownIdx]; // TODO just use above guess via hybrid?
  }*/

  return 0;
}

/// Recursive function to return depth from an atom along bonds
static int atom_depth(int& depth,
                      int at, Topology const& topIn, std::vector<bool>& visited, int maxdepth)
{
  if (depth == maxdepth) return 0;
  depth++;
  visited[at] = true;
  int depthFromHere = 1;
  for (Atom::bond_iterator bat = topIn[at].bondbegin(); bat != topIn[at].bondend(); ++bat)
  {
    if (!visited[*bat])
      depthFromHere += atom_depth( depth, *bat, topIn, visited, maxdepth );
  }
  return depthFromHere;
}

static inline double wrap360(double phi) {
  if (phi > Constants::PI)
    return phi - Constants::TWOPI;
  else if (phi < -Constants::PI)
    return phi + Constants::TWOPI;
  else
    return phi;
}

/** Attempt to assign a reasonable value for phi internal coordinate for atom i
  * given that atoms j k and l have known positions.
  *   j - k
  *  /     \
  * i       l
  */
int Cpptraj::Structure::Model::AssignPhi(double& phi, int ai, int aj, int ak, int al,
                                         Topology const& topIn, Frame const& frameIn,
                                         std::vector<bool> const& atomPositionKnown,
                                         BuildAtom const& AtomJ, int debug)
{
  // Figure out hybridization and chirality of atom j.
  if (debug > 0)
    mprintf("DEBUG: AssignPhi for atom j : %s\n", topIn.AtomMaskName(aj).c_str());

  Atom const& AJ = topIn[aj];
  if (debug > 0) mprintf("DEBUG:\t\tNbonds: %i\n", AJ.Nbonds());
  // If atom J only has 2 bonds, ai-aj-ak-al is the only possibility.
  if (AJ.Nbonds() < 3) {
    if (debug > 0)
      mprintf("DEBUG:\t\tFewer than 3 bonds. Setting phi to -180.\n");
    phi = -180 * Constants::DEGRAD;
    return 0;
  }

  // TODO check that atom i actually ends up on the list?
  std::vector<int> const& priority = AtomJ.Priority();
  if (debug > 0) {
    mprintf("DEBUG: Original chirality around J %s is %s\n", topIn.AtomMaskName(aj).c_str(), chiralStr(AtomJ.Chirality()));
    mprintf("DEBUG:\t\tPriority around J %s(%i):", 
            topIn.AtomMaskName(aj).c_str(), (int)atomPositionKnown[aj]);
    for (int idx = 0; idx < AJ.Nbonds(); idx++)
      mprintf(" %s(%i)", topIn.AtomMaskName(priority[idx]).c_str(), (int)atomPositionKnown[priority[idx]]);
    mprintf("\n");
  }

  // Fill in what values we can for known atoms
  std::vector<double> knownPhi( AJ.Nbonds() );
  int knownIdx = -1;
  for (int idx = 0; idx < AJ.Nbonds(); idx++) {
    int atnum = priority[idx];
    if (atnum != ak && atomPositionKnown[atnum] &&
                       atomPositionKnown[aj] &&
                       atomPositionKnown[ak] &&
                       atomPositionKnown[al])
    {
      knownPhi[idx] = Torsion(frameIn.XYZ(atnum),
                              frameIn.XYZ(aj),
                              frameIn.XYZ(ak),
                              frameIn.XYZ(al));
      if (debug > 0)
        mprintf("DEBUG:\t\tKnown phi for %s = %g\n", topIn.AtomMaskName(atnum).c_str(), knownPhi[idx]*Constants::RADDEG);
      if (knownIdx == -1) knownIdx = idx; // FIXME handle more than 1 known
    }
  }

  // If we have to assign an initial phi, make trans the longer branch
  if (knownIdx == -1) {
    std::vector<bool> visited = atomPositionKnown;
    // TODO: Ensure bonded atoms are not yet visited?
    visited[aj] = true;
    visited[ak] = true;
    std::vector<int> depth( AJ.Nbonds() );
    for (int idx = 0; idx < AJ.Nbonds(); idx++) {
      int atnum = priority[idx];
      if (atnum != ak) {
        int currentDepth = 0;
        depth[idx] = atom_depth(currentDepth, atnum, topIn, visited, 10);
        if (debug > 0)
          mprintf("DEBUG:\t\tAJ %s depth from %s is %i\n",
                  topIn.AtomMaskName(aj).c_str(), topIn.AtomMaskName(atnum).c_str(), depth[idx]);
        if (knownIdx == -1 && depth[idx] < 3) {
          knownIdx = idx;
          knownPhi[idx] = 0;
        }
      }
    }
  }

  // The interval will be 360 / (number of bonds - 1)
  double interval = Constants::TWOPI / (AJ.Nbonds() - 1);

  double startPhi;
  if (knownIdx == -1) {
    startPhi = -180*Constants::DEGRAD;
    if (debug > 0) mprintf("DEBUG:\t\tNo known phi. Setting to %g.\n", startPhi*Constants::RADDEG);
    knownIdx = 0;
  } else
    startPhi = knownPhi[knownIdx];

  if (AtomJ.Chirality() == IS_R) {
    startPhi = -startPhi;
    interval = -interval;
  }
  if (debug > 0) {
    mprintf("DEBUG:\t\tStart phi is %g degrees\n", startPhi*Constants::RADDEG);
    mprintf("DEBUG:\t\tInterval is %g degrees\n", interval * Constants::RADDEG);
  }
    
  // Forward direction
  double currentPhi = startPhi;
  for (int idx = knownIdx; idx < AJ.Nbonds(); idx++) {
    int atnum = priority[idx];
    if (atnum != ak) {
      if (atnum == ai) phi = currentPhi;
      if (debug > 0)
        mprintf("DEBUG:\t\t\t%s phi= %g\n", topIn.AtomMaskName(atnum).c_str(), currentPhi*Constants::RADDEG);
      currentPhi += interval;
      currentPhi = wrap360(currentPhi);
    }
  }
  // Reverse direction
  currentPhi = startPhi - interval;
  for (int idx = knownIdx - 1; idx > -1; idx--) {
    int atnum = priority[idx];
    if (atnum != ak) {
      if (atnum == ai) phi = currentPhi;
      if (debug > 0)
        mprintf("DEBUG:\t\t\t%s phi= %g\n", topIn.AtomMaskName(atnum).c_str(), currentPhi*Constants::RADDEG);
      currentPhi -= interval;
      currentPhi = wrap360(currentPhi);
    }
  }

  return 0;
}
