#include "Model.h"
#include "../CpptrajStdio.h"
#include "../Topology.h"
#include "../TorsionRoutines.h"
#include "../Constants.h"

/** Attempt to assign a reasonable value for theta internal coordinate for
  * atom i given that atoms j and k have known positions.
  */
int Cpptraj::Structure::Model::AssignTheta(double& theta, int ai, int aj, int ak, Topology const& topIn, Frame const& frameIn, std::vector<bool> const& atomPositionKnown)
{
  // Figure out hybridization and chirality of atom j.
  mprintf("DEBUG: AssignTheta for atom j : %s\n", topIn.AtomMaskName(aj).c_str());

  Atom const& AJ = topIn[aj];
  mprintf("DEBUG:\t\tNbonds: %i\n", AJ.Nbonds());
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
    // Assign a theta based on hybridization
    switch (AJ.Nbonds()) {
      case 4 : theta = 109.5 * Constants::DEGRAD; break;
      case 3 : theta = 120.0 * Constants::DEGRAD; break;
      case 2 : theta = 180.0 * Constants::DEGRAD; break;
      default : mprinterr("Internal Error: AssignTheta(): Unhandled # bonds for %s (%i)\n", topIn.AtomMaskName(aj).c_str(), AJ.Nbonds()); return 1;
    }/*
  } else {
    theta = knownTheta[knownIdx]; // TODO just use above guess via hybrid?
  }*/

  return 0;
}

/** Attempt to assign a reasonable value for phi internal coordinate for atom i
  * given that atoms j k and l have known positions.
  *   j - k
  *  /     \
  * i       l
  */
int Cpptraj::Structure::Model::AssignPhi(double& phi, int ai, int aj, int ak, int al, Topology const& topIn, Frame const& frameIn, std::vector<bool> const& atomPositionKnown, std::vector<Cpptraj::Chirality::ChiralType> const& atomChirality)
{
  using namespace Cpptraj::Chirality;
  // Figure out hybridization and chirality of atom j.
  mprintf("DEBUG: AssignPhi for atom j : %s\n", topIn.AtomMaskName(aj).c_str());

  Atom const& AJ = topIn[aj];
  mprintf("DEBUG:\t\tNbonds: %i\n", AJ.Nbonds());
  // If atom J only has 2 bonds, ai-aj-ak-al is the only possibility.
  if (AJ.Nbonds() < 3) {
    mprintf("DEBUG:\t\tFewer than 3 bonds. Setting phi to -180.\n");
    phi = -180 * Constants::DEGRAD;
    return 0;
  }

  // FIXME aj ak and al should be known
  // TODO check that atom i actually ends up on the list?
  std::vector<int> Priority( AJ.Nbonds() );
  int* priority = static_cast<int*>( &Priority[0] );
  double tors = 0;
  // DetermineChirality will always put aj as the last index.
  ChiralType chirality = DetermineChirality(tors, priority, aj, topIn, frameIn, 1); // FIXME debug
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
  mprintf("DEBUG:\t\tPriority around J %s(%i) (tors=%g):", 
          topIn.AtomMaskName(aj).c_str(), (int)atomPositionKnown[aj], tors*Constants::RADDEG);
  for (int idx = 0; idx < AJ.Nbonds(); idx++)
    mprintf(" %s(%i)", topIn.AtomMaskName(priority[idx]).c_str(), (int)atomPositionKnown[priority[idx]]);
  mprintf("\n");
  // Determine if chirality is valid
  bool chirality_is_valid = true;
  for (unsigned int idx = 0; idx < Priority.size(); idx++) {
    if (atomPositionKnown[priority[idx]] != atomPositionKnown[aj]) {
      chirality_is_valid = false;
      break;
    }
  }
  mprintf("DEBUG:\t\tChirality is valid: %i\n", (int)chirality_is_valid);

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
      mprintf("DEBUG:\t\tKnown phi for %s = %g\n", topIn.AtomMaskName(atnum).c_str(), knownPhi[idx]*Constants::RADDEG);
      if (knownIdx == -1) knownIdx = idx; // FIXME handle more than 1 known
    }
  }

  // The interval will be 360 / (number of bonds - 1)
  double interval = Constants::TWOPI / (AJ.Nbonds() - 1);

  double startPhi;
  if (knownIdx == -1) {
    startPhi = -180*Constants::DEGRAD;
    // If all atom statuses match the chirality is valid. Have R start -180.
    //if (chirality_is_valid)
    //{
    //  if (chirality == IS_R) // FIXME just always start -180?
    //    startPhi = -180*Constants::DEGRAD;
    //}
    mprintf("DEBUG:\t\tNo known phi. Setting to %g.\n", startPhi*Constants::RADDEG);
    knownIdx = 0;
  } else
    startPhi = knownPhi[knownIdx];

  if (atomChirality[aj] == IS_R) {
    startPhi = -startPhi;
    interval = -interval;
  }
  mprintf("DEBUG:\t\tStart phi is %g degrees\n", startPhi*Constants::RADDEG);
  mprintf("DEBUG:\t\tInterval is %g degrees\n", interval * Constants::RADDEG);
    
  // Forward direction
  double currentPhi = startPhi;
  for (int idx = knownIdx; idx < AJ.Nbonds(); idx++) {
    int atnum = priority[idx];
    if (atnum != ak) {
      if (atnum == ai) phi = currentPhi;
      mprintf("DEBUG:\t\t\t%s phi= %g\n", topIn.AtomMaskName(atnum).c_str(), currentPhi*Constants::RADDEG);
      currentPhi += interval;
    }
  }
  currentPhi = startPhi - interval;
  for (int idx = knownIdx - 1; idx > -1; idx--) {
    int atnum = priority[idx];
    if (atnum != ak) {
      if (atnum == ai) phi = currentPhi;
      mprintf("DEBUG:\t\t\t%s phi= %g\n", topIn.AtomMaskName(atnum).c_str(), currentPhi*Constants::RADDEG);
      currentPhi -= interval;
    }
  }

  return 0;
}
