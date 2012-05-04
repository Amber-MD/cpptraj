#include "Action_Principal.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Principal::Action_Principal() :
  doRotation_(false)
{
  useMass_ = false;
}


int Action_Principal::init() {
  // Keywords
  doRotation_ = actionArgs.hasKey("dorotation");
  useMass_ = actionArgs.hasKey("mass");

  // Masks
  mask_.SetMaskString( actionArgs.getNextMask() );

  mprintf("    PRINCIPAL");
  if (doRotation_)
    mprintf(" with rotation by");
  else
    mprintf(" without rotation by");
  if (useMass_)
    mprintf(" center of mass");
  else
    mprintf(" center of geometry");
  mprintf(", atoms selected by [%s]\n", mask_.MaskString());

  return 0;
}

int Action_Principal::setup() {
  if (currentParm->SetupIntegerMask(mask_)) return 1;

  if (mask_.None()) {
    mprintf("Warning: No atoms selected for %s [%s].\n",currentParm->c_str(), mask_.MaskString());
    return 1;
  }

  mprintf("\tSelected %i atoms.\n", mask_.Nselected());
  return 0;
}

int Action_Principal::action() {
  double Inertia[9], CXYZ[3], Eval[3];
  double CHECK[3][3];

  currentFrame->CalculateInertia( mask_, Inertia, CXYZ );

  Principal_.Diagonalize( Inertia, Eval );

  // Rotate - since Inertia is already transposed (eigenvectors
  // are returned in rows since Principal_.Diagonalize calls
  // FORTRAN LAPACK routine dsyev just do plain rotation.
  // Ordering of vectors however is ascending rather than descending
  // so swap that.
  Eval[0] = Inertia[0];
  Eval[1] = Inertia[1];
  Eval[2] = Inertia[2];
  Inertia[0] = Inertia[6];
  Inertia[1] = Inertia[7];
  Inertia[2] = Inertia[8];
  Inertia[6] = Eval[0];
  Inertia[7] = Eval[1];
  Inertia[8] = Eval[2];

  // DEBUG
  

  if (doRotation_)
    currentFrame->Rotate( Inertia );

  return 0;
}

