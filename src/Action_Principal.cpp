#include "Action_Principal.h"
#include "CpptrajStdio.h"
// DEBUG
//#incl ude "Constants.h"
#include "vectormath.h"
#include "Matrix_3x3.h"

// CONSTRUCTOR
Action_Principal::Action_Principal() :
  doRotation_(false)
{
  useMass_ = false;
}

// Action_Principal::init()
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

// Action_Principal::setup()
int Action_Principal::setup() {
  if (currentParm->SetupIntegerMask(mask_)) return 1;

  if (mask_.None()) {
    mprintf("Warning: No atoms selected for %s [%s].\n",currentParm->c_str(), mask_.MaskString());
    return 1;
  }

  mprintf("\tSelected %i atoms.\n", mask_.Nselected());
  return 0;
}

// DEBUG
/*static void CheckAngle( char axis, char eigen, double *Evec) {
  double CXYZ[3];
  double* ptr;
  switch (axis) {
    case 'x' : CXYZ[0] = 1; CXYZ[1] = 0; CXYZ[2] = 0; break;
    case 'y' : CXYZ[0] = 0; CXYZ[1] = 1; CXYZ[2] = 0; break;
    case 'z' : CXYZ[0] = 0; CXYZ[1] = 0; CXYZ[2] = 1; break;
    default: mprinterr("CHECKANGLE ERROR\n"); return;
  }
  switch (eigen) {
    case 'x' : ptr = Evec; break;
    case 'y' : ptr = Evec + 3; break;
    case 'z' : ptr = Evec + 6; break;
    default: mprinterr("CHECKANGLE ERROR\n"); return;
  }
  double angle = dot_product_angle(CXYZ, ptr);
  mprintf("PRINCIPAL DEBUG: Angle between %c and eigen(%c) is %lf\n",axis,eigen,angle*RADDEG);
}*/

// Action_Principal::action()
int Action_Principal::action() {
  double Inertia[9], CXYZ[3], Evec[9], Eval[3];

  currentFrame->CalculateInertia( mask_, Inertia, CXYZ );
  //printMatrix_3x3("PRINCIPAL INERTIA", Inertia);

  Matrix_3x3 TEMP( Inertia );
  // NOTE: Diagonalize_Sort_Chirality places sorted eigenvectors in rows.
  TEMP.Diagonalize_Sort_Chirality( Evec, Eval, debug );
  if (debug > 2) {
    printVector("PRINCIPAL EIGENVALUES", Eval );
    //TEMP.Print("GENERAL");
    printMatrix_3x3("PRINCIPAL EIGENVECTORS (Rows)", Evec);
  }
  
  // Rotate - since Evec is already transposed (eigenvectors
  // are returned in rows) just do plain rotation to affect an
  // inverse rotation.
  if (doRotation_)
    currentFrame->Rotate( Evec );

  return 0;
}

