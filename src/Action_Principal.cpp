#include "Action_Principal.h"
#include "CpptrajStdio.h"
// DEBUG
#include "ptraj_actions.h"
#include "Constants.h"
#include "vectormath.h"
#include "Matrix_3x3.h"

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

// DEBUG
static void CheckAngle( char axis, char eigen, double *Evec) {
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
}

int Action_Principal::action() {
  double Inertia[9], CXYZ[3], Eval[3];

  currentFrame->CalculateInertia( mask_, Inertia, CXYZ );
  printMatrix_3x3("CPPTRAJ INERTIA", Inertia);

  // DEBUG - General test
  double Evec[9];
  Matrix_3x3 TEMP( Inertia );
  TEMP.Diagonalize( Eval, Evec );
  printVector("GENERAL EVAL", Eval );
  TEMP.Print("GENERAL");
  printMatrix_3x3("GENERAL EVEC", Evec);
  
  int i1,i2,i3;
  /*double TEMP[9];
  for (int i = 0; i < 9; ++i)
    TEMP[i] = Inertia[i];
  General_.Diagonalize( TEMP, Eval );
  printVector("GENERAL EVAL", Eval );*/
  if (Eval[0] > Eval[1] && Eval[0] > Eval[2]) { // 0 is max
    if (Eval[1] > Eval[2]) {
      i1 = 0; i2 = 1; i3 = 2;
    } else {
      i1 = 0; i2 = 2; i3 = 1;
    }
  } else if (Eval[1] > Eval[0] && Eval[1] > Eval[2]) { // 1 is max
    if (Eval[0] > Eval[2]) {
      i1 = 1; i2 = 0; i3 = 2;
    } else {
      i1 = 1; i2 = 2; i3 = 0;
    }
  } else if (Eval[0] > Eval[1]) { // 2 is max
    i1 = 2; i2 = 0; i3 = 1;
  } else {
    i1 = 2; i2 = 1; i3 = 0;
  }
  mprintf("EIGENVALUE ORDER (0=high, 3=med, 6=low): %i %i %i\n",i1,i2,i3);

  // Swap Eigenvectors
  // NOTE: Eigenvectors are in columns. Place them in rows in order to
  //       effect an inverse rotation.
  Inertia[0] = Evec[i1];
  Inertia[1] = Evec[i1+3];
  Inertia[2] = Evec[i1+6];

  Inertia[3] = Evec[i2];
  Inertia[4] = Evec[i2+3];
  Inertia[5] = Evec[i2+6];

  Inertia[6] = Evec[i3]; 
  Inertia[7] = Evec[i3+3]; 
  Inertia[8] = Evec[i3+6];

  // Swap eigenvalues
  CXYZ[0] = Eval[i1]; 
  CXYZ[1] = Eval[i2]; 
  CXYZ[2] = Eval[i3];

  // Invert eigenvector signs to avoid reflections
  if (i1 == 0 && i2 == 2 && i3 == 1) {
    Inertia[3] = -Inertia[3];
    Inertia[4] = -Inertia[4];
    Inertia[5] = -Inertia[5];
  } else if (i1 == 2 && i2 == 0 && i3 == 1) {
    Inertia[0] = -Inertia[0];
    Inertia[1] = -Inertia[1];
    Inertia[2] = -Inertia[2];
    Inertia[3] = -Inertia[3];
    Inertia[4] = -Inertia[4];
    Inertia[5] = -Inertia[5];
    Inertia[6] = -Inertia[6];
    Inertia[7] = -Inertia[7];
    Inertia[8] = -Inertia[8];
  }

  printVector("SWAPPED EVALs", CXYZ);
  printMatrix_3x3("SWAPPED EVECs (rows)", Inertia);

  // END DEBUG

  //Principal_.Diagonalize( Inertia, Eval );

  // Rotate - since Inertia is already transposed (eigenvectors
  // are returned in rows since Principal_.Diagonalize calls
  // FORTRAN LAPACK routine dsyev just do plain rotation.
  // Ordering of vectors however is ascending rather than descending
  // so swap that.
/*  CXYZ[0] = Inertia[0];
  CXYZ[1] = Inertia[1];
  CXYZ[2] = Inertia[2];
  Inertia[0] = Inertia[6];
  Inertia[1] = Inertia[7];
  Inertia[2] = Inertia[8];
  Inertia[6] = CXYZ[0];
  Inertia[7] = CXYZ[1];
  Inertia[8] = CXYZ[2];*/

  // DEBUG
/*
  double CHECK[3][3];
  CHECK[0][0] = Inertia[0];
  CHECK[0][1] = Inertia[1];
  CHECK[0][2] = Inertia[2];
  CHECK[1][0] = Inertia[3];
  CHECK[1][1] = Inertia[4];
  CHECK[1][2] = Inertia[5];
  CHECK[2][0] = Inertia[6];
  CHECK[2][1] = Inertia[7];
  CHECK[2][2] = Inertia[8];
  if ( jacobiCheckChirality( Eval, CHECK ) !=0 ) {
    mprintf("PRINCIPAL: WARNING!!! CHECK CHIRALITY: vectors swapped!\n");
    Inertia[1] = CHECK[0][1];
    Inertia[4] = CHECK[1][1];
    Inertia[7] = CHECK[2][1];
  } 
*/
  // DEBUG - Calc Angles between eigenvectors and axis
  // 0 1 2 is X eigenvector
  // 3 4 5 is Y eigenvector
  // 6 7 8 is Z eigenvector
  /*CheckAngle('x','x',Inertia);
  CheckAngle('x','y',Inertia);
  CheckAngle('x','z',Inertia);
  CheckAngle('y','x',Inertia);
  CheckAngle('y','y',Inertia);
  CheckAngle('y','z',Inertia);
  CheckAngle('z','x',Inertia);
  CheckAngle('z','y',Inertia);
  CheckAngle('z','z',Inertia);*/
  mprintf("EIGENVECTOR MATRIX ANGLE = %lf\n", matrix_to_angle( Inertia )*RADDEG );
  mprintf("\n");
/*
  CXYZ[0] = 1; CXYZ[1] = 0; CXYZ[2] = 0; // Xvec
  double angle = dot_product_angle(CXYZ, Inertia);
  mprintf("PRINCIPAL DEBUG: Angle between X and eigenX is %lf\n", angle*RADDEG);
  CXYZ[0] = 0; CXYZ[1] = 1; CXYZ[2] = 0; // Yvec
  angle = dot_product_angle(CXYZ, Inertia+3);
  mprintf("PRINCIPAL DEBUG: Angle between Y and eigenY is %lf\n", angle*RADDEG);
  CXYZ[0] = 0; CXYZ[1] = 0; CXYZ[2] = 1; // Zvec
  angle = dot_product_angle(CXYZ, Inertia+6);
  mprintf("PRINCIPAL DEBUG: Angle between Z and eigenZ is %lf\n\n", angle*RADDEG);
*/

  if (doRotation_)
    currentFrame->Rotate( Inertia );

  return 0;
}


