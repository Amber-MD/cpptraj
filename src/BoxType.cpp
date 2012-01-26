/*! \file BoxType.cpp
    \brief Functions related to box information.
 */
#include "BoxType.h"
#include <cstddef>
#include "CpptrajStdio.h"

// CheckBoxType()
/// Determine box type (none/ortho/nonortho) based on box angles.
BoxType CheckBoxType(double *box, int debug) {
  BoxType btype = NOBOX;
  if (box==NULL || box[0]<0 || box[1]<0 || box[2]<0) 
    btype = NOBOX;
  // Determine orthogonal / non-orthogonal from angles
  else if (box[0]==0.0 || box[1]==0.0 || box[2]==0.0)
    btype = NOBOX;
  else if (box[0]==90.0 && box[1]==90.0 && box[2]==90.0)
    btype = ORTHO;
  else
    btype = NONORTHO;
  if (debug>0) mprintf("\tBox type is %i (beta=%lf)\n",btype,box[1]);
  return btype;
}

// AmberIfbox()
/** Based on beta, return Amber IFBOX type:
  *   0: No box
  *   1: Box
  *   2: Truncated octahedral box
  */
int AmberIfbox(double beta) {
  if (beta <= 0.0)
    return 0;
  else if (beta > 109.47 && beta < 109.48) 
    return 2;
  else
    return 1;
}

// SetBoxInfo()
/** Given an angle (beta) and 3 edge lengths, set up a box 
  * (3 lengths + 3 angles) and return box type. 
  * Currently recognized betas:
  *   90.00 - Orthogonal
  *  109.47 - Truncated octahedral
  *   60.00 - Rhombic dodecahedron
  * Any other betas just sets all angles to beta and a warning is printed.
  */
BoxType SetBoxInfo(double *bIn, double *Box, int debug) {
  BoxType btype = NOBOX;
  double beta, bx, by, bz;

  beta = bIn[0];
  bx = bIn[1];
  by = bIn[2];
  bz = bIn[3];
    
  // Determine box type from beta (none, ortho, non-ortho (truncated oct/triclinic)
  if (beta<=0.0) {
    //if (BoxType>0)
    //  mprintf("    %s: Removing box information.\n",parmName);
    btype=NOBOX;
    Box[0]=0.0; Box[1]=0.0; Box[2]=0.0;
    Box[3]=0.0; Box[4]=0.0; Box[5]=0.0;
  } else if (beta == 90.0) {
    btype=ORTHO;
    Box[0]=bx; Box[1]=by; Box[2]=bz;
    Box[3]=90.0; Box[4]=90.0; Box[5]=90.0;
    if (debug>0) mprintf("    Setting box to be orthogonal.\n");
  } else if (beta > 109.47 && beta < 109.48) {
    btype=NONORTHO;
    Box[0]=bx; Box[1]=by; Box[2]=bz;
    //Box[3] = TRUNCOCTBETA;
    Box[3] = beta;
    Box[4]=Box[3]; Box[5]=Box[3];
    if (debug>0) mprintf("    Setting box to be a truncated octahedron, angle is %lf\n",Box[3]);
  } else if (beta == 60.0) {
    btype=NONORTHO;
    Box[0]=bx; Box[1]=by; Box[2]=bz;
    Box[3]=60.0; Box[4]=90.0; Box[5]=60.0;
    if (debug>0)
      mprintf("    Setting box to be a rhombic dodecahedron, alpha=gamma=60.0, beta=90.0\n");
  } else {
    btype=NONORTHO;
    Box[0]=bx; Box[1]=by; Box[2]=bz;
    Box[3]=beta; Box[4]=beta; Box[5]=beta;
    mprintf("    Warning: Unrecognized box type, beta is %lf\n",beta);
  }
  return btype;
}

