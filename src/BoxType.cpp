// BoxType
#include "BoxType.h"
#include <cstddef>
#include "CpptrajStdio.h"

#define TRUNCOCTBETA 109.4712206344906917365733534097672

/* CheckBoxType()
 * Determine box type (none/ortho/nonortho) based on box angles.
 */
BoxType CheckBoxType(double *box, int debug) {
  BoxType btype = NOBOX;
  if (box==NULL) return NOBOX;

  // Determine orthogonal / non-orthogonal from angles
  if (box[0]==0.0 || box[1]==0.0 || box[2]==0.0)
    btype = NOBOX;
  else if (box[0]==90.0 && box[1]==90.0 && box[2]==90.0)
    btype = ORTHO;
  else
    btype = NONORTHO;
  if (debug>0) mprintf("    Box type is %i (beta=%lf)\n",btype,box[0]);
  return btype;
}

/* AmberIfbox()
 * Based on beta, return Amber IFBOX type:
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

/* SetBoxInfo()
 * Given an angle (beta) and 3 edge lengths, set up a box 
 * (3 lengths + 3 angles) and return box type. Also set the
 * ifbox flag if not NULL (0=no box, 1 = std periodic, 
 * 2 = truncated oct.
 * Currently recognized betas:
 *   90.00 - Orthogonal
 *  109.47 - Truncated octahedral
 *   60.00 - Rhombic dodecahedron
 * Any other betas just sets all angles to beta and a warning is printed.
 */
BoxType SetBoxInfo(double *bIn, double *Box, int debug) {
  int ifbox=0;
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
    ifbox=0;
    Box[0]=0.0; Box[1]=0.0; Box[2]=0.0;
    Box[3]=0.0; Box[4]=0.0; Box[5]=0.0;
  } else if (beta == 90.0) {
    btype=ORTHO;
    ifbox=1;
    Box[0]=bx; Box[1]=by; Box[2]=bz;
    Box[3]=90.0; Box[4]=90.0; Box[5]=90.0;
    if (debug>0) mprintf("    Setting box to be orthogonal.\n");
  } else if (beta > 109.47 && beta < 109.48) {
    btype=NONORTHO;
    ifbox=2;
    Box[0]=bx; Box[1]=by; Box[2]=bz;
    //Box[3] = TRUNCOCTBETA;
    Box[3] = beta;
    Box[4]=Box[3]; Box[5]=Box[3];
    if (debug>0) mprintf("    Setting box to be a truncated octahedron, angle is %lf\n",Box[3]);
  } else if (beta == 60.0) {
    btype=NONORTHO;
    ifbox=1;
    Box[0]=bx; Box[1]=by; Box[2]=bz;
    Box[3]=60.0; Box[4]=90.0; Box[5]=60.0;
    if (debug>0)
      mprintf("    Setting box to be a rhombic dodecahedron, alpha=gamma=60.0, beta=90.0\n");
  } else {
    btype=NONORTHO;
    ifbox=1;
    Box[0]=bx; Box[1]=by; Box[2]=bz;
    Box[3]=beta; Box[4]=beta; Box[5]=beta;
    mprintf("    Warning: Unrecognized box type, beta is %lf\n",beta);
  }
  return btype;
}

/* SetChamberBoxInfo()
 * Given the original value of IFBOX pulled from the topology file,
 * this method will set up the Box angles according to
 *   ifbox == 0 --> NOBOX
 *   ifbox == 1 --> ORTHO
 *   ifbox == 2 --> NONORTHO
 * Because the box information doesn't exist in the topology file,
 * there are no box dimensions, so these will just be set to 0
 */
BoxType SetBoxInfoFromIfbox (int orig_ifbox, double *Box, int debug) {

  BoxType btype = NOBOX;

  // Since we don't know the box dimensions, just zero them
  Box[0] = 0.0; Box[1] = 0.0; Box[2] = 0.0;

  if (orig_ifbox == 0) {
    // No box
    btype = NOBOX;
    Box[3] = 0.0; Box[4] = 0.0; Box[5] = 0.0;
  } else if (orig_ifbox == 1) {
    // Orthogonal box
    btype = ORTHO;
    Box[3] = 90.0; Box[4] = 90.0; Box[5] = 90.0;
  } else if (orig_ifbox == 1) {
    // Non-orthogonal box. Assume truncated octahedron, since that's
    // all LEaP produces by default
    btype = NONORTHO;
    Box[3] = TRUNCOCTBETA; Box[4] = TRUNCOCTBETA; Box[5] = TRUNCOCTBETA;
  } else {
    // This should never be hit, so print out a warning and set no box
    mprintf("    Warning: Unrecognized ifbox value %i. Assuming no box\n", orig_ifbox);
    btype = NOBOX;
    Box[3] = 0.0; Box[4] = 0.0; Box[5] = 0.0; 
  }

  return btype;
}
