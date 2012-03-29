#include <cmath> // cos, sin, sqrt
#include "Box.h"
#include "Constants.h" // DEGRAD
#include "CpptrajStdio.h"

// CONSTRUCTOR
Box::Box() :
  debug_(0),
  btype_(NOBOX)
{
  box_[0] = 0;
  box_[1] = 0;
  box_[2] = 0;
  box_[3] = 0;
  box_[4] = 0;
  box_[5] = 0;
}

// COPY CONSTRUCTOR
Box::Box(const Box& rhs) : 
  debug_(rhs.debug_),
  btype_(rhs.btype_)
{
  box_[0] = rhs.box_[0];
  box_[1] = rhs.box_[1];
  box_[2] = rhs.box_[2];
  box_[3] = rhs.box_[3];
  box_[4] = rhs.box_[4];
  box_[5] = rhs.box_[5];
}

// Assignment
Box &Box::operator=(const Box& rhs) {
  if (&rhs == this) return *this;
  debug_ = rhs.debug_;
  btype_ = rhs.btype_;
  box_[0] = rhs.box_[0];
  box_[1] = rhs.box_[1];
  box_[2] = rhs.box_[2];
  box_[3] = rhs.box_[3];
  box_[4] = rhs.box_[4];
  box_[5] = rhs.box_[5];
  return *this;
}

const double Box::TRUNCOCTBETA = 109.4712206344906917365733534097672;

void Box::PrintBoxType() {
  switch (btype_) {
    case NOBOX : mprintf("No Box"); break;
    case ORTHO : mprintf("Ortho. Box"); break;
    case NONORTHO: mprintf("Non-ortho. Box"); break;
  }
}

// Box::SetBetaLengths()
/** Expected array format: {OLDBETA, BOX(1), BOX(2), BOX(3)}
  */
void Box::SetBetaLengths(std::vector<double> &betaXYZ) {
  if (betaXYZ.size() != 4) {
    mprinterr("Error: SetLengthsBeta: Size of array is not 4 (%zu)\n",betaXYZ.size());
    return;
  }
  box_[0] = betaXYZ[1];
  box_[1] = betaXYZ[2];
  box_[2] = betaXYZ[3];
  box_[3] = 0;
  box_[4] = betaXYZ[0];
  box_[5] = 0;
  SetBoxType();
}

// Box::SetTruncOct()
/** Set as truncated octahedron with no lengths. */
void Box::SetTruncOct() {
  box_[0] = 0;
  box_[1] = 0;
  box_[2] = 0;
  box_[3] = TRUNCOCTBETA;
  box_[4] = TRUNCOCTBETA;
  box_[5] = TRUNCOCTBETA;
  btype_ = NONORTHO;
}

// Box::SetBoxType()
/// Determine box type (none/ortho/nonortho) based on box angles.
void Box::SetBoxType() {
  btype_ = NONORTHO;
  // No angles, no box
  if ( box_[3] <= 0 && box_[4] <= 0 && box_[5] <= 0)
    btype_ = NOBOX;
  // Determine orthogonal / non-orthogonal from angles
  else if (box_[3] == 90.0 && box_[4] == 90.0 && box_[5] == 90.0)
    btype_ = ORTHO;
  else if (box_[3] == 0 && box_[4] != 0 && box_[5] == 0) {
    // Only beta angle is set (e.g. from Amber topology).
    if (box_[4] == 90.0) {
      btype_ = ORTHO;
      box_[3] = 90.0;
      box_[5] = 90.0;
      if (debug_>0)
        mprintf("    Setting box to be orthogonal\n");
    } else if (box_[4] > 109.47 && box_[4] < 109.48) {
      btype_ = NONORTHO;
      //Box[3] = TRUNCOCTBETA;
      box_[3] = box_[4];
      box_[5] = box_[4];
      if (debug_>0) 
        mprintf("    Setting box to be a truncated octahedron, angle is %lf\n",box_[3]);
    } else if (box_[4] == 60.0) {
      btype_ = NONORTHO;
      box_[3]=60.0; 
      box_[4]=90.0; 
      box_[5]=60.0;
      if (debug_>0)
        mprintf("    Setting box to be a rhombic dodecahedron, alpha=gamma=60.0, beta=90.0\n");
    }
  } 
  //if (debug_>0) 
    mprintf("\tBox type is %i (beta=%lf)\n",(int)btype_, box_[4]);
}

// Box::AmberIfbox()
/** Based on beta, return Amber IFBOX type:
  *   0: No box
  *   1: Box
  *   2: Truncated octahedral box
  */
int Box::AmberIfbox() {
  if (box_[4] <= 0.0)
    return 0;
  else if (box_[4] > 109.47 && box_[4] < 109.48)
    return 2;
  else
    return 1;
}

// Box::ToRecip()
/** Use box coordinates to calculate reciprocal space conversions for use
  * with imaging routines. Return cell volume.
  */
// NOTE: Move to separate routine in DistRoutines?
double Box::ToRecip(double *ucell, double *recip) {
  double u12x,u12y,u12z;
  double u23x,u23y,u23z;
  double u31x,u31y,u31z;
  double volume,onevolume;

  ucell[0] = box_[0]; // ucell(1,1)
  ucell[1] = 0.0;    // ucell(2,1)
  ucell[2] = 0.0;    // ucell(3,1)
  ucell[3] = box_[1]*cos(DEGRAD*box_[5]); // ucell(1,2)
  ucell[4] = box_[1]*sin(DEGRAD*box_[5]); // ucell(2,2)
  ucell[5] = 0.0;                       // ucell(3,2)
  ucell[6] = box_[2]*cos(DEGRAD*box_[4]);                                         // ucell(1,3)
  ucell[7] = (box_[1]*box_[2]*cos(DEGRAD*box_[3]) - ucell[6]*ucell[3]) / ucell[4]; // ucell(2,3)
  ucell[8] = sqrt(box_[2]*box_[2] - ucell[6]*ucell[6] - ucell[7]*ucell[7]);       // ucell(3,3)

  // Get reciprocal vectors
  u23x = ucell[4]*ucell[8] - ucell[5]*ucell[7];
  u23y = ucell[5]*ucell[6] - ucell[3]*ucell[8];
  u23z = ucell[3]*ucell[7] - ucell[4]*ucell[6];
  u31x = ucell[7]*ucell[2] - ucell[8]*ucell[1];
  u31y = ucell[8]*ucell[0] - ucell[6]*ucell[2];
  u31z = ucell[6]*ucell[1] - ucell[7]*ucell[0];
  u12x = ucell[1]*ucell[5] - ucell[2]*ucell[4];
  u12y = ucell[2]*ucell[3] - ucell[0]*ucell[5];
  u12z = ucell[0]*ucell[4] - ucell[1]*ucell[3];
  volume=ucell[0]*u23x + ucell[1]*u23y + ucell[2]*u23z;
  onevolume = 1.0 / volume;

  recip[0] = u23x*onevolume;
  recip[1] = u23y*onevolume;
  recip[2] = u23z*onevolume;
  recip[3] = u31x*onevolume;
  recip[4] = u31y*onevolume;
  recip[5] = u31z*onevolume;
  recip[6] = u12x*onevolume;
  recip[7] = u12y*onevolume;
  recip[8] = u12z*onevolume;

  return volume;
}

