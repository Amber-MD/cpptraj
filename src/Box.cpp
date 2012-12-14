#include <cmath> // cos, sin, sqrt
#include "Box.h"
#include "Constants.h" // DEGRAD
#include "CpptrajStdio.h"

// CONSTRUCTOR
Box::Box() : btype_(NOBOX) //, debug_(0)
{
  box_[0] = 0;
  box_[1] = 0;
  box_[2] = 0;
  box_[3] = 0;
  box_[4] = 0;
  box_[5] = 0;
}

Box::Box(const double* bIn) //: debug_(0)
{
  SetBox( bIn );
}

// COPY CONSTRUCTOR
Box::Box(const Box& rhs) : btype_(rhs.btype_) //, debug_(rhs.debug_)
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
  //debug_ = rhs.debug_;
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

const char* Box::BoxNames[] = {
  "None", "Orthogonal", "Trunc. Oct.", "Rhombic Dodec.", "Non-orthogonal"
};

// Box::TypeName()
const char* Box::TypeName() const {
  return BoxNames[btype_];
}

// Box::SetBetaLengths()
void Box::SetBetaLengths(double beta, double xin, double yin, double zin) {
  box_[0] = xin;
  box_[1] = yin;
  box_[2] = zin;
  box_[3] = 0;
  box_[4] = beta;
  box_[5] = 0;
  SetBoxType();
}

/** Set box from double[6] array */
void Box::SetBox(const double* xyzabg) {
  if (xyzabg == 0) {
    mprinterr("Error: SetBox: Input array is null\n");
    return;
  }
  box_[0] = xyzabg[0];
  box_[1] = xyzabg[1];
  box_[2] = xyzabg[2];
  box_[3] = xyzabg[3];
  box_[4] = xyzabg[4];
  box_[5] = xyzabg[5];
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
  btype_ = TRUNCOCT;
}

// Box::SetNoBox()
void Box::SetNoBox() {
  box_[0] = 0;
  box_[1] = 0;
  box_[2] = 0;
  box_[3] = 0;
  box_[4] = 0;
  box_[5] = 0;
  btype_ = NOBOX;
}

// Box::SetMissingInfo()
/** Set this box info from rhs if <= 0. */
void Box::SetMissingInfo(const Box& rhs) {
  if (box_[0] <= 0) box_[0] = rhs.box_[0];
  if (box_[1] <= 0) box_[1] = rhs.box_[1];
  if (box_[2] <= 0) box_[2] = rhs.box_[2];
  if (box_[3] <= 0) box_[3] = rhs.box_[3];
  if (box_[4] <= 0) box_[4] = rhs.box_[4];
  if (box_[5] <= 0) box_[5] = rhs.box_[5];
  SetBoxType();
}

static bool IsTruncOct(double angle) {
  if (angle > 109.47 && angle < 109.48) return true;
  return false;
}

// Box::SetBoxType()
/** Determine box type (none/ortho/nonortho) based on box angles. */
void Box::SetBoxType() {
  btype_ = NONORTHO;
  // No angles, no box
  if ( box_[3] <= 0 && box_[4] <= 0 && box_[5] <= 0)
    btype_ = NOBOX;
  // All 90, orthogonal 
  else if (box_[3] == 90.0 && box_[4] == 90.0 && box_[5] == 90.0)
    btype_ = ORTHO;
  // All 109.47, truncated octahedron
  else if ( IsTruncOct( box_[3] ) && IsTruncOct( box_[4] ) && IsTruncOct( box_[5] ) )
    btype_ = TRUNCOCT;
  else if (box_[3] == 0 && box_[4] != 0 && box_[5] == 0) {
    // Only beta angle is set (e.g. from Amber topology).
    if (box_[4] == 90.0) {
      btype_ = ORTHO;
      box_[3] = 90.0;
      box_[5] = 90.0;
      //if (debug_>0) mprintf("\tSetting box to be orthogonal\n");
    } else if ( IsTruncOct( box_[4] ) ) {
      btype_ = TRUNCOCT;
      //Box[3] = TRUNCOCTBETA;
      box_[3] = box_[4];
      box_[5] = box_[4];
      //if (debug_>0) mprintf("\tSetting box to be a truncated octahedron, angle is %lf\n",box_[3]);
    } else if (box_[4] == 60.0) {
      btype_ = RHOMBIC;
      box_[3]=60.0; 
      box_[4]=90.0; 
      box_[5]=60.0;
      //if (debug_>0) mprintf("\tSetting box to be a rhombic dodecahedron, alpha=gamma=60.0, beta=90.0\n");
    } else {
      mprintf("Warning: Box: Unrecognized beta (%lf); setting all angles to beta.\n",box_[4]);
      box_[3] = box_[4];
      box_[5] = box_[4];
    }
  } 
  //if (debug_>0) mprintf("\tBox type is %s (beta=%lf)\n",TypeName(), box_[4]);
}

// Box::ToRecip()
/** Use box coordinates to calculate reciprocal space conversions for use
  * with imaging routines. Return cell volume.
  */
double Box::ToRecip(Matrix_3x3& ucell, Matrix_3x3& recip) const {
  double u12x,u12y,u12z;
  double u23x,u23y,u23z;
  double u31x,u31y,u31z;
  double volume,onevolume;

  ucell[0] = box_[0]; // ucell(1,1)
  ucell[1] = 0.0;     // ucell(2,1)
  ucell[2] = 0.0;     // ucell(3,1)
  ucell[3] = box_[1]*cos(DEGRAD*box_[5]); // ucell(1,2)
  ucell[4] = box_[1]*sin(DEGRAD*box_[5]); // ucell(2,2)
  ucell[5] = 0.0;                         // ucell(3,2)
  ucell[6] = box_[2]*cos(DEGRAD*box_[4]);                                          // ucell(1,3)
  ucell[7] = (box_[1]*box_[2]*cos(DEGRAD*box_[3]) - ucell[6]*ucell[3]) / ucell[4]; // ucell(2,3)
  ucell[8] = sqrt(box_[2]*box_[2] - ucell[6]*ucell[6] - ucell[7]*ucell[7]);        // ucell(3,3)

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
