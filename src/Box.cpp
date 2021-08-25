#include <cmath> // cos, sin, sqrt, fabs
#include "Box.h"
#include "Constants.h" // DEGRAD
#include "CpptrajStdio.h"
#include <algorithm> // std::copy

/** CONSTRUCTOR */
Box::Box() :
  cellVolume_(0)
{
  box_[0] = 0;
  box_[1] = 0;
  box_[2] = 0;
  box_[3] = 0;
  box_[4] = 0;
  box_[5] = 0;
}

/** COPY CONSTRUCTOR */
Box::Box(const Box& rhs) :
  unitCell_(rhs.unitCell_),
  fracCell_(rhs.fracCell_),
  cellVolume_(rhs.cellVolume_) //, debug_(rhs.debug_)
{
  box_[0] = rhs.box_[0];
  box_[1] = rhs.box_[1];
  box_[2] = rhs.box_[2];
  box_[3] = rhs.box_[3];
  box_[4] = rhs.box_[4];
  box_[5] = rhs.box_[5];
}

/** Assignment */
Box& Box::operator=(const Box& rhs) {
  if (&rhs == this) return *this;
  box_[0] = rhs.box_[0];
  box_[1] = rhs.box_[1];
  box_[2] = rhs.box_[2];
  box_[3] = rhs.box_[3];
  box_[4] = rhs.box_[4];
  box_[5] = rhs.box_[5];
  unitCell_ = rhs.unitCell_;
  fracCell_ = rhs.fracCell_;
  cellVolume_ = rhs.cellVolume_;
  return *this;
}

/// Swap double precision
static inline void dswap(double& d1, double& d2) {
  double dtemp = d1;
  d1 = d2;
  d2 = dtemp;
}

/** Swap this box with given box. */
void Box::swap(Box& rhs) {
  dswap( box_[0], rhs.box_[0] );
  dswap( box_[1], rhs.box_[1] );
  dswap( box_[2], rhs.box_[2] );
  dswap( box_[3], rhs.box_[3] );
  dswap( box_[4], rhs.box_[4] );
  dswap( box_[5], rhs.box_[5] );
  for (int i=0; i<9; i++) {
    dswap(unitCell_[i], rhs.unitCell_[i]);
    dswap(fracCell_[i], rhs.fracCell_[i]);
  }
  dswap(cellVolume_, rhs.cellVolume_);
}

#ifdef MPI
/** Broadcast box info from Comm master. */
int Box::BroadcastBox(Parallel::Comm const& commIn) {
  commIn.MasterBcast( box_,    6, MPI_DOUBLE );
  unitCell_.BroadcastMatrix( commIn );
  fracCell_.BroadcastMatrix( commIn );
  commIn.MasterBcast( &cellVolume_, 1, MPI_DOUBLE );
  return 0;
}

/** Send box info to recvrank. */
int Box::SendBox(int recvrank, Parallel::Comm const& commIn) const {
  commIn.Send( box_,         6, MPI_DOUBLE, recvrank, 1801 );
  unitCell_.SendMatrix(recvrank, commIn);
  fracCell_.SendMatrix(recvrank, commIn);
  commIn.Send( &cellVolume_, 1, MPI_DOUBLE, recvrank, 1802 );
  return 0;
}

/** Get box info from recvrank. */
int Box::RecvBox(int sendrank, Parallel::Comm const& commIn) {
  commIn.Recv( box_,         6, MPI_DOUBLE, sendrank, 1801 );
  unitCell_.RecvMatrix(sendrank, commIn);
  fracCell_.RecvMatrix(sendrank, commIn);
  commIn.Recv( &cellVolume_, 1, MPI_DOUBLE, sendrank, 1802 );
  return 0;
}
#endif

// -----------------------------------------------------------------------------
const double Box::TRUNCOCTBETA_ = 2.0*acos(1.0/sqrt(3.0))*Constants::RADDEG;

/** This value is chosen so that TRUNCOCTBETA_ - TruncOctDelta_ is just
  * less than 109.47, about 109.4699, since 109.47 is probably the lowest
  * precision angle we can accept. +0.000000703 to avoid FP round-down.
  */
const double Box::TruncOctDelta_ = 0.001220703;

const double Box::TruncOctMin_ = Box::TRUNCOCTBETA_ - Box::TruncOctDelta_;

const double Box::TruncOctMax_ = Box::TRUNCOCTBETA_ + Box::TruncOctDelta_;

/** Used to detect low precision trunc. oct. angles (e.g. 109.47). */
const double Box::TruncOctEps_ = 0.001;

/** Used in IsEq. Larger value than SMALL in case we are reading from
  * low-precision coordinates.
  */
const double Box::EqEps_ = 0.00001;

/** Correspond to CellShapeType */
const char* Box::CellShapeStr_[] = {
  "None",
  "Cubic",
  "Tetragonal",
  "Orthorhombic",
  "Monoclinic",
  "Triclinic",
  "Hexagonal",
  "Rhombohedral",
  "Truncated octahedron",
  "Rhombic dodecahedron"
};

/** Correspond to ParamType */
const char* Box::ParamStr_[] = { "X", "Y", "Z", "alpha", "beta", "gamma" };

// -----------------------------------------------------------------------------
/** \return True if angle is truncated octahedron within a certain range. */
bool Box::IsTruncOct(double angle) {
  return (angle > TruncOctMin_ && angle < TruncOctMax_);
}

/** \return True if the given truncated octahedral angle will cause imaging issues. */
bool Box::BadTruncOctAngle(double angle) {
  return (fabs( TRUNCOCTBETA_ - angle ) > TruncOctEps_);
}

/** \return True if 'lhs' is approximately equal to 'rhs' */
bool Box::IsEq(double lhs, double rhs) {
  return (fabs(rhs - lhs) < EqEps_);
}

/** \return True if cell "A" axis is aligned along the X-axis (i.e. XYZ ABG reference). */
bool Box::Is_X_Aligned() const {
  if (fabs(unitCell_[1]) > Constants::SMALL) return false;
  if (fabs(unitCell_[2]) > Constants::SMALL) return false;
  if (fabs(unitCell_[5]) > Constants::SMALL) return false;
  return true;
}

/** \return True if cell "A" axis is aligned along the X-axis and cell vectors are
  *         orthogonal.
  */
bool Box::Is_X_Aligned_Ortho() const {
  if (fabs(unitCell_[7]) > Constants::SMALL) return false;
  if (fabs(unitCell_[6]) > Constants::SMALL) return false;
  if (fabs(unitCell_[3]) > Constants::SMALL) return false;
  if (fabs(unitCell_[5]) > Constants::SMALL) return false;
  if (fabs(unitCell_[2]) > Constants::SMALL) return false;
  if (fabs(unitCell_[1]) > Constants::SMALL) return false;
  return true;
}

/** \return True if the matrix has symmetric off-diagonal elements. */
bool Box::Is_Symmetric() const {
  if (!IsEq(unitCell_[1], unitCell_[3])) return false;
  if (!IsEq(unitCell_[2], unitCell_[6])) return false;
  if (!IsEq(unitCell_[5], unitCell_[7])) return false;
  return true;
}

/** For debugging purposes, print XYZ ABG and unit/frac cell matrices. */
void Box::printBoxStatus(const char* desc) const {
  mprintf("DEBUG: [%s] Box: %s  is_x_aligned= %i  ortho= %i\n", desc, CellShapeName(), (int)Is_X_Aligned(), (int)Is_X_Aligned_Ortho());
  mprintf("DEBUG:   XYZ= %12.4f %12.4f %8.3f  ABG= %12.4f %12.4f %12.4f\n",
          box_[0], box_[1], box_[2], box_[3], box_[4], box_[5]);
  unitCell_.Print(desc);
  fracCell_.Print("frac");
}

/** For debugging purposes, more in-depth printout. */
void Box::PrintDebug(const char* desc) const {
  mprintf("DEBUG: %s Box: %s\n"
          "DEBUG: %s XYZ= %12.4f %12.4f %12.4f  ABG= %12.4f %12.4f %12.4f\n"
          "DEBUG: %s Ucell = { %20.10E %20.10E %20.10E\n"
          "DEBUG: %s           %20.10E %20.10E %20.10E\n"
          "DEBUG: %s           %20.10E %20.10E %20.10E\n"
          "DEBUG: %s Frac  = { %20.10E %20.10E %20.10E\n"
          "DEBUG: %s           %20.10E %20.10E %20.10E\n"
          "DEBUG: %s           %20.10E %20.10E %20.10E\n"
          "DEBUG: %s is_x_aligned= %i  is_x_aligned_ortho= %i\n",
          desc, CellShapeName(),
          desc, box_[0], box_[1], box_[2], box_[3], box_[4], box_[5],
          desc, unitCell_[0], unitCell_[1], unitCell_[2],
          desc, unitCell_[3], unitCell_[4], unitCell_[5],
          desc, unitCell_[6], unitCell_[7], unitCell_[8],
          desc, fracCell_[0], fracCell_[1], fracCell_[2],
          desc, fracCell_[3], fracCell_[4], fracCell_[5],
          desc, fracCell_[6], fracCell_[7], fracCell_[8],
          desc, (int)Is_X_Aligned(), (int)Is_X_Aligned_Ortho());
}

/** \return Cell shape based on current XYZ (i.e. ABC) alpha beta gamma. */
Box::CellShapeType Box::CellShape() const {
  if (!HasBox()) return NO_SHAPE;
  bool A_equals_B = IsEq( box_[X], box_[Y] );
  bool Lengths_Equal = A_equals_B && IsEq( box_[X], box_[Z] );
  //mprintf("DEBUG: Lengths_Equal= %i  %g  %g  %g (deltas %g and %g)\n", (int)Lengths_Equal, box_[0], box_[1], box_[2], box_[0]-box_[1], box_[0]-box_[2]);
  bool alpha_90 = IsEq( box_[ALPHA], 90.0 );
  bool beta_90  = IsEq( box_[BETA],  90.0 );
  bool gamma_90 = IsEq( box_[GAMMA], 90.0 );

  if (alpha_90 && beta_90 && gamma_90) {
    if (Lengths_Equal)
      return CUBIC;
    else if ( A_equals_B )
      return TETRAGONAL;
    else
      return ORTHORHOMBIC;
  } else if (alpha_90 && gamma_90) {
    return MONOCLINIC;
  } else if (A_equals_B && alpha_90 && beta_90 && IsEq( box_[GAMMA], 120.0 )) {
    return HEXAGONAL;
  } else if (Lengths_Equal) {
    if ( IsTruncOct( box_[ALPHA] ) && IsTruncOct( box_[BETA] ) && IsTruncOct( box_[GAMMA] ) )
      return OCTAHEDRAL;
    else if (beta_90 && IsEq( box_[ALPHA], 60.0 ) && IsEq( box_[GAMMA], 60.0 ))
      return RHOMBIC_DODECAHEDRON;
    else if ( box_[ALPHA] < 120.0 && IsEq( box_[ALPHA], box_[BETA] ) && IsEq( box_[ALPHA], box_[GAMMA] ) )
      return RHOMBOHEDRAL;
  }
  return TRICLINIC;
}

/** Check the box for potential problems. It is expected that if this routine
  * is called, valid box information is present. If not, this is an error.
  * \return 1 if no box, 0 otherwise.
  */
int Box::CheckBox() const {
  // Check for invalid lengths/angles
  bool hasZeros = false;
  for (int i = 0; i < 3; i++) {
    if (box_[i] < Constants::SMALL) {
      mprintf("Warning: Box %s vector length is zero.\n", ParamStr_[i]);
      hasZeros = true;
    }
  }
  for (int i = 3; i < 6; i++) {
    if (box_[i] < Constants::SMALL) {
      mprintf("Warning: Box %s angle is zero.\n", ParamStr_[i]);
      hasZeros = true;
    }
  }
  if (hasZeros) return 1;

  CellShapeType cellShape = CellShape();
  // Check for low-precision truncated octahedron angles.
  if (cellShape == OCTAHEDRAL) {
    if ( BadTruncOctAngle(box_[3]) || BadTruncOctAngle(box_[4]) || BadTruncOctAngle(box_[5]) )
      mprintf("Warning: Low precision truncated octahedron angles detected (%g vs %g).\n"
              "Warning:   If desired, the 'box' command can be used during processing\n"
              "Warning:   to set higher-precision angles.\n", box_[4], TRUNCOCTBETA_);
  }
  // Check for skewed box.
  const double boxFactor = 0.5005;
  double Xaxis_X = unitCell_[0];
  double Yaxis_X = unitCell_[3];
  double Yaxis_Y = unitCell_[4];
  double Zaxis_X = unitCell_[6];
  double Zaxis_Y = unitCell_[7];
  if ( fabs(Yaxis_X) > boxFactor * Xaxis_X ||
       fabs(Zaxis_X) > boxFactor * Xaxis_X ||
       fabs(Zaxis_Y) > boxFactor * Yaxis_Y )
  {
    mprintf("Warning: Box is too skewed to perform accurate imaging.\n"
            "Warning:  Images and imaged distances may not be the absolute minimum.\n");
    // TODO should this return 1?
  }
  return 0;
}

// Box::SetNoBox()
/** Remove all box information. */
void Box::SetNoBox() {
  box_[0] = 0;
  box_[1] = 0;
  box_[2] = 0;
  box_[3] = 0;
  box_[4] = 0;
  box_[5] = 0;
  unitCell_.Zero();
  fracCell_.Zero();
  cellVolume_ = 0;
}

/** Print box info to STDOUT. */
void Box::PrintInfo() const {
  mprintf("\tBox: '%s' XYZ= { %8.3f %8.3f %8.3f } ABG= { %6.2f %6.2f %6.2f }\n",
          CellShapeName(), box_[0], box_[1], box_[2], box_[3], box_[4], box_[5]);
}

// -----------------------------------------------------------------------------
// Static calculation routines.

/** Calculate fractional matrix from unit cell matrix. 
  * \return Volume of the unit cell.
  */
double Box::CalcFracFromUcell(Matrix_3x3& recip, Matrix_3x3 const& ucell) {
  // Get reciprocal vectors
  double u23x = ucell[4]*ucell[8] - ucell[5]*ucell[7];
  double u23y = ucell[5]*ucell[6] - ucell[3]*ucell[8];
  double u23z = ucell[3]*ucell[7] - ucell[4]*ucell[6];
  double u31x = ucell[7]*ucell[2] - ucell[8]*ucell[1];
  double u31y = ucell[8]*ucell[0] - ucell[6]*ucell[2];
  double u31z = ucell[6]*ucell[1] - ucell[7]*ucell[0];
  double u12x = ucell[1]*ucell[5] - ucell[2]*ucell[4];
  double u12y = ucell[2]*ucell[3] - ucell[0]*ucell[5];
  double u12z = ucell[0]*ucell[4] - ucell[1]*ucell[3];
  double volume = ucell[0]*u23x + ucell[1]*u23y + ucell[2]*u23z;
  double onevolume = 1.0 / volume;

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

/** Calculate unit cell matrix from XYZ ABG array. 
  * NOTE: Since cos(90 deg.) is not numerically stable (i.e. it is never
  *       exactly 0.0), trap 90 deg. angles explicitly to ensure
  *       the unit cell vectors are properly orthogonal.
  */
void Box::CalcUcellFromXyzAbg(Matrix_3x3& ucell, const double* xyzabg) {
  // A vector
  ucell[0] = xyzabg[0]; // u(1,1)
  ucell[1] = 0.0;       // u(2,1)
  ucell[2] = 0.0;       // u(3,1)
  // B vector
  if ( fabs(xyzabg[5] - 90.0) < Constants::SMALL ) {
    ucell[3] = 0.0;
    ucell[4] = xyzabg[1];
  } else {
    ucell[3] = xyzabg[1]*cos(Constants::DEGRAD*xyzabg[5]); // u(1,2)
    ucell[4] = xyzabg[1]*sin(Constants::DEGRAD*xyzabg[5]); // u(2,2)
  }
  ucell[5] = 0.0;                                          // u(3,2)
  // C vector
  if ( fabs(xyzabg[4] - 90.0) < Constants::SMALL )
    ucell[6] = 0.0;
  else
    ucell[6] = xyzabg[2]*cos(Constants::DEGRAD*xyzabg[4]);
  double y_z_cos_alpha;
  if ( fabs(xyzabg[3] - 90.0) < Constants::SMALL )
    y_z_cos_alpha = 0.0;
  else
    y_z_cos_alpha = xyzabg[1]*xyzabg[2]*cos(Constants::DEGRAD*xyzabg[3]);
  ucell[7] = (y_z_cos_alpha - ucell[6]*ucell[3]) / ucell[4];
  ucell[8] = sqrt(xyzabg[2]*xyzabg[2] - ucell[6]*ucell[6] - ucell[7]*ucell[7]);
/*
    //The old code for reference.
    ucell[0] = xyzabg[0]; // u(1,1)
    ucell[1] = 0.0;     // u(2,1)
    ucell[2] = 0.0;     // u(3,1)
    ucell[3] = xyzabg[1]*cos(Constants::DEGRAD*xyzabg[5]); // u(1,2)
    ucell[4] = xyzabg[1]*sin(Constants::DEGRAD*xyzabg[5]); // u(2,2)
    ucell[5] = 0.0;                                    // u(3,2)
    ucell[6] = xyzabg[2]*cos(Constants::DEGRAD*xyzabg[4]);
    ucell[7] = (xyzabg[1]*xyzabg[2]*cos(Constants::DEGRAD*xyzabg[3]) - ucell[6]*ucell[3]) / ucell[4];
    ucell[8] = sqrt(xyzabg[2]*xyzabg[2] - ucell[6]*ucell[6] - ucell[7]*ucell[7]);
*/
}

/** Calculate XYZ ABG array from unit cell matrix */
void Box::CalcXyzAbgFromUcell(double* box, Matrix_3x3 const& ucell) {
  Vec3 x_axis = ucell.Row1();
  Vec3 y_axis = ucell.Row2();
  Vec3 z_axis = ucell.Row3();
  box[0] = x_axis.Normalize(); // A
  box[1] = y_axis.Normalize(); // B
  box[2] = z_axis.Normalize(); // C
  box[3] = y_axis.Angle( z_axis ) * Constants::RADDEG; // alpha
  box[4] = x_axis.Angle( z_axis ) * Constants::RADDEG; // beta
  box[5] = x_axis.Angle( y_axis ) * Constants::RADDEG; // gamma
}

/** Convert symmetric shape matrix data to params x, y, z, a, b, g */
void Box::CalcXyzAbgFromShape(double* box, const double* boxtmp)
{
  double boxX = sqrt( boxtmp[0]*boxtmp[0] + boxtmp[1]*boxtmp[1] + boxtmp[3]*boxtmp[3] );
  double boxY = sqrt( boxtmp[1]*boxtmp[1] + boxtmp[2]*boxtmp[2] + boxtmp[4]*boxtmp[4] );
  double boxZ = sqrt( boxtmp[3]*boxtmp[3] + boxtmp[4]*boxtmp[4] + boxtmp[5]*boxtmp[5] );
  double boxXY = boxtmp[1]*(boxtmp[0] + boxtmp[2]) + boxtmp[3]*boxtmp[4];
  double boxYZ = boxtmp[4]*(boxtmp[2] + boxtmp[5]) + boxtmp[1]*boxtmp[3];
  double boxXZ = boxtmp[3]*(boxtmp[0] + boxtmp[5]) + boxtmp[1]*boxtmp[4];
  double alpha = acos( boxYZ / (boxY*boxZ) ) * Constants::RADDEG;
  double beta  = acos( boxXZ / (boxX*boxZ) ) * Constants::RADDEG;
  double gamma = acos( boxXY / (boxX*boxY) ) * Constants::RADDEG;
  //mprintf("DEBUG: Box XYZ= %g %g %g  ABG= %g %g %g\n", boxX, boxY, boxZ, alpha, beta, gamma);
  box[0] = boxX;
  box[1] = boxY;
  box[2] = boxZ;
  box[3] = alpha;
  box[4] = beta;
  box[5] = gamma;
}

/** Convert unit cell parameters (X, Y, Z, a, b, g) to symmetric shape matrix
  * (S11, S12, S22, S13, S23, S33).
  */
void Box::CalcShapeFromXyzAbg(double* shape, const double* box)
{
  // Calculate metric tensor HtH:
  //   HtH(i,j) = vi * vj
  // where vx are basis vectors i and j. Given that v0 is a, v1 is b, v2 is c:
  //       a^2 a*b a*c
  // HtH = b*a b^2 b*c
  //       c*a c*b c^2
  Matrix_3x3 HtH;

  HtH[0] = box[0] * box[0];
  HtH[4] = box[1] * box[1];
  HtH[8] = box[2] * box[2];

  // Angles near 90 have elements set to 0.0.
  // XY (gamma)
  if (fabs(box[5] - 90.0) > Constants::SMALL)
    HtH[3] = box[0]*box[1]*cos(Constants::DEGRAD*box[5]);
  else
    HtH[3] = 0.0;
  HtH[1] = HtH[3];
  // XZ (beta)
  if (fabs(box[4] - 90.0) > Constants::SMALL)
    HtH[6] = box[0]*box[2]*cos(Constants::DEGRAD*box[4]);
  else
    HtH[6] = 0.0;
  HtH[2] = HtH[6];
  // YZ (alpha)
  if (fabs(box[3] - 90.0) > Constants::SMALL)
    HtH[7] = box[1]*box[2]*cos(Constants::DEGRAD*box[3]);
  else
    HtH[7] = 0.0;
  HtH[5] = HtH[7];

  // Diagonalize HtH
  //HtH.Print("HtH"); // DEBUG
  Vec3 Evals;
  if (HtH.Diagonalize( Evals )) {
    mprinterr("Error: Could not diagonalize metric tensor.\n");
    for (int i=0; i<6; i++) shape[i] = 0.0;
    return;
  }

  if (Evals[0] < Constants::SMALL ||
      Evals[1] < Constants::SMALL ||
      Evals[2] < Constants::SMALL)
  {
    mprinterr("Error: Obtained negative eigenvalues when attempting to"
              " diagonalize metric tensor.\n");
    return;
  }
  //Evals.Print("Cvals"); // DEBUG
  //HtH.Print("Cpptraj"); // DEBUG

  double A = sqrt( Evals[0] );
  double B = sqrt( Evals[1] );
  double C = sqrt( Evals[2] );

  shape[0] = A*HtH[0]*HtH[0] + B*HtH[1]*HtH[1] + C*HtH[2]*HtH[2];
  shape[2] = A*HtH[3]*HtH[3] + B*HtH[4]*HtH[4] + C*HtH[5]*HtH[5];
  shape[5] = A*HtH[6]*HtH[6] + B*HtH[7]*HtH[7] + C*HtH[8]*HtH[8];
  shape[1] = A*HtH[0]*HtH[3] + B*HtH[1]*HtH[4] + C*HtH[2]*HtH[5];
  shape[3] = A*HtH[0]*HtH[6] + B*HtH[1]*HtH[7] + C*HtH[2]*HtH[8];
  shape[4] = A*HtH[3]*HtH[6] + B*HtH[4]*HtH[7] + C*HtH[5]*HtH[8];
}

// -----------------------------------------------------------------------------
// Setup routines

/** Set unit cell, fractional cell, and XYZ ABG array from shape matrix. */
int Box::SetupFromShapeMatrix(const double* shape) {
  unitCell_[0] = shape[0];
  unitCell_[1] = shape[1];
  unitCell_[2] = shape[3];

  unitCell_[3] = shape[1];
  unitCell_[4] = shape[2];
  unitCell_[5] = shape[4];

  unitCell_[6] = shape[3];
  unitCell_[7] = shape[4];
  unitCell_[8] = shape[5];

  CalcXyzAbgFromShape(box_, shape);

  cellVolume_ = CalcFracFromUcell(fracCell_, unitCell_);

# ifdef DEBUG_BOX
  printBoxStatus("SetupFromShapeMatrix");
# endif
  if (CheckBox()) {
    SetNoBox();
    return 1;
  }
  return 0;
}

/** Set up Xyz Abg array and frac cell from unit cell. */
int Box::SetupFromUcell(const double* ucell) {
  std::copy(ucell, ucell+9, unitCell_.Dptr());

  CalcXyzAbgFromUcell(box_, unitCell_);

  cellVolume_ = CalcFracFromUcell(fracCell_, unitCell_);

# ifdef DEBUG_BOX
  printBoxStatus("SetupFromUcell");
# endif
  if (CheckBox()) {
    SetNoBox();
    return 1;
  }
  return 0;
}

/** Set unit cell and fractional cell from XYZ ABG parameters. */
int Box::SetupFromXyzAbg(double bx, double by, double bz, double ba, double bb, double bg) {
  box_[0] = bx;
  box_[1] = by;
  box_[2] = bz;
  box_[3] = ba;
  box_[4] = bb;
  box_[5] = bg;

  CalcUcellFromXyzAbg(unitCell_, box_);

  cellVolume_ = CalcFracFromUcell(fracCell_, unitCell_);

# ifdef DEBUG_BOX
  printBoxStatus("SetupFromXyzAbgIndividual");
# endif
  if (CheckBox()) {
    SetNoBox();
    return 1;
  }
  return 0;
}

/** Set unit cell and fractional cell from XYZ ABG array. */
int Box::SetupFromXyzAbg(const double* xyzabg) {
  box_[0] = xyzabg[0];
  box_[1] = xyzabg[1];
  box_[2] = xyzabg[2];
  box_[3] = xyzabg[3];
  box_[4] = xyzabg[4];
  box_[5] = xyzabg[5];

  CalcUcellFromXyzAbg(unitCell_, xyzabg);

  cellVolume_ = CalcFracFromUcell(fracCell_, unitCell_);

# ifdef DEBUG_BOX
  printBoxStatus("SetupFromXyzAbg");
# endif
  if (CheckBox()) {
    SetNoBox();
    return 1;
  }
  return 0;
}

// -----------------------------------------------------------------------------
// Assign routines

/** Assign Xyz Abg array and frac cell from unit cell. */
void Box::AssignFromUcell(const double* ucell) {
  // Sanity check
  //if (btype_ == NOBOX) {
  //  mprintf("Internal Error: AssignFromUcell(): No box has been set.\n");
  //  return;
  //}

  std::copy(ucell, ucell+9, unitCell_.Dptr());

  CalcXyzAbgFromUcell(box_, unitCell_);

  cellVolume_ = CalcFracFromUcell(fracCell_, unitCell_);
# ifdef DEBUG_BOX
  printBoxStatus("AssignFromUcell");
# endif
}

/** Assign from XYZ ABG parameters. */
void Box::AssignFromXyzAbg(double bx, double by, double bz, double ba, double bb, double bg) {
  // Sanity check
  //if (btype_ == NOBOX) {
  //  mprintf("Internal Error: AssignFromXyzAbgIndividual(): No box has been set.\n");
  //  return;
  //}

  box_[0] = bx;
  box_[1] = by;
  box_[2] = bz;
  box_[3] = ba;
  box_[4] = bb;
  box_[5] = bg;

  CalcUcellFromXyzAbg(unitCell_, box_);

  cellVolume_ = CalcFracFromUcell(fracCell_, unitCell_);
# ifdef DEBUG_BOX
  printBoxStatus("AssignFromXyzAbgIndividual");
# endif
}

/** Assign from XYZ ABG array. */
void Box::AssignFromXyzAbg(const double* xyzabg) {
  // Sanity check
  //if (btype_ == NOBOX) {
  //  mprintf("Internal Error: AssignFromXyzAbg(): No box has been set.\n");
  //  return;
  //}
  // TODO detect orthogonal?
  box_[0] = xyzabg[0];
  box_[1] = xyzabg[1];
  box_[2] = xyzabg[2];
  box_[3] = xyzabg[3];
  box_[4] = xyzabg[4];
  box_[5] = xyzabg[5];

  CalcUcellFromXyzAbg(unitCell_, xyzabg);
  //CalcUcellFromXyzAbg(unitCell_, btype_, box_, 1.0);

  cellVolume_ = CalcFracFromUcell(fracCell_, unitCell_);
# ifdef DEBUG_BOX
  printBoxStatus("AssignFromXyzAbg");
# endif
}

/** Assign from symmetric shape matrix. */
void Box::AssignFromShapeMatrix(const double* shape) {
  // Sanity check
  //if (btype_ == NOBOX) {
  //  mprintf("Internal Error: AssignFromShapeMatrix(): No box has been set.\n");
  //  return;
  //}

  unitCell_[0] = shape[0];
  unitCell_[1] = shape[1];
  unitCell_[2] = shape[3];

  unitCell_[3] = shape[1];
  unitCell_[4] = shape[2];
  unitCell_[5] = shape[4];

  unitCell_[6] = shape[3];
  unitCell_[7] = shape[4];
  unitCell_[8] = shape[5];

  CalcXyzAbgFromShape(box_, shape);

  cellVolume_ = CalcFracFromUcell(fracCell_, unitCell_);
# ifdef DEBUG_BOX
  printBoxStatus("AssignFromShapeMatrix");
# endif
}

// -----------------------------------------------------------------------------
/** Rotate unit cell vectors; recalculate fractional cell vectors. */
void Box::RotateUcell(Matrix_3x3 const& rot) {
  if (!HasBox()) return;
  // Rotate each unit cell vector
  Vec3 ucellx = rot * unitCell_.Row1();
  Vec3 ucelly = rot * unitCell_.Row2();
  Vec3 ucellz = rot * unitCell_.Row3();
  //ucellx.Print("ucellx");
  //ucelly.Print("ucelly");
  //ucellz.Print("ucellz");
  // Set new unit cell matrix
  unitCell_ = Matrix_3x3( ucellx[0], ucellx[1], ucellx[2],
                          ucelly[0], ucelly[1], ucelly[2],
                          ucellz[0], ucellz[1], ucellz[2] );
  //unitCell_.Print("RotatedUcell");
  // Recalculate the fractional cell; volume is unchanged.
  CalcFracFromUcell(fracCell_, unitCell_);
}

/** Inverse rotation of unit cell vectors; recalculate fractional cell vectors. */
void Box::InverseRotateUcell(Matrix_3x3 const& rot) {
  if (!HasBox()) return;
  // Inverse rotate each unit cell vector
  Vec3 ucellx = rot.TransposeMult( unitCell_.Row1() );
  Vec3 ucelly = rot.TransposeMult( unitCell_.Row2() );
  Vec3 ucellz = rot.TransposeMult( unitCell_.Row3() );
  // Set new unit cell matrix
  unitCell_ = Matrix_3x3( ucellx[0], ucellx[1], ucellx[2],
                          ucelly[0], ucelly[1], ucelly[2],
                          ucellz[0], ucellz[1], ucellz[2] );
  // Recalculate the fractional cell; volume is unchanged.
  CalcFracFromUcell(fracCell_, unitCell_);
}

/** Set 'shape' with values for symmetric shape matrix from XYZ ABG array. */
void Box::GetSymmetricShapeMatrix(double* shape) const {
  CalcShapeFromXyzAbg(shape, box_);
}

//  Box::RecipLengths()
/** \return Vector containing reciprocal lengths from fractional cell matrix.
  * Used primarily by the Ewald and PairList routines.
  */
Vec3 Box::RecipLengths() const {
  return Vec3( 1.0/sqrt(fracCell_[0]*fracCell_[0] + fracCell_[1]*fracCell_[1] + fracCell_[2]*fracCell_[2]),
               1.0/sqrt(fracCell_[3]*fracCell_[3] + fracCell_[4]*fracCell_[4] + fracCell_[5]*fracCell_[5]),
               1.0/sqrt(fracCell_[6]*fracCell_[6] + fracCell_[7]*fracCell_[7] + fracCell_[8]*fracCell_[8]) );
}
