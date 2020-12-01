#include <cmath> // cos, sin, sqrt, fabs
#include "Box.h"
#include "Constants.h" // DEGRAD
#include "CpptrajStdio.h"
#include <algorithm> // std::copy

// CONSTRUCTOR
Box::Box() :
  btype_(NOBOX),
  cellVolume_(0) //, debug_(0)
{
  box_[0] = 0;
  box_[1] = 0;
  box_[2] = 0;
  box_[3] = 0;
  box_[4] = 0;
  box_[5] = 0;
}

/*
Box::Box(const double* bIn) //: debug_(0)
{
  SetBox( bIn );
}

Box::Box(const float* bIn) { SetBox( bIn ); }
*/

/** CONSTRUCTOR - Set up from unit cell matrix. */
Box::Box(Matrix_3x3 const& ucell) { 
  SetupFromUcell( ucell.Dptr() );
}

/** COPY CONSTRUCTOR */
Box::Box(const Box& rhs) :
  btype_(rhs.btype_),
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
  //debug_ = rhs.debug_;
  btype_ = rhs.btype_;
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

/// Swap box type
static inline void bswap(Box::BoxType& b1, Box::BoxType& b2) {
  Box::BoxType btemp = b1;
  b1 = b2;
  b2 = btemp;
}

/** Swap this box with given box. */
void Box::swap(Box& rhs) {
  bswap( btype_,  rhs.btype_ );
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
int Box::SyncBox(Parallel::Comm const& commIn) {
  commIn.MasterBcast( &btype_, 1, MPI_INT );
  commIn.MasterBcast( box_,    6, MPI_DOUBLE );
  unitCell_.SyncMatrix( commIn );
  fracCell_.SyncMatrix( commIn );
  commIn.MasterBcast( &cellVolume_, 1, MPI_DOUBLE );
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

const char* Box::BoxNames_[] = {
  "None", "Orthogonal", "Trunc. Oct.", "Rhombic Dodec.", "Non-orthogonal"
};

// -----------------------------------------------------------------------------
bool Box::IsTruncOct(double angle) {
  return (angle > TruncOctMin_ && angle < TruncOctMax_);
}

bool Box::BadTruncOctAngle(double angle) {
  return (fabs( TRUNCOCTBETA_ - angle ) > TruncOctEps_);
}

/// \return True if 'angle' is approximately equal to 'tgt'
bool Box::IsEq(double angle, double tgt) {
  return (fabs(tgt - angle) < Constants::SMALL);
}

/** \return True if cell is aligned along "normal" (i.e. XYZ ABG) reference. */
bool Box::IsNormal() const {
  if (fabs(unitCell_[1]) > 0.0) return false;
  if (fabs(unitCell_[2]) > 0.0) return false;
  if (fabs(unitCell_[5]) > 0.0) return false;
  return true;
}

/** \return True if cell is aligned along "normal" and is orthogonal. */
bool Box::IsOrthoNormal() const {
  if (fabs(unitCell_[7]) > 0.0) return false;
  if (fabs(unitCell_[6]) > 0.0) return false;
  if (fabs(unitCell_[5]) > 0.0) return false;
  if (fabs(unitCell_[3]) > 0.0) return false;
  if (fabs(unitCell_[2]) > 0.0) return false;
  if (fabs(unitCell_[1]) > 0.0) return false;
  return true;
}

// Box::SetBoxType()
/** Determine box type (none/ortho/nonortho) based on box angles. */
void Box::SetBoxType() {
  btype_ = NONORTHO;
  // TODO error if any of box is zero?
  bool noLengths = (box_[0] < Constants::SMALL &&
                    box_[1] < Constants::SMALL &&
                    box_[2] < Constants::SMALL);
  bool noAngles = ( box_[3] <= 0 && box_[4] <= 0 && box_[5] <= 0);
  if ( noLengths ) {
    // No lengths, no box
    btype_ = NOBOX;
    if (!noAngles)
      mprintf("Warning: Box length(s) <= 0.0; setting box to NONE.\n");
  } else if ( noAngles ) {
    // No angles, no box
    mprintf("Warning: Box angle(s) <= 0.0; setting box to NONE.\n");
    btype_ = NOBOX;
  } else if (box_[3] == 90.0 && box_[4] == 90.0 && box_[5] == 90.0)
    // All 90, orthogonal
    btype_ = ORTHO;
  else if ( IsTruncOct( box_[3] ) && IsTruncOct( box_[4] ) && IsTruncOct( box_[5] ) )
    // All 109.47, truncated octahedron
    btype_ = TRUNCOCT;
  else if ( IsEq(box_[3],60.0) && IsEq(box_[4],90.0) && IsEq(box_[5],60.0) )
    // 60/90/60, rhombic dodecahedron
    btype_ = RHOMBIC;
  else
    btype_ = NONORTHO;
  //if (debug_>0) mprintf("\tBox type is %s (beta=%lf)\n",TypeName(), box_[4]);

  if (btype_ == TRUNCOCT) {
    // Check for low-precision truncated octahedron angles.
    if ( BadTruncOctAngle(box_[3]) || BadTruncOctAngle(box_[4]) || BadTruncOctAngle(box_[5]) )
      mprintf("Warning: Low precision truncated octahedron angles detected (%g vs %g).\n"
              "Warning:   If desired, the 'box' command can be used during processing\n"
              "Warning:   to set higher-precision angles.\n", box_[4], TRUNCOCTBETA_);
  } else if (btype_ == NONORTHO) {
    // Check for skewed box.
    const double boxFactor = 0.5005;
    double Xaxis_X = box_[0];
    double Yaxis_X = box_[1]*cos(Constants::DEGRAD*box_[5]);
    double Yaxis_Y = box_[1]*sin(Constants::DEGRAD*box_[5]);
    double Zaxis_X = box_[2]*cos(Constants::DEGRAD*box_[4]);
    double Zaxis_Y = (box_[1]*box_[2]*cos(Constants::DEGRAD*box_[3]) - Zaxis_X*Yaxis_X) / Yaxis_Y;
    if ( fabs(Yaxis_X) > boxFactor * Xaxis_X ||
         fabs(Zaxis_X) > boxFactor * Xaxis_X ||
         fabs(Zaxis_Y) > boxFactor * Yaxis_Y )
    {
      mprintf("Warning: Non-orthogonal box is too skewed to perform accurate imaging.\n"
              "Warning:  Images and imaged distances may not be the absolute minimum.\n");
    }
  }
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
  btype_ = NOBOX;
  unitCell_.Zero();
  fracCell_.Zero();
  cellVolume_ = 0;
}

//  Box::RecipLengths()
/** \return Vector containing reciprocal lengths from given fractional cell matrix.
  * Used primarily by the Ewald and PairList routines.
  */
Vec3 Box::RecipLengths(Matrix_3x3 const& recip) {
  return Vec3( 1.0/sqrt(recip[0]*recip[0] + recip[1]*recip[1] + recip[2]*recip[2]),
               1.0/sqrt(recip[3]*recip[3] + recip[4]*recip[4] + recip[5]*recip[5]),
               1.0/sqrt(recip[6]*recip[6] + recip[7]*recip[7] + recip[8]*recip[8]) );
}

/** Print box info to STDOUT. */
void Box::PrintInfo() const {
  mprintf("\tBox: '%s' XYZ= { %8.3f %8.3f %8.3f } ABG= { %6.2f %6.2f %6.2f }\n",
          BoxNames_[btype_], box_[0], box_[1], box_[2], box_[3], box_[4], box_[5]);
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

/** Calculate unit cell matrix from XYZ ABG array. */
void Box::CalcUcellFromXyzAbg(Matrix_3x3& ucell, const double* xyzabg) {
  // If box lengths are zero no imaging possible
  //if (xyzabg[0]==0.0 || xyzabg[1]==0.0 || xyzabg[2]==0.0) {
  //  ucell.Zero();
  //  recip.Zero();
  //  return -1.0;
  //}
  ucell[0] = xyzabg[0]; // u(1,1)
  ucell[1] = 0.0;     // u(2,1)
  ucell[2] = 0.0;     // u(3,1)
  ucell[3] = xyzabg[1]*cos(Constants::DEGRAD*xyzabg[5]); // u(1,2)
  ucell[4] = xyzabg[1]*sin(Constants::DEGRAD*xyzabg[5]); // u(2,2)
  ucell[5] = 0.0;                                    // u(3,2)
  ucell[6] = xyzabg[2]*cos(Constants::DEGRAD*xyzabg[4]);
  ucell[7] = (xyzabg[1]*xyzabg[2]*cos(Constants::DEGRAD*xyzabg[3]) - ucell[6]*ucell[3]) / ucell[4];
  ucell[8] = sqrt(xyzabg[2]*xyzabg[2] - ucell[6]*ucell[6] - ucell[7]*ucell[7]);
}

// Box::CalcUcellFromXyzAbg()
void Box::CalcUcellFromXyzAbg(Matrix_3x3& ucell, BoxType btype, const double* xyzabg, double scale) {
  double by, bz;
  //switch (btype) {
  //  case NOBOX: ucell.Zero(); break;
  //  case ORTHO:
  if (btype == ORTHO) {
      ucell[0] = xyzabg[0] * scale;
      ucell[1] = 0.0;
      ucell[2] = 0.0;
      ucell[3] = 0.0;
      ucell[4] = xyzabg[1] * scale;
      ucell[5] = 0.0;
      ucell[6] = 0.0;
      ucell[7] = 0.0;
      ucell[8] = xyzabg[2] * scale;
  //    break;
  //  case TRUNCOCT:
  //  case RHOMBIC:
  //  case NONORTHO:
  } else {
      by = xyzabg[1] * scale;
      bz = xyzabg[2] * scale;
      ucell[0] = xyzabg[0] * scale;
      ucell[1] = 0.0;
      ucell[2] = 0.0;
      ucell[3] = by*cos(Constants::DEGRAD*xyzabg[5]);
      ucell[4] = by*sin(Constants::DEGRAD*xyzabg[5]);
      ucell[5] = 0.0;
      ucell[6] = bz*cos(Constants::DEGRAD*xyzabg[4]);
      ucell[7] = (by*bz*cos(Constants::DEGRAD*xyzabg[3]) - ucell[6]*ucell[3]) / ucell[4];
      ucell[8] = sqrt(bz*bz - ucell[6]*ucell[6] - ucell[7]*ucell[7]);
  //    break;
  }
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
void CalcShapeFromXyzAbg(double* shape, const double* box)
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
void Box::SetupFromShapeMatrix(const double* shape) {
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

  SetBoxType();
}

/** Set up Xyz Abg array and frac cell from unit cell. */
void Box::SetupFromUcell(const double* ucell) {
  std::copy(ucell, ucell+9, unitCell_.Dptr());

  CalcXyzAbgFromUcell(box_, unitCell_);

  cellVolume_ = CalcFracFromUcell(fracCell_, unitCell_);

  SetBoxType();
}

/** Set unit cell and fractional cell from XYZ ABG array. */
void Box::SetupFromXyzAbg(const double* xyzabg) {
  box_[0] = xyzabg[0];
  box_[1] = xyzabg[1];
  box_[2] = xyzabg[2];
  box_[3] = xyzabg[3];
  box_[4] = xyzabg[4];
  box_[5] = xyzabg[5];

  CalcUcellFromXyzAbg(unitCell_, xyzabg);

  cellVolume_ = CalcFracFromUcell(fracCell_, unitCell_);

  SetBoxType();
}

// -----------------------------------------------------------------------------
// Assign routines

/** Assign from XYZ ABG array. */
void Box::AssignFromXyzAbg(const double* xyzabg) {
  // Sanity check
  if (btype_ == NOBOX) {
    mprinterr("Internal Error: AssignFromXyzAbg(): No box has been set.\n");
    return;
  }
  box_[0] = xyzabg[0];
  box_[1] = xyzabg[1];
  box_[2] = xyzabg[2];
  box_[3] = xyzabg[3];
  box_[4] = xyzabg[4];
  box_[5] = xyzabg[5];

  //CalcUcellFromXyzAbg(unitCell_, xyzabg);
  CalcUcellFromXyzAbg(unitCell_, btype_, box_, 1.0);

  cellVolume_ = CalcFracFromUcell(fracCell_, unitCell_);
}

/** Assign from symmetric shape matrix. */
void Box::AssignFromShapeMatrix(const double* shape) {
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
}

// -----------------------------------------------------------------------------
/** Set 'shape' with values for symmetric shape matrix from XYZ ABG array. */
void Box::GetSymmetricShapeMatrix(double* shape) const {
  CalcShapeFromXyzAbg(shape, box_);
}

/*
// Box::SetBetaLengths()
void Box::SetBetaLengths(double beta, double xin, double yin, double zin) {
  box_[0] = xin;
  box_[1] = yin;
  box_[2] = zin;
  box_[3] = 0;
  box_[4] = beta;
  box_[5] = 0;
  SetBoxType();
}*/

/** Set box from double[6] array */
/*
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
}*/

/** Set box from float[6] array */
/*
void Box::SetBox(const float* xyzabg) {
  if (xyzabg == 0) {
    mprinterr("Error: Box input array is null\n");
    return;
  }
  box_[0] = (double)xyzabg[0];
  box_[1] = (double)xyzabg[1];
  box_[2] = (double)xyzabg[2];
  box_[3] = (double)xyzabg[3];
  box_[4] = (double)xyzabg[4];
  box_[5] = (double)xyzabg[5];
  SetBoxType();
}
*/


/*
void Box::SetBox(float A, float B, float C, float alpha, float beta, float gamma)
{
  box_[0] = A;
  box_[1] = B;
  box_[2] = C;
  box_[3] = alpha;
  box_[4] = beta;
  box_[5] = gamma;
  SetBoxType();
}*/

// Box::SetTruncOct()
/** Set as truncated octahedron with all lengths set to whatever X is. */
/*
void Box::SetTruncOct() {
  box_[1] = box_[0];
  box_[2] = box_[0];
  box_[3] = TRUNCOCTBETA_;
  box_[4] = TRUNCOCTBETA_;
  box_[5] = TRUNCOCTBETA_;
  btype_ = TRUNCOCT;
  mprintf("Info: Setting box to be perfect truncated octahedron (a=b=g=%g)\n", box_[3]);
}*/

// Box::SetMissingInfo()
/** Set this box info from rhs if <= 0. */
/*
void Box::SetMissingInfo(const Box& rhs) {
  if (box_[0] <= 0) box_[0] = rhs.box_[0];
  if (box_[1] <= 0) box_[1] = rhs.box_[1];
  if (box_[2] <= 0) box_[2] = rhs.box_[2];
  if (box_[3] <= 0) box_[3] = rhs.box_[3];
  if (box_[4] <= 0) box_[4] = rhs.box_[4];
  if (box_[5] <= 0) box_[5] = rhs.box_[5];
  SetBoxType();
}*/

// Box::ToRecip()
/** Use box coordinates to calculate unit cell and fractional matrix for use
  * with imaging routines. Return cell volume.
  */ // FIXME this should eventually be removed
double Box::ToRecip(Matrix_3x3& ucell, Matrix_3x3& recip) const {
  ucell = unitCell_;
  recip = fracCell_;
  return cellVolume_;
/*
  double u12x,u12y,u12z;
  double u23x,u23y,u23z;
  double u31x,u31y,u31z;
  double volume,onevolume;
  // If box lengths are zero no imaging possible
  if (box_[0]==0.0 || box_[1]==0.0 || box_[2]==0.0) {
    ucell.Zero();
    recip.Zero();
    return -1.0;
  }
  ucell[0] = box_[0]; // u(1,1)
  ucell[1] = 0.0;     // u(2,1)
  ucell[2] = 0.0;     // u(3,1)
  ucell[3] = box_[1]*cos(Constants::DEGRAD*box_[5]); // u(1,2)
  ucell[4] = box_[1]*sin(Constants::DEGRAD*box_[5]); // u(2,2)
  ucell[5] = 0.0;                                    // u(3,2)
  ucell[6] = box_[2]*cos(Constants::DEGRAD*box_[4]);
  ucell[7] = (box_[1]*box_[2]*cos(Constants::DEGRAD*box_[3]) - ucell[6]*ucell[3]) / ucell[4];
  ucell[8] = sqrt(box_[2]*box_[2] - ucell[6]*ucell[6] - ucell[7]*ucell[7]);

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

  return volume;*/
}
