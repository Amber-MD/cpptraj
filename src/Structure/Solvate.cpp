#include "Solvate.h"
#include "../ArgList.h"
#include "../CpptrajStdio.h"
#include "../DataSet_Coords.h"
#include "../DataSetList.h"
#include "../Frame.h"
#include "../Topology.h"
#include "../Parm/ParameterSet.h"
#include <algorithm> //std::max
#include <cmath> // cos, sin, sqrt

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
Solvate::Solvate() :
  debug_(0),
  bufferX_(0),
  bufferY_(0),
  bufferZ_(0),
  bufferD_(0),
  closeness_(1.0),
  isotropic_(false),
  clip_(true),
  center_(true)
{
}

const char* Solvate::SetboxKeywords() {
  return "[{buffer <buffer> | bufx <bufx> bufy <bufy> bufz <bufz>}]";
}

const char* Solvate::SolvateKeywords1() {
  return "{buffer <buffer> | bufx <bufx> bufy <bufy> bufz <bufz>}";
}

const char* Solvate::SolvateKeywords2() {
  return "solventbox <unit> [closeness <closeness>] [iso] [nocenter]";
}

/** Get any buffer arguments */
int Solvate::getBufferArg(ArgList& argIn, double defaultBuffer) {
  if (argIn.Contains("buffer")) {
    bufferX_ = argIn.getKeyDouble("buffer", defaultBuffer);
    bufferY_ = bufferX_;
    bufferZ_ = bufferX_;
    bufferD_ = bufferX_;
  } else {
    bufferX_ = argIn.getKeyDouble("bufx", defaultBuffer);
    bufferY_ = argIn.getKeyDouble("bufy", defaultBuffer);
    bufferZ_ = argIn.getKeyDouble("bufz", defaultBuffer);
    bufferD_ = argIn.getKeyDouble("bufd", defaultBuffer);
  }
  if (bufferX_ < 0 || bufferY_ < 0 || bufferZ_ < 0) {
    mprinterr("Error: Either 'buffer' or 'bufx/bufy/bufx' must be specified and >= 0\n");
    return 1;
  }
  return 0;
}

/** Initialize arguments. */
int Solvate::InitSolvate(ArgList& argIn, bool octIn, int debugIn) {
  debug_ = debugIn;
  doTruncatedOct_ = octIn;

  if (getBufferArg(argIn, -1.0)) return 1;

  isotropic_ = argIn.hasKey("iso");

  solventBoxName_ = argIn.GetStringKey("solventbox");
  if (solventBoxName_.empty()) {
    mprinterr("Error: Specify solvent box unit name with 'solventbox'\n");
    return 1;
  }

  closeness_ = argIn.getKeyDouble("closeness", 1.0);
  center_ = !argIn.hasKey("nocenter");

  if (doTruncatedOct_) {
    if (!clip_) {
      mprinterr("Error: Truncated octahedral box currently requires 'clip'.\n");
      return 1;
    }
  }

  return 0;
}

/** Initialize args for setbox */
int Solvate::InitSetbox(ArgList& argIn, int debugIn) {
  debug_ = debugIn;

  if (getBufferArg(argIn, 0.0)) return 1;

  return 0;
}

/** Print info to stdout */
void Solvate::PrintSolvateInfo() const {
  mprintf("\tSolvent buffer XYZ: %g %g %g Ang.\n", bufferX_, bufferY_, bufferZ_);
  mprintf("\tSolvent closeness: %g Ang.\n", closeness_);
  mprintf("\tSolvent isotropic=%i  clip=%i  center_=%i  oct=%i\n", (int)isotropic_, (int)clip_, (int)center_, (int)doTruncatedOct_);
} 

/** Get solvent unit box from DataSetList */
DataSet_Coords* Solvate::GetSolventUnit(DataSetList const& DSL) const {
  if (solventBoxName_.empty()) {
    mprinterr("Internal Error: Solvate::GetSolventUnit() called before solventBoxName_ set.\n");
    return 0;
  }
  DataSetList sets = DSL.SelectGroupSets( "*", DataSet::COORDINATES ); // TODO specific set type for units?
  // First try to match aspect, then match name
  DataSet_Coords* solventUnit = 0;
  // Aspect
  for (DataSetList::const_iterator it = sets.begin(); it != sets.end(); ++it)
  {
    DataSet_Coords* ds = static_cast<DataSet_Coords*>( *it );
    if (!ds->Meta().Aspect().empty()) {
      if (solventBoxName_ == ds->Meta().Aspect()) {
        solventUnit = ds;
        break;
      }
    }
  }
  // Name
  if (solventUnit == 0) {
    for (DataSetList::const_iterator it = sets.begin(); it != sets.end(); ++it)
    {
      DataSet_Coords* ds = static_cast<DataSet_Coords*>( *it );
      if (solventBoxName_ == ds->Meta().Name()) {
        solventUnit = ds;
        break;
      }
    }
  }
  if (solventUnit != 0)
    mprintf("\tSolvent unit: %s\n", solventUnit->legend());
  else
    mprinterr("Error: Could not get solvent unit named %s\n", solventBoxName_.c_str());

  return solventUnit;
}

/** Atom default radius in Angstroms from LEaP */
const double Solvate::ATOM_DEFAULT_RADIUS_ = 1.5;

/** Get radii for atoms in topology */
std::vector<double> Solvate::getAtomRadii(double& maxR, Topology const& topOut,
                                          Cpptraj::Parm::ParameterSet const& set0)
const
{
  using namespace Cpptraj::Parm;
  maxR = 0;
  std::vector<double> Radii;
  Radii.reserve( topOut.Natom() );
  for (int at = 0; at < topOut.Natom(); at++)
  {
    // Get radius
    double atom_radius = 0.0;
    //bool has_vdw = false;
    if (topOut[at].HasType()) {
      ParmHolder<AtomType>::const_iterator it = set0.AT().GetParam( TypeNameHolder(topOut[at].Type()) );
      if (it != set0.AT().end() && it->second.HasLJ()) {
        atom_radius = it->second.LJ().Radius();
        //has_vdw = true;
      }
    }
    if (atom_radius < 0.1) {
      if (topOut[at].Element() == Atom::HYDROGEN)
        atom_radius = 1.0;
      else
        atom_radius = ATOM_DEFAULT_RADIUS_;
    }
    Radii.push_back( atom_radius );
    maxR = std::max(maxR, atom_radius);
  }
  return Radii;
}

/** Set VDW bounding box. */
int Solvate::setVdwBoundingBox(double& boxX, double& boxY, double& boxZ,
                               std::vector<double> const& Radii,
                               Frame& frameOut, bool orient)
const
{
  if (orient) {
    // Center on origin and align principal axes.
    Vec3 ctr = frameOut.VGeometricCenter();
    ctr.Neg();
    frameOut.Translate(ctr);

    Matrix_3x3 Inertia;
    Vec3 Eval;
    frameOut.CalculateInertia( AtomMask(0, frameOut.Natom()), Inertia );
    Inertia.Diagonalize( Eval );
    //mprintf("Eigenvalues: %f %f %f\n", Eval[0], Eval[1], Eval[2]);

    // Check the handedness of the diagonalized matrix
    // If it is the wrong hand then it will transform the unit into a mirror image.
    Vec3 vPos = Inertia.Row1().Cross( Inertia.Row2() );
    //mprintf("VPOS= %f %f %f\n", vPos[0], vPos[1], vPos[2]);
    double dDot = vPos * Inertia.Row3();
    //mprintf( "The handedness of the transformation is (+1=Right): %f\n", dDot );

    // If the handedness of the matrix is wrong then change it
    if ( dDot < 0.0 ) {
      Inertia[6] *= -1.0;
      Inertia[7] *= -1.0;
      Inertia[8] *= -1.0;
    }

    //frameOut.Rotate( Inertia );
    // CPPTRAJ inverse rotate consistent with MatrixTimesVector from LEAP
    frameOut.InverseRotate( Inertia );
  }

  // Set vdw bounding box
  double Xmin = 0;
  double Ymin = 0;
  double Zmin = 0;
  double Xmax = 0;
  double Ymax = 0;
  double Zmax = 0;

  for (int at = 0; at < frameOut.Natom(); at++)
  {
    // Get radius
    double atom_radius = Radii[at];

    const double* XYZ = frameOut.XYZ(at);
    //mprintf("DBG: %12.4f %12.4f %12.4f %12.4f\n", XYZ[0], XYZ[1], XYZ[2], atom_radius);
    double dXp = XYZ[0] + atom_radius;
    double dYp = XYZ[1] + atom_radius;
    double dZp = XYZ[2] + atom_radius;
    double dXm = XYZ[0] - atom_radius;
    double dYm = XYZ[1] - atom_radius;
    double dZm = XYZ[2] - atom_radius;
    if (at == 0) {
      Xmin = dXm;
      Ymin = dYm;
      Zmin = dZm;
      Xmax = dXp;
      Ymax = dYp;
      Zmax = dZp;
    } else {
      if (dXm < Xmin) Xmin = dXm;
      if (dYm < Ymin) Ymin = dYm;
      if (dZm < Zmin) Zmin = dZm;
      if (dXp > Xmax) Xmax = dXp;
      if (dYp > Ymax) Ymax = dYp;
      if (dZp > Zmax) Zmax = dZp;
    }
  }
  //mprintf("Min= %f %f %f  Max= %f %f %f\n", Xmin, Ymin, Zmin, Xmax, Ymax, Zmax);
  // Define box
  boxX = Xmax - Xmin;
  boxY = Ymax - Ymin;
  boxZ = Zmax - Zmin;

  // Define center
  Vec3 toCenter( -(Xmin + 0.5 * boxX),
                 -(Ymin + 0.5 * boxY),
                 -(Zmin + 0.5 * boxZ) );
  //mprintf("ToolCenterUnitByRadii translate vector is %f %f %f\n", toCenter[0], toCenter[1], toCenter[2]);
  // Translate to origin
  frameOut.Translate(toCenter);

  return 0;
}

/** Set VDW bounding box */
int Solvate::SetVdwBoundingBox(Topology& topOut, Frame& frameOut, Cpptraj::Parm::ParameterSet const& set0)
{
  // Set vdw box
  double soluteMaxR;
  std::vector<double> soluteRadii = getAtomRadii(soluteMaxR, topOut, set0) ;
  double boxX, boxY, boxZ;
  if (setVdwBoundingBox(boxX, boxY, boxZ, soluteRadii, frameOut, isotropic_)) {
    mprinterr("Error: Setting vdw bounding box for %s failed.\n", topOut.c_str());
    return 1;
  }
  mprintf("  Solute vdw bounding box:              %-5.3f %-5.3f %-5.3f\n", boxX, boxY, boxZ);

  double dXWidth = boxX + bufferX_ * 2;
  double dYWidth = boxY + bufferY_ * 2;
  double dZWidth = boxZ + bufferZ_ * 2;

    mprintf("  Total bounding box for atom centers:  %5.3f %5.3f %5.3f\n", 
            dXWidth, dYWidth, dZWidth );

  // Setup box
  frameOut.ModifyBox().SetupFromXyzAbg(boxX, boxY, boxZ, 90.0, 90.0, 90.0);
  frameOut.BoxCrd().PrintInfo();
  topOut.SetParmBox( frameOut.BoxCrd() );
  mprintf("  Total vdw box size:%s%5.3f %5.3f %5.3f angstroms.\n", "                   ",
          frameOut.BoxCrd().Param(Box::X),
          frameOut.BoxCrd().Param(Box::Y),
          frameOut.BoxCrd().Param(Box::Z));
  mprintf("  Volume: %5.3lf A^3\n", frameOut.BoxCrd().CellVolume());
  return 0;
}

/** Scale buffer if needed to meet diagonal clearance. */
void Solvate::octBoxCheck(Frame const& frameOut) {
  if (frameOut.Natom() < 1) return;
  double dPBuf[4];

  dPBuf[0] = bufferX_;
  dPBuf[1] = bufferY_;
  dPBuf[2] = bufferZ_;
  dPBuf[3] = bufferD_;
  //mprintf("dPBuf %f %f %f %f\n", dPBuf[0], dPBuf[1], dPBuf[2], dPBuf[3]);

  const double* XYZ = frameOut.XYZ(0);
  double dXmax = XYZ[0];
  double dYmax = XYZ[1];
  double dZmax = XYZ[2];
  double dXmin = XYZ[0];
  double dYmin = XYZ[1];
  double dZmin = XYZ[2];

  for (int at = 1; at < frameOut.Natom(); at++)
  {
    XYZ = frameOut.XYZ(at);
    if      ( XYZ[0] > dXmax )
            dXmax = XYZ[0];
    else if ( XYZ[0] < dXmin )
            dXmin = XYZ[0];
    if      ( XYZ[1] > dYmax )
            dYmax = XYZ[1];
    else if ( XYZ[1] < dYmin )
            dYmin = XYZ[1];
    if      ( XYZ[2] > dZmax )
            dZmax = XYZ[2];
    else if ( XYZ[2] < dZmin )
            dZmin = XYZ[2];
  }

  // calc halfbox on centers
  double dXhalf = 0.5 * (dXmax - dXmin) + dPBuf[0];
  double dYhalf = 0.5 * (dYmax - dYmin) + dPBuf[1];
  double dZhalf = 0.5 * (dZmax - dZmin) + dPBuf[2];

  // find unit vector of diagonal
  double dX = dYhalf * dZhalf;
  double dY = dXhalf * dZhalf;
  double dZ = dXhalf * dYhalf;

  double dTmp = 1.0 / sqrt( dX*dX + dY*dY + dZ*dZ );

  double dXunit = dX * dTmp;
  double dYunit = dY * dTmp;
  double dZunit = dZ * dTmp;

  // find max atom distance from origin along the diagonal
  double dMax = 0.0;
  for (int at = 0; at < frameOut.Natom(); at++)
  {
    XYZ = frameOut.XYZ(at);
    dTmp = fabs(XYZ[0]) * dXunit +
           fabs(XYZ[1]) * dYunit +
           fabs(XYZ[2]) * dZunit;
    if ( dTmp > dMax )
      dMax = dTmp;
  }

  // calc distance of diagonal face from origin
  double dBmax = 0.5 * sqrt( dXhalf*dXhalf + dYhalf*dYhalf + dZhalf*dZhalf );

  // see if diagonal clearance is satisfied
  dTmp = dMax + dPBuf[3];
  if ( dTmp <= dBmax ) {
    if ( dPBuf[3] == 0.0 )
      mprintf("  (Diagonal clearance is %f)\n", dMax );
    return;
  }

  //  not satisfied: scale up box
  dTmp /= dBmax;
  mprintf("  Scaling up box by a factor of %f to meet diagonal cut criterion\n", dTmp );

  bufferX_ *= dTmp;
  bufferY_ *= dTmp;
  bufferZ_ *= dTmp;
}

/** Rotate 45 deg. around z axis, (90-tetra/2) around y axis, 90 around x axis.
  *    (1  0  0)    (cos2  0 -sin2)    (cos1 -sin1  0)
  *    (0  0 -1)    (   0  1     0)    (sin1  cos1  0)
  *    (0  1  0)    (sin2  0  cos2)    (   0     0  1)
  *
  *    cntr-clk       clock              clock
  *    Looking down + axis of rotation toward origin
  *
  * For isotropic truncated octahedral box only.
  */
void Solvate::ewald_rotate(Frame& frameOut, double& dPAngle)
{
  double tetra_angl = 2 * acos( 1. / sqrt(3.) );
  double pi = 3.1415927; // FIXME use Constants
  double phi = pi / 4.;
  double cos1 = cos(phi);
  double sin1 = sin(phi);
         phi = pi/2. - tetra_angl/2.;
  double cos2 = sqrt(2.)/sqrt(3.);
  double sin2=1./sqrt(3.);

  double t11= cos2*cos1;
  double t12=-cos2*sin1;
  double t13=-sin2;
  double t21=-sin2*cos1;
  double t22= sin2*sin1;
  double t23=-cos2;
  double t31= sin1;
  double t32= cos1;
  double t33=0;

  //lAtoms = lLoop( (OBJEKT)uUnit, ATOMS );
  //while ( (aAtom = (ATOM)oNext(&lAtoms)) != NULL ) {
  //      double  dX, dY, dZ;
  double* xyz = frameOut.xAddress();
  for (int at = 0; at != frameOut.Natom(); at++, xyz += 3)
  {
    double dX = xyz[0];
    double dY = xyz[1];
    double dZ = xyz[2];

    xyz[0] = t11*dX + t12*dY + t13*dZ;
    xyz[1] = t21*dX + t22*dY + t23*dZ;
    xyz[2] = t31*dX + t32*dY + t33*dZ;
  }
  dPAngle = tetra_angl*180./pi;
}

/** Create box, fill with solvent */
int Solvate::SolvateBox(Topology& topOut, Frame& frameOut, Cpptraj::Parm::ParameterSet const& set0,
                        DataSet_Coords& SOLVENTBOX)
{
  // Sanity check
  if (topOut.Natom() != frameOut.Natom()) {
    mprinterr("Internal Error: Solvate::SolvateBox(): Topology %s #atoms %i != frame #atoms %i\n",
              topOut.c_str(), topOut.Natom(), frameOut.Natom());
    return 1;
  }

  // TODO Remove any existing box info?
  //if (frameOut.BoxCrd.HasBox())
  // TODO principal align

  // Set vdw box
  double soluteMaxR;
  std::vector<double> soluteRadii = getAtomRadii(soluteMaxR, topOut, set0) ;
  double boxX, boxY, boxZ;
  if (setVdwBoundingBox(boxX, boxY, boxZ, soluteRadii, frameOut, isotropic_)) {
    mprinterr("Error: Setting vdw bounding box for %s failed.\n", topOut.c_str());
    return 1;
  }
  // Check if buffers need to be increased for trunc oct.
  if (doTruncatedOct_)
    octBoxCheck( frameOut );
  mprintf("  Solute vdw bounding box:              %-5.3f %-5.3f %-5.3f\n", boxX, boxY, boxZ);

  double dXWidth = boxX + bufferX_ * 2;
  double dYWidth = boxY + bufferY_ * 2;
  double dZWidth = boxZ + bufferZ_ * 2;

  if (isotropic_) {
    double dTemp = dXWidth * dYWidth * dZWidth;

    double dMax = std::max(dXWidth, dYWidth);
    dMax = std::max(dMax, dZWidth);
    dXWidth = dYWidth = dZWidth = dMax;

    dTemp = (dMax * dMax * dMax - dTemp ) / dTemp;

    mprintf("  Total bounding box for atom centers:  %5.3f %5.3f %5.3f\n", 
            dXWidth, dYWidth, dZWidth );
    mprintf("      (box expansion for 'iso' is %5.1lf%%)\n", dTemp * 100.0 );

     // To make the actual clip right, 'iso' the solute box
    dTemp = std::max(boxX, boxY);
    dTemp = std::max(dTemp, boxZ);
    //dXBox = dYBox = dZBox = dTemp;
    boxX = boxY = boxZ = dTemp;
  } else
    mprintf("  Total bounding box for atom centers:  %5.3f %5.3f %5.3f\n", 
            dXWidth, dYWidth, dZWidth );

  if (clip_) {
    //  If the solvated system should be clipped to the exact
    //      size the user specified then note the criterion
    //      & dimensions (for 0,0,0-centered system)
    //iCriteria |= TOOLOUTSIDEOFBOX;
    clipX_ = 0.5 * boxX + bufferX_;
    clipY_ = 0.5 * boxY + bufferY_;
    clipZ_ = 0.5 * boxZ + bufferZ_;
    //mprintf("cCriteria: %f %f %f\n", clipX_, clipY_, clipZ_);
  }

  // TODO check COORDS size
  Frame solventFrame = SOLVENTBOX.AllocateFrame();
  SOLVENTBOX.GetFrame(0, solventFrame);

  // Set vdw box for solvent
  double solventMaxR;
  std::vector<double> solventRadii = getAtomRadii(solventMaxR, SOLVENTBOX.Top(), set0);
  double solventX, solventY, solventZ;
  if (solventFrame.BoxCrd().HasBox()) {
    // Use input box lengths
    // TODO check ortho?
    solventX = solventFrame.BoxCrd().Param(Box::X);
    solventY = solventFrame.BoxCrd().Param(Box::Y);
    solventZ = solventFrame.BoxCrd().Param(Box::Z);
  } else {
    if (setVdwBoundingBox(solventX, solventY, solventZ, solventRadii, solventFrame, false)) {
      mprinterr("Error: Setting vdw bounding box for %s failed.\n", topOut.c_str());
      return 1;
    }
  }
  mprintf("  Solvent unit box:                     %5.3f %5.3f %5.3f\n", solventX, solventY, solventZ);

  // See how many solvent boxes required in each dimension

  int iX = (int)( dXWidth / solventX ) + 1;
  int iY = (int)( dYWidth / solventY ) + 1;
  int iZ = (int)( dZWidth / solventZ ) + 1;

  //  Calculate the center of the first solvent box 
  //  (the one that goes in the max XYZ corner), given
  //  that the solute is centered at 0,0,0
  double dXStart = 0.5 * solventX * (double) (iX-1);
  double dYStart = 0.5 * solventY * (double) (iY-1);
  double dZStart = 0.5 * solventZ * (double) (iZ-1);

             /* If the caller wants a solvent shell then */
             /* make sure that the box used to find interesting solute */
             /* spheres takes into account the dFarness parameter */
             /* so that there are at least some solute spheres in */
             /* the interesting list to check against solvent */

 //if ( bShell ) 
 //    dBuffer = dFarness;
 //else 
 //    dBuffer = 0.0;

  addSolventUnits(iX, iY, iZ, soluteMaxR, dXStart, dYStart, dZStart, solventX, solventY, solventZ,
                  solventFrame, SOLVENTBOX.Top(), frameOut, topOut,
                  soluteRadii, solventRadii);

  // Define the size of the new solvent/solute system
  soluteRadii = getAtomRadii(soluteMaxR, topOut, set0) ;
  if (setVdwBoundingBox(boxX, boxY, boxZ, soluteRadii, frameOut, false)) {
    mprinterr("Error: Setting vdw bounding box for solute/solvent system failed.\n", topOut.c_str());
    return 1;
  }
  if (doTruncatedOct_) {
    double dAngle = 0;
    ewald_rotate(frameOut, dAngle);
    //mprintf("EwaldRotate: %f\n", dAngle);
    // Add an angstrom to the desired box size rather than using the bounding box size
    dAngle = clipX_ + .5;
    boxX = boxY = boxZ = dAngle * sqrt(3.0) * 0.5;
  }

  // Setup box
  double boxBeta = 90.0;
  if (doTruncatedOct_) {
    boxBeta = Box::TruncatedOctAngle();
    boxX *= 2.0;
    boxY *= 2.0;
    boxZ *= 2.0;
  }
  //mprintf("Max: %f %f %f\n", boxX, boxY, boxZ);
  frameOut.ModifyBox().SetupFromXyzAbg(boxX, boxY, boxZ, boxBeta, boxBeta, boxBeta);
  frameOut.BoxCrd().PrintInfo();
  topOut.SetParmBox( frameOut.BoxCrd() );
  mprintf("  Total vdw box size:%s%5.3f %5.3f %5.3f angstroms.\n", "                   ",
          frameOut.BoxCrd().Param(Box::X),
          frameOut.BoxCrd().Param(Box::Y),
          frameOut.BoxCrd().Param(Box::Z));
  mprintf("  Volume: %5.3lf A^3\n", frameOut.BoxCrd().CellVolume());
  // Sum mass
  double sumMass = 0.0;
  for (int at = 0; at < topOut.Natom(); at++) {
    if (topOut[at].HasType()) {
      Cpptraj::Parm::ParmHolder<AtomType>::const_iterator it = set0.AT().GetParam( TypeNameHolder(topOut[at].Type()) );
      if (it != set0.AT().end()) {
        sumMass += it->second.Mass();
      }
    }
  }
  if (sumMass > 0.0) {
    mprintf("  Total mass %5.3f amu,  Density %5.3lf g/cc\n", sumMass, sumMass / (frameOut.BoxCrd().CellVolume() * 0.602204));
  } else {
    mprintf("Warning: Mass could not be determined, so density unknown (i.e. type of all atoms could not be found)\n");
  }
  // Center if needed
  if (center_) {
    double dX2 = boxX * 0.5;
    double dY2 = boxY * 0.5;
    double dZ2 = boxZ * 0.5;
    double* xptr = frameOut.xAddress();
    for (int at = 0; at < frameOut.Natom(); at++, xptr += 3)
    {
      xptr[0] += dX2;
      xptr[1] += dY2;
      xptr[2] += dZ2;
    }
  }

  return 0;
}

/** \return a list of solute atoms that might clash with solvent. */
int Solvate::findCloseSoluteAtoms(std::vector<int>& closeSoluteAtoms, double soluteMaxR,
                                  int firstSolventAtom, Frame const& frameOut, Vec3 const& vCenter,
                                  double dXWidth, double dYWidth, double dZWidth)
const
{
  closeSoluteAtoms.clear();
//  mprintf( "Searching for close solute atoms, buffer zone %f solute max %f closeness %f\n", 0.0, soluteMaxR, closeness_);
  // Determine clearance from the box for testing whether
  // a solute atom might contact an atom in the box. Assumes
  // solvent box includes vdw.
  double dTemp = soluteMaxR * 2.5 * closeness_;

  // FIXME do buffer zone?
  // when building a shell, solute atoms within the
  // shell distance are also of interest
  //  dTemp += dBufferZone;

  double dXmin = vCenter[0] - dXWidth/2.0 - dTemp;
  double dXmax = vCenter[0] + dXWidth/2.0 + dTemp;
  double dYmin = vCenter[1] - dYWidth/2.0 - dTemp;
  double dYmax = vCenter[1] + dYWidth/2.0 + dTemp;
  double dZmin = vCenter[2] - dZWidth/2.0 - dTemp;
  double dZmax = vCenter[2] + dZWidth/2.0 + dTemp;
//  mprintf("Search clearances %f Min= %f %f %f Max= %f %f %f\n", dTemp, dXmin, dYmin, dZmin, dXmax, dYmax, dZmax);

  // Loop over solute atoms
  for (int at = 0; at < firstSolventAtom; at++)
  {
    const double* XYZ = frameOut.XYZ(at);
    if ( XYZ[0] < dXmin ) continue;
    if ( XYZ[0] > dXmax ) continue;
    if ( XYZ[1] < dYmin ) continue;
    if ( XYZ[1] > dYmax ) continue;
    if ( XYZ[2] < dZmin ) continue;
    if ( XYZ[2] > dZmax ) continue;

    // all atom.coords inside solvent limit, so add to list
//    mprintf("Found an interesting sphere %f %f %f\n", XYZ[0], XYZ[1], XYZ[2] );
    closeSoluteAtoms.push_back( at );
  }
  return 0;
}

/** LEaP's closeness modifier. Was just set to 1. Included here just in case. */
const double Solvate::CLOSENESSMODIFIER_ = 1.0;

/** Determine which solvent residues will not clash with existing solute atoms. */
int Solvate::determineValidSolventResidues(std::vector<int>& validSolventResidues,
                                           std::vector<int> const& closeSoluteAtoms,
                                           Frame const& solventFrame, Topology const& solventTop,
                                           Frame const& frameOut,
                                           std::vector<double> const& soluteRadii,
                                           std::vector<double> const& solventRadii)
const
{
  validSolventResidues.clear();
  for (int vres = 0; vres < solventTop.Nres(); vres++)
  {
    Residue const& solventRes = solventTop.Res(vres);
    bool collision = false;
    // Check each atom in the solvent residue against list of close solute atoms
    for (int vat = solventRes.FirstAtom(); vat != solventRes.LastAtom(); vat++)
    {
      const double* VXYZ = solventFrame.XYZ(vat);
      // First check for clipping
      if (clip_) {
        double dXabs = fabs(VXYZ[0]);
        double dYabs = fabs(VXYZ[1]);
        double dZabs = fabs(VXYZ[2]);
        //if ( dXabs >= clipX_ ) { collision = true; mprintf("CLIP %12.4f %12.4f %12.4f\n",VXYZ[0],VXYZ[1],VXYZ[2]); break; }
        //if ( dYabs >= clipY_ ) { collision = true; mprintf("CLIP %12.4f %12.4f %12.4f\n",VXYZ[0],VXYZ[1],VXYZ[2]); break; }
        //if ( dZabs >= clipZ_ ) { collision = true; mprintf("CLIP %12.4f %12.4f %12.4f\n",VXYZ[0],VXYZ[1],VXYZ[2]); break; }
        if ( dXabs >= clipX_ ) { collision = true; break; }
        if ( dYabs >= clipY_ ) { collision = true; break; }
        if ( dZabs >= clipZ_ ) { collision = true; break; }
        // Check if atom falls outside of oct/diagonal clip
        if (doTruncatedOct_) {
          double dX = 0.5 * (clipX_ - dXabs) / clipX_;
                 dX = fabs( dX - 0.5 );
          double dY = 0.5 * (clipY_ - dYabs) / clipY_;
                 dY = fabs( dY - 0.5 );
          double dZ = 0.5 * (clipZ_ - dZabs) / clipZ_;
                 dZ = fabs( dZ - 0.5 );
          if ( ( dX + dY + dZ )  >  0.75 ) { collision = true; break; }
        }
      }
      double dR = solventRadii[vat] * closeness_ * CLOSENESSMODIFIER_;
      // Loop over close solute atoms, check for clash
      for (std::vector<int>::const_iterator uat = closeSoluteAtoms.begin(); uat != closeSoluteAtoms.end(); ++uat)
      {
        const double* UXYZ = frameOut.XYZ(*uat);
        double dX = VXYZ[0] - UXYZ[0];
        double dY = VXYZ[1] - UXYZ[1];
        double dZ = VXYZ[2] - UXYZ[2];

        double dist2 = dX*dX + dY*dY + dZ*dZ;
//        mprintf("RADIUS %f\n", soluteRadii[*uat]);
        double dRadii = dR + soluteRadii[*uat];
        dRadii *= dRadii;

        if (dist2 < dRadii) {
          collision = true;
          //mprintf("OVERLAP %12.4f %12.4f %12.4f %12.4f %12.4f\n", VXYZ[0],VXYZ[1],VXYZ[2], dist2, dRadii);
          break;
        }
      } // END loop over close solute atoms
      if (collision) break;
    } // END loop over solvent residue atoms
    if (!collision)
      validSolventResidues.push_back( vres );
  } // END loop over solvent residues
  return 0;
}

/** Add as many solvent units as needed to complete solvation. */
int Solvate::addSolventUnits(int numX, int numY, int numZ, double soluteMaxR,
                             double dXStart, double dYStart, double dZStart,
                             double dXSolvent, double dYSolvent, double dZSolvent,
                             Frame& solventFrame, Topology const& solventTop,
                             Frame& frameOut, Topology& topOut,
                             std::vector<double> const& soluteRadii,
                             std::vector<double> const& solventRadii)
const
{
  if (debug_ > 0) {
    mprintf( "Max R = %f\n", soluteMaxR);
    mprintf( "The number of boxes:  x=%2d  y=%2d  z=%2d\n", numX, numY, numZ );
  }
  int NboxesToAdd = numX * numY * numZ;
  int NatomsToAdd = NboxesToAdd * solventTop.Natom();
  mprintf("  Will add %i solvent boxes, max %i atoms.\n", NboxesToAdd, NatomsToAdd);
  frameOut.IncreaseX( NatomsToAdd );

  // Current solvent unit center
  Vec3 currentSolventCenter = solventFrame.VGeometricCenter();

  int firstSolventAtom = topOut.Natom();
  if (debug_ > 0)
    mprintf("DEBUG: First solvent atom is %i\n", firstSolventAtom+1);

  std::vector<int> bondedAtoms;
  bondedAtoms.reserve(12); // Reserve for 6 bonds

  std::vector<int> closeSoluteAtoms;
  std::vector<int> validSolventResidues;

  double dX = dXStart;
  for ( int ix=0; ix < numX; ix++, dX -= dXSolvent ) {
    double dY = dYStart;
    for ( int iy=0; iy < numY; iy++, dY -= dYSolvent ) {
      double dZ = dZStart;
      for ( int iz=0; iz < numZ; iz++, dZ -= dZSolvent ) {
        //mprintf( "Adding box at: x=%d  y=%d  z=%d\n", ix, iy, iz);
        Vec3 vPos(dX, dY, dZ);

        findCloseSoluteAtoms(closeSoluteAtoms, soluteMaxR, firstSolventAtom, frameOut, vPos,
                             dXSolvent, dYSolvent, dZSolvent);

        //mprintf( "Center of solvent box is: %lf, %lf, %lf\n", dX, dY, dZ );
        Vec3 trans( dX - currentSolventCenter[0],
                    dY - currentSolventCenter[1],
                    dZ - currentSolventCenter[2] );
        solventFrame.Translate(trans);
        //Vec3 debugVec = solventFrame.VGeometricCenter();
        //debugVec.Print("DEBUG: check solvent center");

        determineValidSolventResidues(validSolventResidues, closeSoluteAtoms, solventFrame, solventTop,
                                      frameOut, soluteRadii, solventRadii); 

        // Update the current solvent center
        currentSolventCenter[0] = dX;
        currentSolventCenter[1] = dY;
        currentSolventCenter[2] = dZ;

        // Add valid residues from solvent unit to output topology for this cube
        int currentResNum = topOut.Nres();
        //for (int ires = 0; ires != solventTop.Nres(); ires++)
        for (std::vector<int>::const_iterator ires = validSolventResidues.begin();
                                              ires != validSolventResidues.end(); ++ires)
        {
          int atomOffset = topOut.Natom();
          Residue solventRes = solventTop.Res(*ires);
          solventRes.SetOriginalNum( currentResNum++ );
          bondedAtoms.clear();
          for (int iat = solventRes.FirstAtom(); iat != solventRes.LastAtom(); iat++)
          {
            Atom solventAtom = solventTop[iat];
            // Save solvent bonds
            for (Atom::bond_iterator bat = solventAtom.bondbegin(); bat != solventAtom.bondend(); ++bat) {
              if (*bat > iat) {
                bondedAtoms.push_back(  iat + atomOffset - solventRes.FirstAtom() );
                bondedAtoms.push_back( *bat + atomOffset - solventRes.FirstAtom() );
              }
            }
            solventAtom.ClearBonds(); // FIXME AddTopAtom should clear
            topOut.AddTopAtom( solventAtom, solventRes );
            // Add PDB info if the topology already has it.
            if (!topOut.Bfactor().empty()) topOut.AddBfactor( 0 );
            if (!topOut.Occupancy().empty()) topOut.AddOccupancy( 1 );
            if (!topOut.PdbSerialNum().empty()) topOut.AddPdbSerialNum( topOut.Natom() );
            const double* VXYZ = solventFrame.XYZ(iat);
            frameOut.AddXYZ( VXYZ );
          } // END loop over solvent atoms
          // Add bonds
          for (std::vector<int>::const_iterator it = bondedAtoms.begin(); it != bondedAtoms.end(); ++it) {
            int at0 = *it;
            ++it;
            topOut.AddBond( at0, *it );
          }
        } // END loop over solvent unit residues
        // Append solvent frame
        //frameOut.AppendFrame( solventFrame );
      } // END loop over Z
    } // END loop over Y
  } // END loop over X

  return 0;
}
