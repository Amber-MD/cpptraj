#include "Solvate.h"
#include "../ArgList.h"
#include "../CpptrajStdio.h"
#include "../DataSet_Coords.h"
#include "../DataSetList.h"
#include "../Frame.h"
#include "../Topology.h"
#include "../Parm/ParameterSet.h"
#include <algorithm> //std::max

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
Solvate::Solvate() :
  debug_(0),
  bufferX_(0),
  bufferY_(0),
  bufferZ_(0),
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
  } else {
    bufferX_ = argIn.getKeyDouble("bufx", defaultBuffer);
    bufferY_ = argIn.getKeyDouble("bufy", defaultBuffer);
    bufferZ_ = argIn.getKeyDouble("bufz", defaultBuffer);
  }
  if (bufferX_ < 0 || bufferY_ < 0 || bufferZ_ < 0) {
    mprinterr("Error: Either 'buffer' or 'bufx/bufy/bufx' must be specified and >= 0\n");
    return 1;
  }
  return 0;
}

/** Initialize arguments. */
int Solvate::InitSolvate(ArgList& argIn, int debugIn) {
  debug_ = debugIn;

  if (getBufferArg(argIn, -1.0)) return 1;

  isotropic_ = argIn.hasKey("iso");

  solventBoxName_ = argIn.GetStringKey("solventbox");
  if (solventBoxName_.empty()) {
    mprinterr("Error: Specify solvent box unit name with 'solventbox'\n");
    return 1;
  }

  closeness_ = argIn.getKeyDouble("closeness", 1.0);
  center_ = !argIn.hasKey("nocenter");

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
  mprintf("\tSolvent isotropic=%i  clip=%i  center_=%i\n", (int)isotropic_, (int)clip_, (int)center_);
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
                               Frame& frameOut)
const
{
  if (isotropic_) {
    // Center on origin and align principal axes.
    Vec3 ctr = frameOut.VGeometricCenter();
    ctr.Neg();
    frameOut.Translate(ctr);

    Matrix_3x3 Inertia;
    Vec3 Eval;
    frameOut.CalculateInertia( AtomMask(0, frameOut.Natom()), Inertia );
    Inertia.Diagonalize( Eval );
    printf("Eigenvalues: %f %f %f\n", Eval[0], Eval[1], Eval[2]);

    frameOut.Rotate( Inertia );
    //frameOut.ModifyBox().RotateUcell( Inertia );
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
  if (setVdwBoundingBox(boxX, boxY, boxZ, soluteRadii, frameOut)) {
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
  if (setVdwBoundingBox(boxX, boxY, boxZ, soluteRadii, frameOut)) {
    mprinterr("Error: Setting vdw bounding box for %s failed.\n", topOut.c_str());
    return 1;
  }
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
/*
    if ( bOct ) {
        // maybe allow oct clip on integer boxes someday.. but for now:
        if ( !bClip )
                DFATAL(( "oct but no clip\n" ));
        iCriteria |= TOOLOUTSIDEOFOCTBOX;
    }
*/
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
    if (setVdwBoundingBox(solventX, solventY, solventZ, solventRadii, solventFrame)) {
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
  if (setVdwBoundingBox(boxX, boxY, boxZ, soluteRadii, frameOut)) {
    mprinterr("Error: Setting vdw bounding box for solute/solvent system failed.\n", topOut.c_str());
    return 1;
  }
  // Setup box
  frameOut.ModifyBox().SetupFromXyzAbg(boxX, boxY, boxZ, 90.0, 90.0, 90.0);
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
        //if ( fabs(VXYZ[0]) >= clipX_ ) { collision = true; mprintf("CLIP %12.4f %12.4f %12.4f\n",VXYZ[0],VXYZ[1],VXYZ[2]); break; }
        //if ( fabs(VXYZ[1]) >= clipY_ ) { collision = true; mprintf("CLIP %12.4f %12.4f %12.4f\n",VXYZ[0],VXYZ[1],VXYZ[2]); break; }
        //if ( fabs(VXYZ[2]) >= clipZ_ ) { collision = true; mprintf("CLIP %12.4f %12.4f %12.4f\n",VXYZ[0],VXYZ[1],VXYZ[2]); break; }
        if ( fabs(VXYZ[0]) >= clipX_ ) { collision = true; break; }
        if ( fabs(VXYZ[1]) >= clipY_ ) { collision = true; break; }
        if ( fabs(VXYZ[2]) >= clipZ_ ) { collision = true; break; }
      }
      double dR = solventRadii[vat] * closeness_ * CLOSENESSMODIFIER_;
      // Loop over close solute atoms, check fir ckasg
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
  mprintf( "Max R = %f\n", soluteMaxR);
  mprintf( "The number of boxes:  x=%2d  y=%2d  z=%2d\n", numX, numY, numZ );
  int NboxesToAdd = numX * numY * numZ;
  int NatomsToAdd = NboxesToAdd * solventTop.Natom();
  mprintf("Will add %i boxes, %i atoms.\n", NboxesToAdd, NatomsToAdd);
  frameOut.IncreaseX( NatomsToAdd );

  // Current solvent unit center
  Vec3 currentSolventCenter = solventFrame.VGeometricCenter();

  int firstSolventAtom = topOut.Natom();
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
