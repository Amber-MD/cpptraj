#include <cmath> // floor
#include "ImageRoutines.h"
#include "DistRoutines.h"
#include "CpptrajStdio.h"
#include "Image_List.h"
#include "Image_List_Pair_CoM.h"
#include "Image_List_Pair_Geom.h"
#include "Image_List_Pair_First.h"
#include "Image_List_Unit_CoM.h"
#include "Image_List_Unit_Geom.h"
#include "Image_List_Unit_First.h"
#include "Image_List_Mask.h"

/** \return Empty list for imaging. */
Image::List* Image::CreateImageList(Mode modeIn, bool useMass, bool center) {
  Image::List* listOut = 0;
  switch (modeIn) {
    case BYMOL :
      if (center) {
        if (useMass)
          listOut = new Image::List_Unit_CoM();
        else
          listOut = new Image::List_Unit_Geom();
      } else
        listOut = new Image::List_Unit_First();
      break;
    case BYRES :
      if (center) {
        if (useMass)
          listOut = new Image::List_Pair_CoM();
        else
          listOut = new Image::List_Pair_Geom();
      } else
        listOut = new Image::List_Pair_First();
      break;
    case BYATOM :
      listOut = new Image::List_Mask();
      break;
  }
  if (listOut == 0) {
    mprinterr("Internal Error: Could not allocate image list over %ss\n", ModeString(modeIn));
  }
  return listOut;
}

// Image::CreateImageList() 
/** \return list of entities to be imaged based on given mode.
  */
Image::List* Image::CreateImageList(Topology const& Parm, Mode modeIn,
                                    std::string const& maskExpression,
                                    bool useMass, bool center)
{
  Image::List* listOut = CreateImageList(modeIn, useMass, center);
  if (listOut != 0) {
    if (listOut->SetupList(Parm, maskExpression)) {
      mprinterr("Error: Could not set up image list for '%s'\n", maskExpression.c_str());
      delete listOut;
      return 0;
    }
    //mprintf("DEBUG: Image list for '%s' over %u %ss.\n",
    //        maskExpression.c_str(), listOut->nEntities(), ModeString(modeIn));
  }
  return listOut;
}

// -----------------------------------------------------------------------------
// Image::SetupTruncoct()
/** Set up centering if putting nonortho cell into familiar trunc. oct. shape.
  * \param frameIn Frame to set up for.
  * \param ComMask If not null center is calcd w.r.t. center of atoms in mask.
  * \param useMass If true calculate COM, otherwise calc geometric center.
  * \param origin If true and ComMask is null use origin, otherwise use box center.
  * \return Coordinates of center.
  */
Vec3 Image::SetupTruncoct( Frame const& frameIn, AtomMask* ComMask, bool useMass, bool origin)
{
  if (ComMask!=0) {
    // Use center of atoms in mask
    if (useMass)
      return frameIn.VCenterOfMass( *ComMask );
    else
      return frameIn.VGeometricCenter( *ComMask );
  } else if (!origin) {
    // Use box center
    return frameIn.BoxCrd().Center(); 
  }
  //fprintf(stdout,"DEBUG: fcom = %lf %lf %lf\n",fcom[0],fcom[1],fcom[2]);
  return Vec3(0.0, 0.0, 0.0); // Default is origin {0,0,0}
}

// Image::Nonortho()
/** \param frameIn Frame to image.
  * \param origin If true image w.r.t. coordinate origin.
  * \param fcom If truncoct is true, calc distance w.r.t. this coordinate.
  * \param ucell Unit cell matrix.
  * \param recip Reciprocal coordinates matrix.
  * \param truncoct If true imaging will occur using truncated octahedron shape.
  * \param center If true image w.r.t. center coords, otherwise use first atom coords.
  * \param useMass If true use COM, otherwise geometric center.
  * \param AtomPairs Atom pairs to image.
  */
void Image::Nonortho(Frame& frameIn, bool origin, Vec3 const& fcom, Vec3 const& offIn, 
                     Matrix_3x3 const& ucell, Matrix_3x3 const& recip,
                     bool truncoct, List const& AtomPairs)
{
  Vec3 Coord;
  Vec3 offset = ucell.TransposeMult( offIn );
  double min = -1.0;

  if (truncoct)
    min = 100.0 * (frameIn.BoxCrd().Param(Box::X)*frameIn.BoxCrd().Param(Box::X)+
                   frameIn.BoxCrd().Param(Box::Y)*frameIn.BoxCrd().Param(Box::Y)+
                   frameIn.BoxCrd().Param(Box::Z)*frameIn.BoxCrd().Param(Box::Z));

  // Loop over atom pairs
  for (unsigned int idx = 0; idx != AtomPairs.nEntities(); idx++)
  {
    Coord = AtomPairs.GetCoord( idx, frameIn );

    // boxTrans will hold calculated translation needed to move atoms back into box
    Vec3 boxTrans = Nonortho(Coord, truncoct, origin, ucell, recip, fcom, min) + offset;

    AtomPairs.DoTranslation( frameIn, idx, boxTrans );
  } // END loop over atom pairs
}

// Image::Nonortho()
/** \param Coord Coordinate to image.
  * \param truncoct If true, image in truncated octahedral shape.
  * \param origin If true, image w.r.t. coordinate origin.
  * \param ucell Unit cell matrix.
  * \param recip Reciprocal coordinates matrix.
  * \param fcom If truncoct, image translated coordinate w.r.t. this coord.
  * \return Vector containing image translation.
  */
Vec3 Image::Nonortho(Vec3 const& Coord, bool truncoct, 
                     bool origin, Matrix_3x3 const& ucell, Matrix_3x3 const& recip, 
                     Vec3 const& fcom, double min)
{
  int ixyz[3];

  Vec3 fc = recip * Coord;

  if ( origin )
    fc += 0.5; 

  Vec3 boxTransOut = ucell.TransposeMult( Vec3(floor(fc[0]), floor(fc[1]), floor(fc[2])) );
  boxTransOut.Neg();

  // Put into familiar trunc. oct. shape
  if (truncoct) {
    Vec3 TransCoord = recip * (Coord + boxTransOut);
    Vec3 f2 = recip * fcom;

    if (origin) {
      TransCoord += 0.5;
      f2 += 0.5;
    }

    DIST2_ImageNonOrthoRecip(TransCoord, f2, min, ixyz, ucell);
    if (ixyz[0] != 0 || ixyz[1] != 0 || ixyz[2] != 0) {
      boxTransOut += ucell.TransposeMult( ixyz );
      //if (debug > 2)
      //  mprintf( "  IMAGING, FAMILIAR OFFSETS ARE %i %i %i\n", 
      //          ixyz[0], ixyz[1], ixyz[2]);
    }
  }
  return boxTransOut;
}

// -----------------------------------------------------------------------------
// Image::SetupOrtho()
/** \param boxIn Box coordinates of Frame to image.
  * \param bp Output: Box + boundary.
  * \param bm Output: Box - boundary.
  * \param origin If true, image w.r.t. coordinate origin, otherwise box center.
  * \return 1 if box lengths are zero, 0 if setup completed successfully.
  */
int Image::SetupOrtho(Box const& boxIn, Vec3& bp, Vec3& bm, bool origin) {
  // Set up boundary information for orthorhombic cell
  if (origin) {
    bp = boxIn.Center();
    bm.SetVec( -bp[0], -bp[1], -bp[2] );
  } else {
    bp.SetVec( boxIn.Param(Box::X), boxIn.Param(Box::Y), boxIn.Param(Box::Z)  );
    bm.Zero();
  }
  if (bp.IsZero()) return 1;
  return 0;
}

// Image::Ortho()
/** \param frameIn Frame to image.
  * \param bp Box + boundary.
  * \param bm Box - boundary.
  * \param center If true image w.r.t. center of atoms, otherwise first atom.
  * \param useMass If true calc center of mass, otherwise geometric center.
  */
void Image::Ortho(Frame& frameIn, Vec3 const& bp, Vec3 const& bm, Vec3 const& offIn,
                  List const& AtomPairs)
{
  Vec3 Coord;
  Vec3 offset(offIn[0] * frameIn.BoxCrd().Param(Box::X),
              offIn[1] * frameIn.BoxCrd().Param(Box::Y),
              offIn[2] * frameIn.BoxCrd().Param(Box::Z));
  // Loop over atom pairs
  for (unsigned int idx = 0; idx != AtomPairs.nEntities(); idx++)
  {
    Coord = AtomPairs.GetCoord(idx, frameIn);

    // boxTrans will hold calculated translation needed to move atoms back into box
    Vec3 boxTrans = Ortho(Coord, bp, bm, frameIn.BoxCrd()) + offset;

    // Translate atoms according to Coord
    AtomPairs.DoTranslation(frameIn, idx, boxTrans);
  } // END loop over atom pairs
}

// Image::Ortho()
/** \param Coord Coordinate to image
  * \param bp Box + boundary
  * \param bm Box - boundary
  * \param BoxVec box lengths.
  * \return Vector containing image translation
  */
Vec3 Image::Ortho(Vec3 const& Coord, Vec3 const& bp, Vec3 const& bm, Box const& BoxVec)
{
  Vec3 trans;
  // Determine how far Coord is out of box
  // Note Box::Param 0 1 2 is X Y Z
  for (int i = 0; i < 3; ++i) {
    trans[i] = 0.0;
    double crd = Coord[i];
    while (crd < bm[i]) {
      crd += BoxVec.Param((Box::ParamType)i);
      trans[i] += BoxVec.Param((Box::ParamType)i);
    }
    while (crd > bp[i]) {
      crd -= BoxVec.Param((Box::ParamType)i);
      trans[i] -= BoxVec.Param((Box::ParamType)i);
    }
  }
  return trans;
}

// -----------------------------------------------------------------------------
// Image::UnwrapNonortho()
void Image::UnwrapNonortho(Frame& tgtIn, Frame& refIn, List const& AtomPairs,
                           Unit const& allEntities,
                           Matrix_3x3 const& ucell, Matrix_3x3 const& recip)
{
  Vec3 vtgt, vref, boxTrans;
  // Loop over atom pairs
  for (unsigned int idx = 0; idx != AtomPairs.nEntities(); idx++)
  {
    vtgt = AtomPairs.GetCoord(idx, tgtIn);
    vref = AtomPairs.GetCoord(idx, refIn);

    boxTrans.Zero();
    // Calculate original distance from the ref (previous) position. 
    Vec3 vd = vtgt - vref; // dx dy dz
    double minDistanceSquare = vd.Magnitude2();
    // Reciprocal coordinates
    vd = recip * vd ; // recip * dxyz
    double cx = floor(vd[0]);
    double cy = floor(vd[1]);
    double cz = floor(vd[2]);
    // Loop over all possible translations 
    for (int ix = -1; ix < 2; ++ix) {
      for (int iy = -1; iy < 2; ++iy) {
        for (int iz = -1; iz < 2; ++iz) {
          // Calculate the translation.
          Vec3 vcc = ucell.TransposeMult( Vec3( cx+(double)ix, 
                                                cy+(double)iy, 
                                                cz+(double)iz ) ); // ucell^T * ccxyz
          // Calc. the potential new coordinate for tgt
          Vec3 vnew = vtgt - vcc; 
          // Calc. the new distance from the ref (previous) position
          Vec3 vr = vref - vnew; 
          double distanceSquare = vr.Magnitude2();
          // If the orig. distance is greater than the new distance, unwrap. 
          if ( minDistanceSquare > distanceSquare ) {
              minDistanceSquare = distanceSquare;
              boxTrans = vcc;
          }
        }
      }
    }
    // Translate tgt atoms
    boxTrans.Neg();
    AtomPairs.DoTranslation( tgtIn, idx, boxTrans );
  } // END loop over atom pairs 
  // Save new ref positions
  refIn.CopyFrom(tgtIn, allEntities);
}

// Image::UnwrapOrtho()
void Image::UnwrapOrtho(Frame& tgtIn, Frame& refIn, List const& AtomPairs,
                        Unit const& allEntities)
{
  Vec3 vtgt, vref, boxTrans;
  Vec3 boxVec = tgtIn.BoxCrd().Lengths();
  // Loop over atom pairs
  for (unsigned int idx = 0; idx != AtomPairs.nEntities(); idx++)
  {
    vtgt = AtomPairs.GetCoord(idx, tgtIn);
    vref = AtomPairs.GetCoord(idx, refIn);

    Vec3 dxyz = vtgt - vref;
    boxTrans[0] = -floor( dxyz[0] / boxVec[0] + 0.5 ) * boxVec[0];
    boxTrans[1] = -floor( dxyz[1] / boxVec[1] + 0.5 ) * boxVec[1];
    boxTrans[2] = -floor( dxyz[2] / boxVec[2] + 0.5 ) * boxVec[2];
    // Translate atoms from first to last
    AtomPairs.DoTranslation( tgtIn, idx, boxTrans );
  } // END loop over atom pairs
  // Save new ref positions
  refIn.CopyFrom(tgtIn, allEntities);
}

// -----------------------------------------------------------------------------
void Image::WrapToCell0(std::vector<double>& CoordsIn, Frame const& frmIn,
                        AtomMask const& maskIn,
                        Matrix_3x3 const& ucell, Matrix_3x3 const& recip)
{
  double* uFrac = &CoordsIn[0];
  int nUatoms = maskIn.Nselected();
  int idx;
  double* result;
  const double* XYZ;
# ifdef _OPENMP
# pragma omp parallel private(idx, result, XYZ)
  {
# pragma omp for
# endif
  for (idx = 0; idx < nUatoms; idx++)
  {
    result = uFrac + idx*3;
    XYZ = frmIn.XYZ( maskIn[idx] );
    // Convert to frac coords
    recip.TimesVec( result, XYZ );
    // Wrap to primary unit cell
    result[0] = result[0] - floor(result[0]);
    result[1] = result[1] - floor(result[1]);
    result[2] = result[2] - floor(result[2]);
    // Convert back to Cartesian
    ucell.TransposeMult( result, result );
  }
# ifdef _OPENMP
  } // END pragma omp parallel
# endif
}
