#include <cmath> // floor
#include "ImageRoutines.h"
#include "DistRoutines.h"

// SetupImageTruncoct()
/** Set up centering if putting nonortho cell into familiar trunc. oct. shape.
  * \param frameIn Frame to set up for.
  * \param ComMask If not null center is calcd w.r.t. center of atoms in mask.
  * \param useMass If true calculate COM, otherwise calc geometric center.
  * \param origin If true and ComMask is null use origin, otherwise use box center.
  * \return Coordinates of center.
  */
Vec3 SetupImageTruncoct( Frame const& frameIn, AtomMask* ComMask, bool useMass, bool origin)
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

// ImageNonortho()
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
void ImageNonortho(Frame& frameIn, bool origin, Vec3 const& fcom, 
                   Matrix_3x3 const& ucell, Matrix_3x3 const& recip,
                   bool truncoct, bool center,
                   bool useMass, std::vector<int> const& AtomPairs)
{
  Vec3 Coord;
  double min = -1.0;

  if (truncoct)
    min = 100.0 * (frameIn.BoxCrd().BoxX()*frameIn.BoxCrd().BoxX()+
                   frameIn.BoxCrd().BoxY()*frameIn.BoxCrd().BoxY()+
                   frameIn.BoxCrd().BoxZ()*frameIn.BoxCrd().BoxZ());

  // Loop over atom pairs
  for (std::vector<int>::const_iterator atom = AtomPairs.begin();
                                        atom != AtomPairs.end(); ++atom)
  {
    int firstAtom = *atom;
    ++atom;
    int lastAtom = *atom;
    //if (debug>2)
    //  mprintf( "  IMAGE processing atoms %i to %i\n", firstAtom+1, lastAtom);
    // Set up Coord with position to check for imaging based on first atom or 
    // center of mass of atoms first to last.
    if (center) {
      if (useMass)
        Coord = frameIn.VCenterOfMass(firstAtom,lastAtom);
      else
        Coord = frameIn.VGeometricCenter(firstAtom,lastAtom);
    } else 
      Coord = frameIn.XYZ( firstAtom );

    // boxTrans will hold calculated translation needed to move atoms back into box
    Vec3 boxTrans = ImageNonortho(Coord, truncoct, origin, ucell, recip, fcom, min);

    frameIn.Translate(boxTrans, firstAtom, lastAtom);

  } // END loop over atom pairs
}

// ImageNonortho()
/** \param Coord Coordinate to image.
  * \param truncoct If true, image in truncated octahedral shape.
  * \param origin If true, image w.r.t. coordinate origin.
  * \param ucell Unit cell matrix.
  * \param recip Reciprocal coordinates matrix.
  * \param fcom If truncoct, image translated coordinate w.r.t. this coord.
  * \return Vector containing image translation.
  */
Vec3 ImageNonortho(Vec3 const& Coord, bool truncoct, 
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

// SetupImageOrtho()
/** \param boxIn Box coordinates of Frame to image.
  * \param bp Output: Box + boundary.
  * \param bm Output: Box - boundary.
  * \param origin If true, image w.r.t. coordinate origin, otherwise box center.
  * \return 1 if box lengths are zero, 0 if setup completed successfully.
  */
int SetupImageOrtho(Box const& boxIn, Vec3& bp, Vec3& bm, bool origin) {
  // Set up boundary information for orthorhombic cell
  if (origin) {
    bp = boxIn.Center();
    bm.SetVec( -bp[0], -bp[1], -bp[2] );
  } else {
    bp.SetVec( boxIn.BoxX(), boxIn.BoxY(), boxIn.BoxZ()  );
    bm.Zero();
  }
  if (bp.IsZero()) return 1;
  return 0;
}

// ImageOrtho()
/** \param frameIn Frame to image.
  * \param bp Box + boundary.
  * \param bm Box - boundary.
  * \param center If true image w.r.t. center of atoms, otherwise first atom.
  * \param useMass If true calc center of mass, otherwise geometric center.
  */
void ImageOrtho(Frame& frameIn, Vec3 const& bp, Vec3 const& bm, bool center, bool useMass,
                std::vector<int> const& AtomPairs)
{
  Vec3 Coord;

  // Loop over atom pairs
  for (std::vector<int>::const_iterator atom = AtomPairs.begin();
                                        atom != AtomPairs.end(); atom++)
  {
    int firstAtom = *atom;
    ++atom;
    int lastAtom = *atom;
    //if (debug>2)
    //  mprintf( "  IMAGE processing atoms %i to %i\n", firstAtom+1, lastAtom);
    // Set up Coord with position to check for imaging based on first atom or 
    // center of mass of atoms first to last.
    if (center) {
      if (useMass)
        Coord = frameIn.VCenterOfMass(firstAtom,lastAtom);
      else
        Coord = frameIn.VGeometricCenter(firstAtom,lastAtom);
    } else 
      Coord = frameIn.XYZ( firstAtom );

    // boxTrans will hold calculated translation needed to move atoms back into box
    Vec3 boxTrans = ImageOrtho(Coord, bp, bm, frameIn.BoxCrd());

    // Translate atoms according to Coord
    frameIn.Translate(boxTrans, firstAtom, lastAtom);
  } // END loop over atom pairs
}

// ImageOrtho()
/** \param Coord Coordinate to image
  * \param bp Box + boundary
  * \param bm Box - boundary
  * \param BoxVec box lengths.
  * \return Vector containing image translation
  */
Vec3 ImageOrtho(Vec3 const& Coord, Vec3 const& bp, Vec3 const& bm, Box const& BoxVec)
{
  Vec3 trans;
  // Determine how far Coord is out of box
  for (int i = 0; i < 3; ++i) {
    trans[i] = 0.0;
    double crd = Coord[i];
    while (crd < bm[i]) {
      crd += BoxVec[i];
      trans[i] += BoxVec[i];
    }
    while (crd > bp[i]) {
      crd -= BoxVec[i];
      trans[i] -= BoxVec[i];
    }
  }
  return trans;
}

// UnwrapNonortho()
void UnwrapNonortho( Frame& frameIn, Frame& ref, AtomMask const& mask,
                     Matrix_3x3 const& ucell, Matrix_3x3 const& recip ) 
{
  for (AtomMask::const_iterator atom = mask.begin();
                                atom != mask.end(); ++atom)
  {
    int i3 = *atom * 3;
    Vec3 vtgt = frameIn.CRD( i3 );
    double minX = vtgt[0];
    double minY = vtgt[1];
    double minZ = vtgt[2];
    Vec3 vref = ref.CRD( i3 );

    Vec3 vd = vtgt - vref; // dx dy dz
    double minDistanceSquare = vd.Magnitude2();

    vd = recip * vd ; // recip * dxyz

    double cx = floor(vd[0]);
    double cy = floor(vd[1]);
    double cz = floor(vd[2]);
    
    for (int ix = -1; ix < 2; ++ix) {
      for (int iy = -1; iy < 2; ++iy) {
        for (int iz = -1; iz < 2; ++iz) {
          Vec3 vcc = ucell.TransposeMult( Vec3( cx+(double)ix, 
                                                cy+(double)iy, 
                                                cz+(double)iz ) ); // ucell^T * ccxyz

          Vec3 vnew = vtgt - vcc; 
 
          Vec3 vr = vref - vnew; 
  
          double distanceSquare = vr.Magnitude2();

          if ( minDistanceSquare > distanceSquare ) {
              minDistanceSquare = distanceSquare;
              minX = vnew[0];
              minY = vnew[1];
              minZ = vnew[2];
          }
        }
      }
    }
    ref[i3  ] = frameIn[i3  ] = minX; 
    ref[i3+1] = frameIn[i3+1] = minY;
    ref[i3+2] = frameIn[i3+2] = minZ;

  } // END loop over selected atoms
}

// UnwrapOrtho()
void UnwrapOrtho( Frame& frameIn, Frame& ref, AtomMask const& mask ) {
  double boxX = frameIn.BoxCrd().BoxX();
  double boxY = frameIn.BoxCrd().BoxY();
  double boxZ = frameIn.BoxCrd().BoxZ();
  for (AtomMask::const_iterator atom = mask.begin();
                                atom != mask.end(); ++atom)
  {
    int i3 = *atom * 3;
    double dx = frameIn[i3  ] - ref[i3  ];
    double dy = frameIn[i3+1] - ref[i3+1];
    double dz = frameIn[i3+2] - ref[i3+2];

    ref[i3  ] = frameIn[i3  ] = frameIn[i3  ] - floor( dx / boxX + 0.5 ) * boxX;
    ref[i3+1] = frameIn[i3+1] = frameIn[i3+1] - floor( dy / boxY + 0.5 ) * boxY;
    ref[i3+2] = frameIn[i3+2] = frameIn[i3+2] - floor( dz / boxZ + 0.5 ) * boxZ;
  }
}
