#include <cmath> // floor
#include "ImageRoutines.h"
#include "DistRoutines.h"

// SetupImageTruncoct()
/** Set up centering if putting nonortho cell into familiar trunc. oct. shape.
  * \param frameIn Frame to set up for.
  * \param ComMask If not NULL center is calcd w.r.t. center of atoms in mask.
  * \param useMass If true calculate COM, otherwise calc geometric center.
  * \param origin If true and ComMask is NULL use origin, otherwise use box center.
  * \return Coordinates of center.
  */
Vec3 SetupImageTruncoct( Frame& frameIn, AtomMask* ComMask, bool useMass, bool origin)
{
  if (ComMask!=NULL) {
    // Use center of atoms in mask
    if (useMass)
      return frameIn.VCenterOfMass( *ComMask );
    else
      return frameIn.VGeometricCenter( *ComMask );
  } else if (!origin) {
    // Use box center
    return Vec3( frameIn.BoxX() / 2, frameIn.BoxY() / 2, frameIn.BoxZ() / 2 );
  }
  //fprintf(stdout,"DEBUG: fcom = %lf %lf %lf\n",fcom[0],fcom[1],fcom[2]);
  return Vec3(); // Default is origin {0,0,0}
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
                   Matrix_3x3 ucell, Matrix_3x3 recip, // TODO: Make const &
                   bool truncoct, bool center,
                   bool useMass, std::vector<int> const& AtomPairs)
{
  Vec3 Coord;
  double min = -1.0;

  if (truncoct)
    min = 100.0 * (frameIn.BoxX()*frameIn.BoxX()+
                   frameIn.BoxY()*frameIn.BoxY()+
                   frameIn.BoxZ()*frameIn.BoxZ());

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

    frameIn.Translate(boxTrans.Dptr(), firstAtom, lastAtom);

  } // END loop over atom pairs
}

// ImageNonortho()
Vec3 ImageNonortho(Vec3 const& Coord, bool truncoct, 
                   bool origin, const Matrix_3x3& ucell, const Matrix_3x3& recip, 
                   Vec3 const& fcom, double min)
{
  int ixyz[3];

  Vec3 fc = recip * Coord;

  if ( origin )
    fc += 0.5; 

  Vec3 boxTransOut = ucell.TransposeMult( Vec3( floor(fc[0]), floor(fc[1]), floor(fc[2])  ) );
  boxTransOut.Neg();

  // Put into familiar trunc. oct. shape
  if (truncoct) {
    Vec3 TransCoord = Coord;
    TransCoord += boxTransOut;
    // ----------
    //MinImageNonOrtho2(Coord, fcom, box_, (int)origin, ixyz, ucell, recip);
    Vec3 f  = recip * TransCoord;
    Vec3 f2 = recip * fcom;

    if (origin) {
      f += 0.5;
      f2 += 0.5;
    }

    min = DIST2_ImageNonOrthoRecip(f.Dptr(), f2.Dptr(), min, ixyz, ucell.Dptr());
    // ----------
    if (ixyz[0] != 0 || ixyz[1] != 0 || ixyz[2] != 0) {
      boxTransOut += ucell.TransposeMult( ixyz );

      //if (debug > 2)
      //  mprintf( "  IMAGING, FAMILIAR OFFSETS ARE %i %i %i\n", 
      //          ixyz[0], ixyz[1], ixyz[2]);
    }
  }
  return boxTransOut;
}

void SetupImageOrtho(Frame& frameIn, Vec3& bp, Vec3& bm, bool origin) {
  // Set up boundary information for orthorhombic cell
  if (origin) {
    bp.SetVec( frameIn.BoxX() / 2,
               frameIn.BoxY() / 2,
               frameIn.BoxZ() / 2 );
    bm.SetVec( -bp[0], -bp[1], -bp[2] );
  } else {
    bp.SetVec( frameIn.BoxX(),
               frameIn.BoxY(),
               frameIn.BoxZ()  );
    bm.Zero();
  }
}

void ImageOrtho(Frame& frameIn, Vec3 const& bp, Vec3 const& bm, bool center, bool useMass,
                std::vector<int> const& AtomPairs)
{
  Vec3 Coord;
  Vec3 BoxVec( frameIn.BoxX(), frameIn.BoxY(), frameIn.BoxZ() );

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
    Vec3 boxTrans = ImageOrtho(Coord, bp, bm, BoxVec);

    // Translate atoms according to Coord
    frameIn.Translate(boxTrans.Dptr(), firstAtom, lastAtom);
  } // END loop over atom pairs
}

// Frame::ImageOrtho()
Vec3 ImageOrtho(Vec3 const& Coord, Vec3 const& bp, Vec3 const& bm, Vec3 const& BoxVec)
{
  double trans[3];
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
  return Vec3( trans );
}

