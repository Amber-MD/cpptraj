#include <cmath> // sqrt
#include <cstring> // memcpy, memset
#include "Frame.h"
#include "Constants.h" // SMALL
#include "vectormath.h" // CROSS_PRODUCT, normalize
#include "CpptrajStdio.h"

const size_t Frame::COORDSIZE_ = 3 * sizeof(double);
const size_t Frame::BOXSIZE_ = 6 * sizeof(double);

// ---------- CONSTRUCTION/DESTRUCTION/ASSIGNMENT ------------------------------
/// CONSTRUCTOR
Frame::Frame( ) : 
  natom_(0),
  maxnatom_(0),
  ncoord_(0),
  T_(0.0),
  X_(NULL),
  V_(NULL)
{
  memset(box_, 0, BOXSIZE_);
}

/// DESTRUCTOR
// Defined since this class can be inherited.
Frame::~Frame( ) { 
  if (X_ != NULL) delete[] X_;
  if (V_ != NULL) delete[] V_;
}

// CONSTRUCTOR
/// Set up for natom
Frame::Frame(int natomIn) :
  natom_(natomIn),
  maxnatom_(natomIn),
  ncoord_(natomIn*3), 
  T_(0.0),
  X_(NULL),
  V_(NULL)
{
  memset(box_, 0, BOXSIZE_);
  if (ncoord_ > 0)
    X_ = new double[ ncoord_ ];
}

// CONSTRUCTOR
/// Set up for the given atom array (including mass).
Frame::Frame(std::vector<Atom> const& atoms) :
  natom_(atoms.size()),
  maxnatom_(natom_),
  ncoord_(natom_*3),
  T_(0.0),
  X_(NULL),
  V_(NULL)
{
  memset(box_, 0, BOXSIZE_);
  if (ncoord_ > 0) {
    X_ = new double[ ncoord_ ];
    Mass_.reserve( natom_ );
    for (std::vector<Atom>::const_iterator atom = atoms.begin(); atom != atoms.end(); ++atom)
      Mass_.push_back( (*atom).Mass() );
  }
}

// CONSTRUCTOR
/// Copy frameIn according to maskIn
Frame::Frame(Frame const& frameIn, AtomMask const& maskIn) : 
  natom_( maskIn.Nselected() ),
  maxnatom_(natom_),
  ncoord_(natom_*3),
  T_( frameIn.T_ ),
  X_(NULL),
  V_(NULL)
{
  if (ncoord_ > 0) {
    X_ = new double[ ncoord_ ];
    double* newX = X_;
    if ( !frameIn.Mass_.empty() ) {
      if ( frameIn.V_ != NULL ) {
        // Copy coords/mass/velo
        V_ = new double[ ncoord_ ];
        double* newV = V_;
        for (AtomMask::const_iterator atom = maskIn.begin(); atom != maskIn.end(); ++atom)
        {
          int oldcrd = ((*atom) * 3);
          memcpy(newX, frameIn.X_ + oldcrd, COORDSIZE_);
          newX += 3;
          memcpy(newV, frameIn.V_ + oldcrd, COORDSIZE_);
          newV += 3;
          Mass_.push_back( frameIn.Mass_[*atom] );
        }
      } else {
        // Copy coords/mass
        for (AtomMask::const_iterator atom = maskIn.begin(); atom != maskIn.end(); ++atom)
        {
          memcpy(newX, frameIn.X_ + ((*atom) * 3), COORDSIZE_);
          newX += 3;
          Mass_.push_back( frameIn.Mass_[*atom] );
        }
      }
    } else {
      // Copy coords only
      for (AtomMask::const_iterator atom = maskIn.begin(); atom != maskIn.end(); ++atom)
      {
        memcpy(newX, frameIn.X_ + ((*atom) * 3), COORDSIZE_);
        newX += 3;
      }
    }
  }
  // Copy box coords
  memcpy(box_, frameIn.box_, BOXSIZE_);
}

// COPY CONSTRUCTOR
Frame::Frame(const Frame& rhs) :
  natom_(rhs.natom_),
  maxnatom_(rhs.maxnatom_),
  ncoord_(rhs.ncoord_),
  T_(rhs.T_),
  X_(NULL),
  V_(NULL),
  Mass_(rhs.Mass_) 
{
  // Copy box coords
  memcpy(box_, rhs.box_, BOXSIZE_);
  // Copy coords/velo; allocate for maxnatom but copy natom
  int maxncoord = maxnatom_ * 3;
  if (rhs.X_!=NULL) {
    X_ = new double[ maxncoord ];
    memcpy(X_, rhs.X_, natom_ * COORDSIZE_);
  }
  if (rhs.V_!=NULL) {
    V_ = new double[ maxncoord ];
    memcpy(V_, rhs.V_, natom_ * COORDSIZE_);
  }
}

// Frame::swap()
void Frame::swap(Frame &first, Frame &second) {
  using std::swap;
  swap(first.natom_, second.natom_);
  swap(first.maxnatom_, second.maxnatom_);
  swap(first.ncoord_, second.ncoord_);
  swap(first.T_, second.T_);
  swap(first.X_, second.X_);
  swap(first.V_, second.V_);
  first.Mass_.swap(second.Mass_);
  swap( first.box_[0], second.box_[0] );
  swap( first.box_[1], second.box_[1] );
  swap( first.box_[2], second.box_[2] );
  swap( first.box_[3], second.box_[3] );
  swap( first.box_[4], second.box_[4] );
  swap( first.box_[5], second.box_[5] );
}

// Frame::operator=()
/** Assignment operator using copy/swap idiom. */
Frame &Frame::operator=(Frame rhs) {
  swap(*this, rhs);
  return *this;
}

// ---------- CONVERT TO/FROM ARRAYS -------------------------------------------
/** Assign float array to this frame. */ 
Frame& Frame::operator=(std::vector<float> const& farray) {
  int f_ncoord = (int)farray.size();
  if (f_ncoord > maxnatom_*3) {
    mprinterr("Error: Float array size %zu > max #coords in frame %zu\n",
              farray.size(), maxnatom_*3);
    return *this;
  }
  double* Xptr = X_;
  for (std::vector<float>::const_iterator Fptr = farray.begin(); 
                                          Fptr != farray.end(); ++Fptr)
  {
    *Xptr = (double)(*Fptr);
    ++Xptr;
  }
  ncoord_ = f_ncoord;
  natom_ = ncoord_ / 3;
  return *this;
}

// Frame::ConvertToFloat()
/** Place atom coordinates according to maskIn into a float array. */
std::vector<float> Frame::ConvertToFloat(AtomMask const& maskIn) {
  std::vector<float> farray;

  farray.reserve( maskIn.Nselected() * 3 );
  for (AtomMask::const_iterator atom = maskIn.begin(); atom != maskIn.end(); ++atom)
  {
    int crdidx = (*atom) * 3; 
    farray.push_back( (float)X_[crdidx  ] );
    farray.push_back( (float)X_[crdidx+1] );
    farray.push_back( (float)X_[crdidx+2] );
  }
  return farray;
}

// ---------- ACCESS INTERNAL DATA ---------------------------------------------
// Frame::printAtomCoord()
/** Print XYZ coords of given atom */
void Frame::printAtomCoord(int atom) {
  int atmidx = atom * 3;
  if (atmidx >= ncoord_) return;
  mprintf("ATOM %i: %lf %lf %lf\n",atom,X_[atmidx],X_[atmidx+1],X_[atmidx+2]);
}

// Frame::Info()
// For debugging
void Frame::Info(const char *msg) {
  if (msg!=NULL)
    mprintf("\tFrame [%s]:",msg);
  else
    mprintf("\tFrame:");
  mprintf("%i atoms, %i coords",natom_, ncoord_);
  if (V_!=NULL) mprintf(" with Velocities");
  if (!Mass_.empty()) mprintf(" with Masses");
  mprintf("\n");
}

// Frame::GetAtomXYZ()
void Frame::GetAtomXYZ(double *Coord, int atom) {
  double* Xptr = X_ + (atom * 3);
  Coord[0] = *(Xptr++);
  Coord[1] = *(Xptr++);
  Coord[2] = *(Xptr);
}

// Frame::AddXYZ()
/** Append the given XYZ coord to this frame. */
void Frame::AddXYZ(const double *XYZin) {
  if (XYZin == NULL) return;
  if (natom_ >= maxnatom_) {
    // Reallocate
    maxnatom_ += 500;
    double *newX = new double[ maxnatom_ * 3 ];
    if (X_!=NULL) {
      memcpy(newX, X_, natom_ * COORDSIZE_);
      delete[] X_;
    }
    X_ = newX;
  }
  memcpy(X_ + ncoord_, XYZin, COORDSIZE_);
  ++natom_;
  ncoord_ += 3;
}

// ---------- FRAME MEMORY ALLOCATION/REALLOCATION -----------------------------
// Frame::SetupFrame()
/** Set up frame for given number of atoms, no mass or velocity information.
  */
int Frame::SetupFrame(int natomIn) {
  natom_ = natomIn;
  ncoord_ = natom_ * 3;
  if (natom_ > maxnatom_) {
    // Reallocate
    if (X_ != NULL) delete[] X_;
    X_ = new double[ ncoord_ ];
    maxnatom_ = natom_;
  }
  if (V_ != NULL) delete[] V_;
  Mass_.clear();
  return 0;
}

// Frame::SetupFrameM()
/** Set up frame for given atom array (mass info included), no velocities. */
int Frame::SetupFrameM(std::vector<Atom> const& atoms) {
  bool reallocate = false;
  natom_ = (int)atoms.size();
  ncoord_ = natom_ * 3;
  // First check if reallocation must occur. If so reallocate coords. Masses
  // will be reallocated as well, or allocated if not already.
  if (natom_ > maxnatom_) {
    reallocate = true;
    if (X_ != NULL) delete[] X_;
    X_ = new double[ ncoord_ ];
    maxnatom_ = natom_;
  }
  if (reallocate || Mass_.empty())
    Mass_.resize(maxnatom_);
  // Copy masses
  Darray::iterator mass = Mass_.begin();
  for (std::vector<Atom>::const_iterator atom = atoms.begin();
                                         atom != atoms.end(); ++atom)
    *(mass++) = (*atom).Mass();
  if (V_ != NULL) delete[] V_;
  return 0;
}

// Frame::SetupFrameV()
/** Set up frame for given atom array. Set up for velocity info if 
  * hasVelocity is true.
  */
int Frame::SetupFrameV(std::vector<Atom> const& atoms, bool hasVelocity) {
  bool reallocate = false;
  natom_ = (int)atoms.size();
  ncoord_ = natom_ * 3;
  if (natom_ > maxnatom_) {
    reallocate = true;
    if (X_ != NULL) delete[] X_;
    X_ = new double[ ncoord_ ];
    maxnatom_ = natom_;
  }
  if (hasVelocity) {
    if (reallocate || V_ == NULL) {
      if (V_ != NULL) delete[] V_;
      V_ = new double[ maxnatom_*3 ];
      // Since velocity might not be read in, initialize it to 0.
      memset(V_, 0, maxnatom_ * COORDSIZE_);
    }
  } else {
    if (V_ != NULL) delete[] V_;
  }
  if (reallocate || Mass_.empty())
    Mass_.resize(maxnatom_);
  // Copy masses
  Darray::iterator mass = Mass_.begin();
  for (std::vector<Atom>::const_iterator atom = atoms.begin();
                                         atom != atoms.end(); ++atom)
    *(mass++) = (*atom).Mass();
  return 0;
}

// Frame::SetupFrameFromMask()
/** Set up frame to hold # selected atoms in given mask. If mass 
  * information is passed in store the masses corresponding to
  * selected atoms in the mask. No velocity info. Only reallocate
  * memory if Nselected > maxnatom.
  */
int Frame::SetupFrameFromMask(AtomMask const& maskIn, std::vector<Atom> const& atoms) {
  bool reallocate = false;
  natom_ = maskIn.Nselected();
  ncoord_ = natom_ * 3;
  if (natom_ > maxnatom_) {
    reallocate = true;
    if (X_ != NULL) delete[] X_;
    X_ = new double[ ncoord_ ];
    maxnatom_ = natom_;
  }
  if (reallocate || Mass_.empty())
    Mass_.resize(maxnatom_);
  // Copy masses according to maskIn
  Darray::iterator mass = Mass_.begin();
  for (AtomMask::const_iterator atom = maskIn.begin(); atom != maskIn.end(); ++atom) 
    *(mass++) = atoms[ *atom ].Mass();
  return 0; 
}

// ---------- FRAME SETUP OF COORDINATES ---------------------------------------
// Frame::SetCoordinatesByMask()
/** Copy raw double array into frame by mask. */
void Frame::SetCoordinatesByMask(const double* Xin, AtomMask const& maskIn) {
  // Copy only coords from Xin according to maskIn
  if (Xin==NULL) {
    mprinterr("Internal Error: SetCoordinatesByMask called with NULL pointer.\n");
    return;
  }
  if (maskIn.Nselected() > maxnatom_) {
    mprinterr("Internal Error: SetCoordinatesByMask: Selected=%i, maxnatom = %i\n",
              maskIn.Nselected(), maxnatom_);
    return;
  }
  natom_ = maskIn.Nselected();
  ncoord_ = natom_ * 3;
  double* newXptr = X_;
  for (AtomMask::const_iterator atom = maskIn.begin(); atom != maskIn.end(); ++atom)
  {
    memcpy( newXptr, Xin + ((*atom) * 3), COORDSIZE_);
    newXptr += 3;
  }
}

// Frame::SetCoordinates()
/** Copy only coordinates, box, and T from input frame to this frame based
  * on selected atoms in mask.
  */
void Frame::SetCoordinates(Frame const& frameIn, AtomMask const& maskIn) {
  if (maskIn.Nselected() > maxnatom_) {
    mprinterr("Error: SetCoordinates: Mask [%s] selected (%i) > max natom (%i)\n",
              maskIn.MaskString(), maskIn.Nselected(), maxnatom_);
    return;
  }
  natom_ = maskIn.Nselected();
  ncoord_ = natom_ * 3;
  memcpy(box_, frameIn.box_, BOXSIZE_);
  T_ = frameIn.T_;
  double* newXptr = X_;
  for (AtomMask::const_iterator atom = maskIn.begin(); atom != maskIn.end(); ++atom)
  {
    memcpy( newXptr, frameIn.X_ + ((*atom) * 3), COORDSIZE_);
    newXptr += 3;
  }
}

// Frame::SetCoordinates()
/** Copy only coordinates from input frame to this frame.
  */
void Frame::SetCoordinates(Frame const& frameIn) {
  if (frameIn.natom_ > maxnatom_) {
    mprinterr("Error: Frame::SetCoordinates: Input frame atoms (%i) > max natom (%i)\n",
              frameIn.natom_, maxnatom_);
    return;
  }
  natom_ = frameIn.natom_;
  ncoord_ = natom_ * 3;
  memcpy(X_, frameIn.X_, natom_ * COORDSIZE_);
}

// Frame::SetFrame()
/** Copy entire input frame (including velocity if defined in both) to this 
  * frame based on selected atoms in mask. Assumes coords and mass exist. 
  */
void Frame::SetFrame(Frame const& frameIn, AtomMask const& maskIn) {
  if (maskIn.Nselected() > maxnatom_) {
    mprinterr("Error: SetFrame: Mask [%s] selected (%i) > max natom (%i)\n",
              maskIn.MaskString(), maskIn.Nselected(), maxnatom_);
    return;
  }
  natom_ = maskIn.Nselected();
  ncoord_ = natom_ * 3;
  // Copy T/box
  memcpy(box_, frameIn.box_, BOXSIZE_);
  T_ = frameIn.T_;
  double* newXptr = X_;
  Darray::iterator mass = Mass_.begin();
  if (frameIn.V_ != NULL && V_ != NULL) {
    // Copy Coords/Mass/Velo
    double *newVptr = V_;
    for (AtomMask::const_iterator atom = maskIn.begin(); atom != maskIn.end(); ++atom)
    {
      int oldcrd = ((*atom) * 3);
      memcpy( newXptr, frameIn.X_ + oldcrd, COORDSIZE_);
      newXptr += 3;
      memcpy( newVptr, frameIn.V_ + oldcrd, COORDSIZE_);
      newVptr += 3;
      *mass = frameIn.Mass_[*atom];
      ++mass;
    }
  } else {
    // Copy coords/mass only
    for (AtomMask::const_iterator atom = maskIn.begin(); atom != maskIn.end(); ++atom)
    {
      memcpy( newXptr, frameIn.X_ + ((*atom) * 3), COORDSIZE_);
      newXptr += 3;
      *mass = frameIn.Mass_[*atom];
      ++mass;
    }
  }
}

// Frame::FrameCopy()
/** Return a copy of the frame */
// TODO: Obsolete
Frame *Frame::FrameCopy( ) {
  Frame* newFrame = new Frame( *this );
  return newFrame;
}

// ---------- FRAME SETUP WITH ATOM MAPPING ------------------------------------
// NOTE: Should these map routines be part of a child class used only by AtomMap?
// Frame::SetCoordinatesByMap()
/** Reorder atoms in the given target frame so that it matches the given
  * atom map with format (tgtAtom = Map[refAtom]). End result is that
  * this frame will be in the same order as the original reference. Assumes
  * that the map is complete (i.e. no unmapped atoms).
  */
void Frame::SetCoordinatesByMap(Frame const& tgtIn, std::vector<int> const& mapIn) {
  if (tgtIn.natom_ > maxnatom_) {
    mprinterr("Error: SetCoordinatesByMap: # Input map frame atoms (%i) > max atoms (%i)\n",
              tgtIn.natom_, maxnatom_);
    return;
  }
  if ((int)mapIn.size() != tgtIn.natom_) {
    mprinterr("Error: SetCoordinatesByMap: Input map size (%zu) != input frame natom (%i)\n",
              mapIn.size(), tgtIn.natom_);
    return;
  }
  natom_ = tgtIn.natom_;
  ncoord_ = natom_ * 3;
  double* newXptr = X_;
  for (std::vector<int>::const_iterator refatom = mapIn.begin(); 
                                        refatom != mapIn.end(); ++refatom)
  {
    memcpy( newXptr, tgtIn.X_ + ((*refatom) * 3), COORDSIZE_ );
    newXptr += 3;
  }
}

// Frame::SetReferenceByMap()
/** Set this frame to include only atoms from the given reference that are 
  * mapped.
  */
void Frame::SetReferenceByMap(Frame const& refIn, std::vector<int> const& mapIn) {
  if (refIn.natom_ > maxnatom_) {
    mprinterr("Error: SetReferenceByMap: # Input map frame atoms (%i) > max atoms (%i)\n",
              refIn.natom_, maxnatom_);
    return;
  }
  if ((int)mapIn.size() != refIn.natom_) {
    mprinterr("Error: SetReferenceByMap: Input map size (%zu) != input frame natom (%i)\n",
              mapIn.size(), refIn.natom_);
    return;
  }
  double* newXptr = X_;
  double* refptr = refIn.X_;
  for (std::vector<int>::const_iterator refatom = mapIn.begin(); 
                                        refatom != mapIn.end(); ++refatom)
  {
    if (*refatom != -1) {
      memcpy( newXptr, refptr, COORDSIZE_ );
      newXptr += 3;
    }
    refptr += 3;
  }
  ncoord_ = (int)(newXptr - X_);
  natom_ = ncoord_ / 3;
}

// Frame::SetTargetByMap()
/** Set this frame to include only atoms from the given target frame, remapped
  * according to the given atom map.
  */
void Frame::SetTargetByMap(Frame const& tgtIn, std::vector<int> const& mapIn) {
  if (tgtIn.natom_ > maxnatom_) {
    mprinterr("Error: SetTargetByMap: # Input map frame atoms (%i) > max atoms (%i)\n",
              tgtIn.natom_, maxnatom_);
    return;
  }
  if ((int)mapIn.size() != tgtIn.natom_) {
    mprinterr("Error: SetTargetByMap: Input map size (%zu) != input frame natom (%i)\n",
              mapIn.size(), tgtIn.natom_);
    return;
  }
  double* newXptr = X_;
  for (std::vector<int>::const_iterator refatom = mapIn.begin(); 
                                        refatom != mapIn.end(); ++refatom)
  {
    if (*refatom != -1) {
      memcpy( newXptr, tgtIn.X_ + ((*refatom) * 3), COORDSIZE_ );
      newXptr += 3;
    }
  }
  ncoord_ = (int)(newXptr - X_);
  natom_ = ncoord_ / 3;
}

// ---------- BASIC ARITHMETIC ------------------------------------------------- 
// Frame::ZeroCoords()
/** Set all coords to 0.0 */
void Frame::ZeroCoords( ) {
  memset(X_, 0, natom_ * COORDSIZE_);
}

// Frame::operator+=()
Frame &Frame::operator+=(const Frame& rhs) {
  // For now ensure same natom
  if (natom_ != rhs.natom_) {
    mprinterr("Error: Frame::operator+=: Frames have different natom.\n");
    return *this;
  }
  for (int i = 0; i < ncoord_; ++i)
    X_[i] += rhs.X_[i];
  return *this;
}

// Frame::operator-=()
Frame &Frame::operator-=(const Frame& rhs) {
  // For now ensure same natom
  if (natom_ != rhs.natom_) {
    mprinterr("Error: Frame::operator-=: Frames have different natom.\n");
    return *this;
  }
  for (int i = 0; i < ncoord_; ++i)
    X_[i] -= rhs.X_[i];
  return *this;
}

// Frame::operator*=()
Frame &Frame::operator*=(const Frame& rhs) {
  // For now ensure same natom
  if (natom_ != rhs.natom_) {
    mprinterr("Error: Frame::operator*=: Frames have different natom.\n");
    return *this;
  }
  for (int i = 0; i < ncoord_; ++i)
    X_[i] *= rhs.X_[i];
  return *this;
}

// Frame::operator*()
// NOTE: return const instance and not const reference to disallow
//       nonsense statements like (a * b) = c
const Frame Frame::operator*(const Frame& rhs) const {
  return (Frame(*this) *= rhs);
}

// Frame::Divide()
/** Divide all coord values of this frame by divisor and store in dividend.
  */
int Frame::Divide(Frame const& dividend, double divisor) {
  if (divisor < SMALL) {
    mprinterr("Error: Frame::Divide(Frame,divisor): Detected divide by 0.\n");
    return 1;
  }
  // For now ensure same natom
  if (natom_ != dividend.natom_) {
    mprinterr("Error: Frame::Divide: Frames have different natom.\n");
    return 1;
  }
  for (int i=0; i < ncoord_; ++i)
    X_[i] = dividend.X_[i] / divisor;
  return 0;
}

// Frame::Divide()
/** Divide all coord values of this frame by divisor.
  */
void Frame::Divide(double divisor) {
  if (divisor < SMALL) {
    mprinterr("Error: Frame::Divide(divisor): Detected divide by 0.\n");
    return;
  }
  for (int i = 0; i < ncoord_; ++i)
    X_[i] /= divisor;
}

// Frame::AddByMask()
/** Increment atoms in this frame by the selected atoms in given frame.
  */
void Frame::AddByMask(Frame const& frameIn, AtomMask const& maskIn) {
  if (maskIn.Nselected() > maxnatom_) {
    mprinterr("Error: AddByMask: Input mask #atoms (%i) > frame #atoms (%i)\n",
              maskIn.Nselected(), maxnatom_);
    return;
  }
  unsigned int xidx = 0;
  for (AtomMask::const_iterator atom = maskIn.begin(); atom != maskIn.end(); ++atom)
  {
    unsigned int fidx = (*atom) * 3;
    X_[xidx  ] += frameIn.X_[fidx  ];
    X_[xidx+1] += frameIn.X_[fidx+1];
    X_[xidx+2] += frameIn.X_[fidx+2];
    xidx += 3;
  }
}

// ---------- CENTER OF MASS / GEOMETRIC CENTER --------------------------------
// Frame::CenterOfMass()
/** Given an AtomMask put center of mass of atoms in mask into Coord. Return 
  * sum of masses in Mask.
  */
double Frame::CenterOfMass(double* Coord, AtomMask const& Mask) 
{
  double Coord0 = 0.0;
  double Coord1 = 0.0;
  double Coord2 = 0.0;
  double sumMass = 0.0;

  for (AtomMask::const_iterator atom = Mask.begin(); atom != Mask.end(); ++atom)
  {
    unsigned int xidx = (*atom) * 3;
    double mass = Mass_[*atom];
    sumMass += mass;
    Coord0 += ( X_[xidx  ] * mass );
    Coord1 += ( X_[xidx+1] * mass );
    Coord2 += ( X_[xidx+2] * mass );
  }
  Coord[0] = Coord0;
  Coord[1] = Coord1;
  Coord[2] = Coord2;
  // NOTE: Not using == since it is unreliable for floating point numbers.
  // Should NEVER have a mass smaller than SMALL.
  if (sumMass < SMALL) return 0;
  Coord[0] /= sumMass;
  Coord[1] /= sumMass;
  Coord[2] /= sumMass;
  return sumMass;
}

// Frame::GeometricCenter()
/** Given an AtomMask put geometric center of atoms in mask into Coord. Return 
  * #atoms in Mask.
  */
double Frame::GeometricCenter(double* Coord, AtomMask const& Mask) {
  double Coord0 = 0.0;
  double Coord1 = 0.0;
  double Coord2 = 0.0;

  for (AtomMask::const_iterator atom = Mask.begin(); atom != Mask.end(); ++atom)
  {
    unsigned int xidx = (*atom) * 3;
    Coord0 += X_[xidx  ];
    Coord1 += X_[xidx+1];
    Coord2 += X_[xidx+2];
  }
  Coord[0] = Coord0;
  Coord[1] = Coord1;
  Coord[2] = Coord2;
  double sumMass = (double)Mask.Nselected();
  // NOTE: Not using == since it is unreliable for floating point numbers.
  // Should NEVER have a mass smaller than SMALL.
  if (sumMass < SMALL) return 0;
  Coord[0] /= sumMass;
  Coord[1] /= sumMass;
  Coord[2] /= sumMass;
  return sumMass;
}

// Frame::CenterOfMass()
/** Put center of mass of all atoms between start and stop in frame into 
  * Coord. Return sum of masses.
  */
double Frame::CenterOfMass(double *Coord, int startAtom, int stopAtom) 
{  
  double Coord0 = 0.0;
  double Coord1 = 0.0;
  double Coord2 = 0.0;
  double sumMass = 0.0;

  Darray::iterator mass = Mass_.begin() + startAtom;
  int startAtom3 = startAtom * 3;
  int stopAtom3 = stopAtom * 3;
  for (int i = startAtom3; i < stopAtom3; i += 3) {
    sumMass += (*mass);
    Coord0 += ( X_[i  ] * (*mass) );
    Coord1 += ( X_[i+1] * (*mass) );
    Coord2 += ( X_[i+2] * (*mass) );
    ++mass;
  }
  Coord[0] = Coord0;
  Coord[1] = Coord1;
  Coord[2] = Coord2;
   // NOTE: Not using == since it is unreliable for floating point numbers.
  // Should NEVER have a mass smaller than SMALL.
  if (sumMass < SMALL) return 0;
  Coord[0] /= sumMass;
  Coord[1] /= sumMass;
  Coord[2] /= sumMass;
  return sumMass;
}

// Frame::GeometricCenter()
/** Put geometric center of all atoms between start and stop in frame into 
  * Coord. Return #atoms used in calc.
  */
double Frame::GeometricCenter(double *Coord, int startAtom, int stopAtom) {  
  double Coord0 = 0.0;
  double Coord1 = 0.0;
  double Coord2 = 0.0;

  int startAtom3 = startAtom * 3;
  int stopAtom3 = stopAtom * 3;
  for (int i = startAtom3; i < stopAtom3; i += 3) {
    Coord0 += X_[i  ];
    Coord1 += X_[i+1];
    Coord2 += X_[i+2];
  }
  Coord[0] = Coord0;
  Coord[1] = Coord1;
  Coord[2] = Coord2;
  double sumMass = (double)(stopAtom - startAtom);
  // NOTE: Not using == since it is unreliable for floating point numbers.
  // Should NEVER have a mass smaller than SMALL.
  if (sumMass < SMALL) return 0;
  Coord[0] /= sumMass;
  Coord[1] /= sumMass;
  Coord[2] /= sumMass;
  return sumMass;
}

// ---------- COORDINATE MANIPULATION ------------------------------------------
// Frame::SCALE()
void Frame::SCALE(AtomMask const& maskIn, double sx, double sy, double sz) {
  for (AtomMask::const_iterator atom = maskIn.begin();
                                atom != maskIn.end(); atom++)
  {
    unsigned int xidx = *atom * 3;
    X_[xidx  ] *= sx;
    X_[xidx+1] *= sy;
    X_[xidx+2] *= sz;
  }
}

// Frame::Translate()
/** Translate all coords by Vec.  */
void Frame::Translate(const double* Vec) {
  double Vec0 = Vec[0];
  double Vec1 = Vec[1];
  double Vec2 = Vec[2];
  for (int i = 0; i < ncoord_; i += 3) {
    X_[i  ] += Vec0;
    X_[i+1] += Vec1;
    X_[i+2] += Vec2;
  }
}

// Frame::Translate()
/** Translate atoms in range by Vec. */
// NOTE: SHOULD CHECK BOUNDS! 
void Frame::Translate(const double *Vec, int firstAtom, int lastAtom) {
  double Vec0 = Vec[0];
  double Vec1 = Vec[1];
  double Vec2 = Vec[2];
  int startatom3 = firstAtom * 3;
  int stopatom3 = lastAtom * 3;
  for (int i = startatom3; i < stopatom3; i += 3) {
    X_[i  ] += Vec0;
    X_[i+1] += Vec1;
    X_[i+2] += Vec2;
  }
}

// Frame::Translate()
/** Translate atom by Vec */
void Frame::Translate(const double* Vec, int atom) {
  int xidx = atom * 3;
  X_[xidx  ] += Vec[0];
  X_[xidx+1] += Vec[1];
  X_[xidx+2] += Vec[2];
}

// Frame::Trans_Rot_Trans()
/** Given an array Vec of size 6 containing two translations:
  *   T0x T0y T0z T1x T1y T1z
  * and a rotation matrix T, apply the first translation, then
  * the rotation, then the second translation.
  */
void Frame::Trans_Rot_Trans(const double *Vec, const double *T) {
  double Vec0 = Vec[0];
  double Vec1 = Vec[1];
  double Vec2 = Vec[2];
  double Vec3 = Vec[3];
  double Vec4 = Vec[4];
  double Vec5 = Vec[5];

  double T0 = T[0]; 
  double T1 = T[1]; 
  double T2 = T[2]; 
  double T3 = T[3]; 
  double T4 = T[4]; 
  double T5 = T[5]; 
  double T6 = T[6]; 
  double T7 = T[7]; 
  double T8 = T[8]; 

  for (int i = 0; i < ncoord_; i += 3) {
    double *Xptr = X_ + i;
    double *Yptr = Xptr + 1;
    double *Zptr = Xptr + 2;
 
    double x = *Xptr + Vec0;
    double y = *Yptr + Vec1;
    double z = *Zptr + Vec2;

    *Xptr = (x*T0) + (y*T1) + (z*T2) + Vec3;
    *Yptr = (x*T3) + (y*T4) + (z*T5) + Vec4;
    *Zptr = (x*T6) + (y*T7) + (z*T8) + Vec5;
  }
}

// Frame::Rotate()
/** Multiply natomx3 matrix X by 3x3 matrix T. If T is a rotation matrix
  * this rotates the coords in X. 
  */
void Frame::Rotate(const double* T) {
  double T0 = T[0];
  double T1 = T[1];
  double T2 = T[2];
  double T3 = T[3];
  double T4 = T[4];
  double T5 = T[5];
  double T6 = T[6];
  double T7 = T[7];
  double T8 = T[8];

  for (int i = 0; i < ncoord_; i += 3) {
    double *Xptr = X_ + i;
    double *Yptr = Xptr + 1;
    double *Zptr = Xptr + 2;

    double x = *Xptr;
    double y = *Yptr;
    double z = *Zptr;

    *Xptr = (x*T0) + (y*T1) + (z*T2);
    *Yptr = (x*T3) + (y*T4) + (z*T5);
    *Zptr = (x*T6) + (y*T7) + (z*T8);
  }
} 

// Frame::InverseRotate()
/** Multiply natomx3 matrix X by transpose of 3x3 matrix T. If T is a rotation
  * matrix this rotates the coords in X in the opposite direction.
  */
void Frame::InverseRotate(const double* T) {
  double T0 = T[0];
  double T1 = T[1];
  double T2 = T[2];
  double T3 = T[3];
  double T4 = T[4];
  double T5 = T[5];
  double T6 = T[6];
  double T7 = T[7];
  double T8 = T[8];
  
  for (int i = 0; i < ncoord_; i += 3) {
    double *Xptr = X_ + i;
    double *Yptr = Xptr + 1;
    double *Zptr = Xptr + 2;

    double x = *Xptr;
    double y = *Yptr;
    double z = *Zptr;

    *Xptr = (x*T0) + (y*T3) + (z*T6);
    *Yptr = (x*T1) + (y*T4) + (z*T7);
    *Zptr = (x*T2) + (y*T5) + (z*T8);
  }
}

// Frame::Center()
/** Center coordinates to center of coordinates in Mask w.r.t. given XYZ in
  * boxcoord. When called from Action_Center boxcoord will be either origin 
  * or box center. Use geometric center if mass is NULL, otherwise center 
  * of mass will be used.
  */
void Frame::Center(AtomMask const& Mask, bool origin, bool useMassIn) 
{
  double center[3];

  if (useMassIn)
    this->CenterOfMass(center, Mask);
  else
    this->GeometricCenter(center, Mask);
  //mprinterr("  FRAME CENTER: %lf %lf %lf\n",center[0],center[1],center[2]); //DEBUG

  if (origin) {
    // Shift to coordinate origin (0,0,0)
    center[0] = -center[0];
    center[1] = -center[1];
    center[2] = -center[2];
  } else {
    // Shift to box center
    center[0] = (box_[0] / 2) - center[0];
    center[1] = (box_[1] / 2) - center[1];
    center[2] = (box_[2] / 2) - center[2];
  }

  this->Translate(center);
}

// Frame::CenterReference()
/** Center coordinates to origin in preparation for RMSD calculation. Store
  * translation vector from origin to reference in Trans.
  */
void Frame::CenterReference(double *Trans, bool useMassIn)
{
  double center[3];
  if (useMassIn)
    this->CenterOfMass(Trans,0,natom_);
  else
    this->GeometricCenter(Trans,0,natom_);
  //mprinterr("  REF FRAME CENTER: %lf %lf %lf\n",Trans[0],Trans[1],Trans[2]); //DEBUG
  // Trans now contains translation from origin -> Ref
  center[0] = -Trans[0];
  center[1] = -Trans[1];
  center[2] = -Trans[2];
  // Center now contains translation from Ref -> origin.
  this->Translate(center);
}

// Frame::ShiftToGeometricCenter()
/** Shift geometric center of coordinates in frame to origin. */
void Frame::ShiftToGeometricCenter( ) {
  double frameCOM[3];

  this->GeometricCenter(frameCOM, 0, natom_);
  //mprinterr("  FRAME COM: %lf %lf %lf\n",frameCOM[0],frameCOM[1],frameCOM[2]); //DEBUG
  
  // Shift to common COM
  frameCOM[0] = -frameCOM[0]; 
  frameCOM[1] = -frameCOM[1]; 
  frameCOM[2] = -frameCOM[2];
  this->Translate(frameCOM);
}

// ---------- COORDINATE CALCULATION ------------------------------------------- 
// Frame::BoxToRecip()
/** Use box coordinates to calculate reciprocal space conversions for use
  * with imaging routines. Return cell volume.
  */
// NOTE: Move to separate routine in DistRoutines?
double Frame::BoxToRecip(double *ucell, double *recip) {
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

// Frame::RADGYR()
/** Return the radius of gyration of atoms in mask. Also set the maximum 
  * distance from center. Use center of mass if useMassIn is true.
  */
double Frame::RADGYR(AtomMask& Mask, bool useMassIn, double *max) 
{
  double mid[3];
  double total_mass = 0.0;
  double maxMass = 1.0;
  double sumDist2 = 0.0;
  *max = 0.0;

  if (useMassIn) {
    total_mass = this->CenterOfMass(mid, Mask);
    for (AtomMask::const_iterator atom = Mask.begin(); atom != Mask.end(); ++atom)
    {
      int xidx = *atom * 3;
      double dx = X_[xidx  ] - mid[0];
      double dy = X_[xidx+1] - mid[1];
      double dz = X_[xidx+2] - mid[2];
      double dist2 = (dx*dx) + (dy*dy) + (dz*dz);
      double mass = Mass_[ *atom ];
      dist2 *= mass;
      if (dist2 > *max) {
        *max = dist2;
        maxMass = mass;
      }
      sumDist2 += dist2;
    }
  } else {
    total_mass = this->GeometricCenter(mid, Mask);
    for (AtomMask::const_iterator atom = Mask.begin(); atom != Mask.end(); ++atom)
    {
      int xidx = *atom * 3;
      double dx = X_[xidx  ] - mid[0];
      double dy = X_[xidx+1] - mid[1];
      double dz = X_[xidx+2] - mid[2];
      double dist2 = (dx*dx) + (dy*dy) + (dz*dz);
      if (dist2 > *max)
        *max = dist2;
      sumDist2 += dist2;
    }
  }
  // NOTE: Not using == since it is unreliable for floating point numbers.
  // Should NEVER have a mass smaller than SMALL (vectormath.h)
  if (total_mass < SMALL) return 0;

  *max = sqrt(*max / maxMass);

  return ( sqrt(sumDist2 / total_mass) ); // Radius of Gyration
}

// Frame::RMSD()
/** Get the RMSD of this Frame to Ref Frame. Ref frame must contain the same
  * number of atoms as this Frame - should be checked for before this routine
  * is called. Put the best-fit rotation matrix in U and the COM translation 
  * vectors in Trans. The translation is composed of two XYZ vectors; the first
  * is the shift of the XYZ coords to origin, and the second is the shift to Ref 
  * origin. To reproduce the fit perform the first translation (Trans[0...2]), 
  * then rotate (U), then the second translation (Trans[3...5]).
  */
double Frame::RMSD( Frame& Ref, double *U, double *Trans, bool useMassIn) 
{
  double refCOM[3];

  // Rotation will occur around geometric center/center of mass
  // Coords are shifted to common CoM first (Trans0-2), then
  // to original reference location (Trans3-5).
  if (useMassIn) 
    Ref.CenterOfMass(Trans+3, 0, natom_);
  else 
    Ref.GeometricCenter(Trans+3, 0, natom_);
  //fprintf(stderr,"  REF   COM: %lf %lf %lf\n",refCOM[0],refCOM[1],refCOM[2]); //DEBUG

  // Shift to common COM
  refCOM[0] = -Trans[3];
  refCOM[1] = -Trans[4];
  refCOM[2] = -Trans[5];
  Ref.Translate(refCOM);

  double rmsval = RMSD_CenteredRef( Ref, U, Trans, useMassIn );

  //DEBUG
  //printRotTransInfo(U,Trans);
  //fprintf(stdout,"RMS is %lf\n",rms_return);

  return rmsval;
}

// Frame::RMSD_CenteredRef()
/** Get the RMSD of this Frame to given Reference Frame. Ref frame must contain 
  * the same number of atoms as this Frame and should have already been 
  * translated to coordinate origin (neither is checked for in the interest
  * of speed). Put the best-fit rotation matrix in U and the COM translation 
  * vector for Frame in Trans[0-2]. The translation is composed of two XYZ 
  * vectors; the first is the shift of the XYZ coords to origin, and the second 
  * is the shift to Reference origin (should already be set). To reproduce the 
  * fit perform the first translation (Trans[0...2]), rotate (U), then the second
  * translation (Trans[3...5]).
  */
double Frame::RMSD_CenteredRef( Frame const& Ref, double U[9], double Trans[6], bool useMassIn)
{
  double frameCOM[3], rms_return, total_mass;
  double mwss, rot[9], rtr[9];
  double xt,yt,zt,xr,yr,zr;
  double Evector[9], Eigenvalue[3];
  double b[9];
  double cp[3], sig3;

  U[0]=0.0; U[1]=0.0; U[2]=0.0;
  U[3]=0.0; U[4]=0.0; U[5]=0.0;
  U[6]=0.0; U[7]=0.0; U[8]=0.0;
  // Trans is set at the end
 
  // Rotation will occur around geometric center/center of mass
  if (useMassIn)
    total_mass = this->CenterOfMass(frameCOM,0,natom_);
  else
    total_mass = this->GeometricCenter(frameCOM,0,natom_);

  if (total_mass<SMALL) {
    mprinterr("Error: Frame::RMSD: Divide by zero.\n");
    return -1;
  }
  //fprintf(stderr,"  FRAME COM: %lf %lf %lf\n",frameCOM[0],frameCOM[1],frameCOM[2]); //DEBUG

  // Shift to common COM
  frameCOM[0] = -frameCOM[0]; 
  frameCOM[1] = -frameCOM[1]; 
  frameCOM[2] = -frameCOM[2];
  this->Translate(frameCOM);
  //mprintf("  SHIFTED FRAME 0: %lf %lf %lf\n",X[0],X[1],X[2]); //DEBUG
  //mprintf("  SHIFTED REF 0  : %lf %lf %lf\n",Ref.X[0],Ref.X[1],Ref.X[2]); //DEBUG

  // Use Kabsch algorithm to calculate optimum rotation matrix.
  // U = [(RtR)^.5][R^-1]
  mwss=0.0;
  rot[0]=0.0; rot[1]=0.0; rot[2]=0.0;
  rot[3]=0.0; rot[4]=0.0; rot[5]=0.0;
  rot[6]=0.0; rot[7]=0.0; rot[8]=0.0;
  // rtr is set below
  // Calculate covariance matrix of Coords and Reference (R = Xt * Ref)
  Darray::iterator mass = Mass_.begin();
  double atom_mass = 1.0;
  for (int i = 0; i < ncoord_; i += 3)
  {
    xt = X_[i  ];
    yt = X_[i+1];
    zt = X_[i+2];
    xr = Ref.X_[i  ];
    yr = Ref.X_[i+1];
    zr = Ref.X_[i+2];

    // Use atom_mass to hold mass for this atom if specified
    if (useMassIn) 
      atom_mass = *(mass++);

    mwss += atom_mass * ( (xt*xt)+(yt*yt)+(zt*zt)+(xr*xr)+(yr*yr)+(zr*zr) );

    // Calculate the Kabsch matrix: R = (rij) = Sum(yni*xnj)
    rot[0] += atom_mass*xt*xr;
    rot[1] += atom_mass*xt*yr;
    rot[2] += atom_mass*xt*zr;

    rot[3] += atom_mass*yt*xr;
    rot[4] += atom_mass*yt*yr;
    rot[5] += atom_mass*yt*zr;

    rot[6] += atom_mass*zt*xr;
    rot[7] += atom_mass*zt*yr;
    rot[8] += atom_mass*zt*zr;
  }
  mwss *= 0.5;    // E0 = 0.5*Sum(xn^2+yn^2) 

  //DEBUG
  //fprintf(stderr,"ROT:\n%lf %lf %lf\n%lf %lf %lf\n%lf %lf %lf\n",
  //        rot[0],rot[1],rot[2],rot[3],rot[4],rot[5],rot[6],rot[7],rot[8]);
  //fprintf(stderr,"MWSS: %lf\n",mwss);

  // calculate Kabsch matrix multiplied by its transpose: RtR 
  rtr[0] = rot[0]*rot[0] + rot[1]*rot[1] + rot[2]*rot[2];
  rtr[1] = rot[0]*rot[3] + rot[1]*rot[4] + rot[2]*rot[5];
  rtr[2] = rot[0]*rot[6] + rot[1]*rot[7] + rot[2]*rot[8];
  rtr[3] = rot[3]*rot[0] + rot[4]*rot[1] + rot[5]*rot[2];
  rtr[4] = rot[3]*rot[3] + rot[4]*rot[4] + rot[5]*rot[5];
  rtr[5] = rot[3]*rot[6] + rot[4]*rot[7] + rot[5]*rot[8];
  rtr[6] = rot[6]*rot[0] + rot[7]*rot[1] + rot[8]*rot[2];
  rtr[7] = rot[6]*rot[3] + rot[7]*rot[4] + rot[8]*rot[5];
  rtr[8] = rot[6]*rot[6] + rot[7]*rot[7] + rot[8]*rot[8];

  // Diagonalize
  Matrix_3x3 TEMP(rtr);
  if (TEMP.Diagonalize_Sort( Evector, Eigenvalue))
    return 0;
  //if (!Diagonalize_Sort( rtr, Evector, Eigenvalue)) return 0;

  // a3 = a1 x a2
  CROSS_PRODUCT(Evector[6], Evector[7], Evector[8],
                Evector[0], Evector[1], Evector[2],
                Evector[3], Evector[4], Evector[5]);

  // Evector dot transpose rot: b = R . ak
  b[0] = Evector[0]*rot[0] + Evector[1]*rot[3] + Evector[2]*rot[6];
  b[1] = Evector[0]*rot[1] + Evector[1]*rot[4] + Evector[2]*rot[7];
  b[2] = Evector[0]*rot[2] + Evector[1]*rot[5] + Evector[2]*rot[8];
  normalize(b);
  b[3] = Evector[3]*rot[0] + Evector[4]*rot[3] + Evector[5]*rot[6];
  b[4] = Evector[3]*rot[1] + Evector[4]*rot[4] + Evector[5]*rot[7];
  b[5] = Evector[3]*rot[2] + Evector[4]*rot[5] + Evector[5]*rot[8];
  normalize(b+3);
  b[6] = Evector[6]*rot[0] + Evector[7]*rot[3] + Evector[8]*rot[6];
  b[7] = Evector[6]*rot[1] + Evector[7]*rot[4] + Evector[8]*rot[7];
  b[8] = Evector[6]*rot[2] + Evector[7]*rot[5] + Evector[8]*rot[8];
  normalize(b+6);
  /*matrix_multiply_3x3(b, Evector, rot);
  normalize(b);
  normalize(b+3);
  normalize(b+6);*/

 /* b3 = b1 x b2 */
  CROSS_PRODUCT(cp[0], cp[1], cp[2],
                 b[0],  b[1],  b[2],
                 b[3],  b[4],  b[5]);

  if ( (cp[0]*b[6] + cp[1]*b[7] + cp[2]*b[8]) < 0.0 )
    sig3 = -1.0;
  else
    sig3 = 1.0;

  b[6] = cp[0];
  b[7] = cp[1];
  b[8] = cp[2];

  // U has the best rotation 
  U[0] = (Evector[0]*b[0]) + (Evector[3]*b[3]) + (Evector[6]*b[6]);  
  U[1] = (Evector[1]*b[0]) + (Evector[4]*b[3]) + (Evector[7]*b[6]);
  U[2] = (Evector[2]*b[0]) + (Evector[5]*b[3]) + (Evector[8]*b[6]);

  U[3] = (Evector[0]*b[1]) + (Evector[3]*b[4]) + (Evector[6]*b[7]);
  U[4] = (Evector[1]*b[1]) + (Evector[4]*b[4]) + (Evector[7]*b[7]);
  U[5] = (Evector[2]*b[1]) + (Evector[5]*b[4]) + (Evector[8]*b[7]);

  U[6] = (Evector[0]*b[2]) + (Evector[3]*b[5]) + (Evector[6]*b[8]);
  U[7] = (Evector[1]*b[2]) + (Evector[4]*b[5]) + (Evector[7]*b[8]);
  U[8] = (Evector[2]*b[2]) + (Evector[5]*b[5]) + (Evector[8]*b[8]);

  // E=E0-sqrt(mu1)-sqrt(mu2)-sig3*sqrt(mu3) 
  rms_return = mwss;
  rms_return -= sqrt(fabs(Eigenvalue[0]));
  rms_return -= sqrt(fabs(Eigenvalue[1]));
  rms_return -= (sig3*sqrt(fabs(Eigenvalue[2])));

  if (rms_return<0) {
    //fprintf(stderr,"RMS returned is <0 before sqrt, setting to 0 (%lf)\n",rms_return);
    rms_return=0.0;
  } else
    rms_return = sqrt((2.0*rms_return)/total_mass);

  // Translation vectors: Coords are shifted to common CoM first (origin), then
  // to original reference location.
  // frameCOM was negated above to facilitate translation to COM.
  // Reference translation should already be set
  Trans[0] = frameCOM[0];
  Trans[1] = frameCOM[1];
  Trans[2] = frameCOM[2];

  //DEBUG
  //printRotTransInfo(U,Trans);
  //fprintf(stdout,"RMS is %lf\n",rms_return);

  return rms_return;
}

// Frame::RMSD()
/** Calculate RMSD of Frame to Ref with no fitting. Frames must contain
  * same # atoms.
  */
double Frame::RMSD( Frame const& Ref, bool useMass) {
  double rms_return = 0.0;
  double total_mass = 0.0;
  
  Darray::iterator mass = Mass_.begin();
  double atom_mass = 1.0;
  for (int i = 0; i < ncoord_; i += 3)
  {
    double xx = Ref.X_[i  ] - X_[i  ];
    double yy = Ref.X_[i+1] - X_[i+1];
    double zz = Ref.X_[i+2] - X_[i+2];
    if (useMass) 
      atom_mass = *(mass++);
    total_mass += atom_mass;
    rms_return += (atom_mass * (xx*xx + yy*yy + zz*zz));
  }

  if (total_mass<SMALL) {
    mprinterr("Error: no-fit RMSD: Divide by zero.\n");
    return -1;
  }
  if (rms_return < 0) {
    //mprinterr("Error: Frame::RMSD: Negative RMS. Coordinates may be corrupted.\n");
    //return -1;
    //mprinterr("RMS returned is <0 before sqrt, setting to 0 (%lf)\n",rms_return);
    return 0;
  }

  return (sqrt(rms_return / total_mass));
}

// Frame::DISTRMSD()
/** Calcuate the distance RMSD of Frame to Ref. Frames must contain
  * same # of atoms. Should not be called for 0 atoms.
  */
double Frame::DISTRMSD( Frame& Ref ) {
  double TgtDist, RefDist;
  double diff, rms_return;
  double x,y,z;
  int a10,a11,a12;
  int a20,a21,a22; 
  double Ndistances = ((natom_ * natom_) - natom_) / 2;

  rms_return = 0;
  a10 = 0;
  for (int atom1 = 0; atom1 < natom_-1; ++atom1) {
    a20 = a10 + 3;
    for (int atom2 = atom1+1; atom2 < natom_; ++atom2) {
      a11 = a10 + 1;
      a12 = a10 + 2;
      a21 = a20 + 1;
      a22 = a20 + 2;
      // Tgt
      x = X_[a10] - X_[a20];
      x = x * x;
      y = X_[a11] - X_[a21];
      y = y * y;
      z = X_[a12] - X_[a22];
      z = z * z;
      TgtDist = sqrt(x + y + z);
      // Ref
      x = Ref.X_[a10] - Ref.X_[a20];
      x = x * x;
      y = Ref.X_[a11] - Ref.X_[a21];
      y = y * y;
      z = Ref.X_[a12] - Ref.X_[a22];
      z = z * z;
      RefDist = sqrt(x + y + z);
      // DRMSD
      diff = TgtDist - RefDist;
      diff *= diff;
      rms_return += diff;

      a20 += 3;
    } 
    a10 += 3;
  }

  TgtDist = rms_return / Ndistances;
  rms_return = sqrt(TgtDist);

  return rms_return;
}

// Frame::SetAxisOfRotation()
/** Given the central two atoms of a dihedral, calculate
  * a vector (U) which will be the axis for rotating the system around that 
  * dihedral and translate the coordinates (X) to the origin of the new axis.
  */
void Frame::SetAxisOfRotation(double *U, int atom1, int atom2) {
  double A1[3], A2[3];
  int a1 = atom1 * 3;
  int a2 = atom2 * 3;
  
  A1[0] = X_[a1  ];
  A1[1] = X_[a1+1];
  A1[2] = X_[a1+2];
  A2[0] = X_[a2  ];
  A2[1] = X_[a2+1];
  A2[2] = X_[a2+2];

  // Calculate vector of dihedral axis, which will be the new rot. axis
  U[0] = A2[0] - A1[0];
  U[1] = A2[1] - A1[1];
  U[2] = A2[2] - A1[2];

  // Normalize Vector for axis of rotation or scaling will occur!
  normalize(U);

  // Now the rest of the coordinates need to be translated to match the new 
  // rotation axis.
  A1[0] = -A1[0];
  A1[1] = -A1[1];
  A1[2] = -A1[2];
  Translate(A1);
}

// Frame::RotateAroundAxis()
/** Given a vector representing an axis and a magnitude, rotate all 
  * coordinates in the given mask around the axis.
  */
void Frame::RotateAroundAxis(double *T, AtomMask &Rmask) {
  double x,y,z;
  double T0,T1,T2,T3,T4,T5,T6,T7,T8;
  //double T[9];
  // Setup rotation matrix for this axis and given theta
  //calcRotationMatrix(T, U, theta);
  // Rotate
  T0=T[0];
  T1=T[1];
  T2=T[2];
  T3=T[3];
  T4=T[4];
  T5=T[5];
  T6=T[6];
  T7=T[7];
  T8=T[8];
  for (AtomMask::const_iterator atom = Rmask.begin();
                                atom != Rmask.end();
                                atom++)
  {
    int i0 = (*atom) * 3;
    int i1 = i0 + 1;
    int i2 = i1 + 1;
    x=X_[i0]; y=X_[i1]; z=X_[i2];

    X_[i0]=(x*T0) + (y*T1) + (z*T2);
    X_[i1]=(x*T3) + (y*T4) + (z*T5);
    X_[i2]=(x*T6) + (y*T7) + (z*T8);
  }

}

// Frame::CalculateInertia()
void Frame::CalculateInertia(AtomMask& Mask, double* Inertia, double* CXYZ)
{
  double Ivec[6]; // Ivec = xx, yy, zz, xy, xz, yz
  //if (useMass)
    CenterOfMass(CXYZ, Mask);
  //else
  //  GeometricCenter(&Mask, CXYZ);

  // Calculate moments of inertia and products of inertia
  Ivec[0] = 0; // xx
  Ivec[1] = 0; // yy
  Ivec[2] = 0; // zz
  Ivec[3] = 0; // xy
  Ivec[4] = 0; // xz
  Ivec[5] = 0; // yz
  //double *crd = X_;
  for (AtomMask::const_iterator atom = Mask.begin(); atom != Mask.end(); ++atom) {
    int xidx = (*atom) * 3;
    double cx = X_[xidx  ] - CXYZ[0];
    double cy = X_[xidx+1] - CXYZ[1];
    double cz = X_[xidx+2] - CXYZ[2];
    double mass = Mass_[*atom];

    Ivec[0] += mass * ( cy * cy + cz * cz );
    Ivec[1] += mass * ( cx * cx + cz * cz );
    Ivec[2] += mass * ( cx * cx + cy * cy );
    Ivec[3] -= mass * cx * cy;
    Ivec[5] -= mass * cy * cz;
    Ivec[4] -= mass * cx * cz;
  }
  Inertia[0] = Ivec[0];
  Inertia[1] = Ivec[3];
  Inertia[2] = Ivec[4];
  Inertia[3] = Ivec[3];
  Inertia[4] = Ivec[1];
  Inertia[5] = Ivec[5];
  Inertia[6] = Ivec[4];
  Inertia[7] = Ivec[5];
  Inertia[8] = Ivec[2];
}

