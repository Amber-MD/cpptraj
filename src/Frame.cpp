#include <cmath>
#include <cstring>
//#include <algorithm> // std::copy, std::swap
#include "Frame.h"
#include "Constants.h"
#include "vectormath.h"
#include "DistRoutines.h"
#include "TorsionRoutines.h"
#include "CpptrajStdio.h"

const size_t Frame::COORDSIZE_ = 3 * sizeof(double);
const size_t Frame::BOXSIZE_ = 6 * sizeof(double);

/// CONSTRUCTOR
Frame::Frame( ) : 
  natom_(0),
  maxnatom_(0),
  Ncoord_(0),
  X_(NULL),
  V_(NULL),
  Mass_(NULL),
  T_(0.0)
{
  memset(box_, 0, BOXSIZE_);
}

// CONSTRUCTOR
/// Set up for natom
Frame::Frame(int natomIn) :
  natom_(natomIn), 
  maxnatom_(natomIn),
  Ncoord_(natomIn*3),
  X_(NULL),
  V_(NULL),
  Mass_(NULL),
  T_(0)
{
  memset(box_, 0, BOXSIZE_);
  if (Ncoord_ > 0) 
    X_ = new double[ Ncoord_ ];
}

// CONSTRUCTOR
/// Copy given coords
Frame::Frame(double *Xin, int natomIn) :
  natom_(natomIn), 
  maxnatom_(natomIn),
  Ncoord_(natomIn*3),
  X_(NULL),
  V_(NULL),
  Mass_(NULL),
  T_(0) 
{
  memset(box_, 0, BOXSIZE_);
  if (Ncoord_ > 0) {
    X_ = new double[ Ncoord_ ];
    memcpy(X_, Xin, Ncoord_*sizeof(double));
  }
}

// CONSTRUCTOR
/// Set up for natom, copy massIn
Frame::Frame(int natomIn, double *massIn) :
  natom_(natomIn), 
  maxnatom_(natomIn),
  Ncoord_(natomIn*3),
  X_(NULL),
  V_(NULL),
  Mass_(NULL),
  T_(0)
{
  memset(box_, 0, BOXSIZE_);
  if (Ncoord_ > 0) {
    X_ = new double[ Ncoord_ ];
    if (massIn != NULL) {
      Mass_ = new double[ natom_ ];
      memcpy(Mass_, massIn, natom_*sizeof(double));
    }
  }
}

// CONSTRUCTOR
/// Set up for nselected, copy mass according to maskIn
Frame::Frame(AtomMask &maskIn, double *massIn) :
  natom_(maskIn.Nselected()),
  maxnatom_(natom_),
  Ncoord_(natom_*3),
  X_(NULL),
  V_(NULL),
  Mass_(NULL),
  T_(0)
{
  memset(box_, 0, BOXSIZE_);
  if (Ncoord_ > 0) {
    X_ = new double[ Ncoord_ ];
    Mass_ = new double[ natom_ ];
    double *newM = Mass_;
    for (AtomMask::const_iterator atom = maskIn.begin();
                                  atom != maskIn.end(); atom++)
    {
      *newM = massIn[ *atom ];
      ++newM;
    }
  }
}

// CONSTRUCTOR
/// Copy frameIn according to maskIn
Frame::Frame(Frame &frameIn, AtomMask &maskIn) {
  natom_ = maskIn.Nselected();
  Ncoord_ = natom_ * 3;
  if (frameIn.maxnatom_ > natom_)
    maxnatom_ = frameIn.maxnatom_;
  else
    maxnatom_ = natom_;
  // maxnatom may be different than natom. Use maxNcoord for allocation
  int maxNcoord = maxnatom_ * 3;
  memcpy(box_, frameIn.box_, BOXSIZE_);
  T_ = frameIn.T_;
  V_ = NULL;
  Mass_ = NULL;
  X_ = NULL;
  if (frameIn.X_!=NULL) {
    X_ = new double[ maxNcoord ];
    double *newX = X_;
    for (AtomMask::const_iterator atom = maskIn.begin();
                                  atom != maskIn.end(); atom++)
    {
      int oldatom3 = (*atom) * 3;
      memcpy(newX, frameIn.X_ + oldatom3, COORDSIZE_);
      newX += 3;
    }
  }
  if (frameIn.V_ != NULL) {
    V_ = new double[ maxNcoord ];
    double *newV = V_;
    for (AtomMask::const_iterator atom = maskIn.begin();
                                  atom != maskIn.end(); atom++)
    {
      int oldatom3 = (*atom) * 3;
      memcpy(newV, frameIn.V_ + oldatom3, COORDSIZE_);
      newV += 3;
    }
  }
  if (frameIn.Mass_ != NULL) {
    Mass_ = new double[ maxnatom_ ];
    double *newM = Mass_;
    for (AtomMask::const_iterator atom = maskIn.begin();
                                  atom != maskIn.end(); atom++)
    {
      *newM = frameIn.Mass_[ *atom ];
      ++newM;
    }
  }
}

/// DESTRUCTOR
Frame::~Frame( ) {
  if (X_!=NULL) delete[] X_;
  if (V_!=NULL) delete[] V_;
  if (Mass_!=NULL) delete[] Mass_;
}

/// COPY CONSTRUCTOR
// NOTE: The order of variables in the class definition is important
//       when using initializer lists.
/*Frame::Frame(const Frame &rhs) :
  natom_(rhs.natom_),
  maxnatom_(rhs.maxnatom_),
  Ncoord_(rhs.Ncoord_),
  X_(Ncoord_ ? new double[ Ncoord_ ] : NULL),
  V_(rhs.V_!=NULL ? new double[ Ncoord_ ] : NULL),
  Mass_(rhs.Mass_!=NULL ? new double[ natom_ ] : NULL),
  T_(rhs.T_)
{
  std::copy(rhs.X_, rhs.X_ + Ncoord_, X_);
  std::copy(rhs.V_, rhs.V_ + Ncoord_, V_);
  std::copy(rhs.Mass_, rhs.Mass_ + natom_, Mass_);
  std::copy(rhs.box_, rhs.box_ + 6, box_);
}*/
// COPY CONSTRUCTOR
Frame::Frame(const Frame &rhs) :
  natom_(rhs.natom_), 
  maxnatom_(rhs.maxnatom_), 
  Ncoord_(rhs.Ncoord_),
  T_(rhs.T_)
{
  if (natom_ > maxnatom_)
    maxnatom_ = natom_;
  // rhs.maxnatom may not be equal to natom. Use maxNcoord for alloc.
  int maxNcoord = maxnatom_ * 3;
  X_ = NULL;
  V_ = NULL;
  Mass_ = NULL;
  if (rhs.X_!=NULL) {
    X_ = new double[ maxNcoord ];
    memcpy(X_, rhs.X_, Ncoord_*sizeof(double));
  }
  if (rhs.V_!=NULL) {
    V_ = new double[ maxNcoord ];
    memcpy(V_, rhs.V_, Ncoord_*sizeof(double));
  }
  if (rhs.Mass_!=NULL) {
    Mass_ = new double[ maxnatom_ ];
    memcpy(Mass_, rhs.Mass_, natom_*sizeof(double));
  }
}

// Frame::Info()
// For debugging
void Frame::Info(const char *msg) {
  if (msg!=NULL)
    mprintf("\tFrame [%s]:",msg);
  else
    mprintf("\tFrame:");
  mprintf("%i atoms, %i coords, maxnatom=%i",natom_,Ncoord_,maxnatom_);
  if (X_!=NULL) mprintf(" X");
  if (V_!=NULL) mprintf(" V");
  if (Mass_!=NULL) mprintf(" M");
  mprintf("\n");
}

// SWAP
/*void Frame::swap(Frame &first, Frame &second) {
  double firstbox[6];
  // Argument-dependent lookup; bring swap into this scope.
  using std::swap;
  swap(first.natom_, second.natom_);
  swap(first.maxnatom_, second.maxnatom_);
  swap(first.Ncoord_, second.Ncoord_);
  swap(first.X_, second.X_);
  swap(first.V_, second.V_);
  swap(first.Mass_, second.Mass_);
  // Swap box
  std::copy(first.box_, first.box_ + 6, firstbox);
  std::copy(second.box_, second.box_ + 6, first.box_);
  std::copy(firstbox, firstbox + 6, second.box_);
  swap(first.T_, second.T_);
}*/

// Frame::operator=()
/** Assignment operator. Cannot assume that *this has been previously
  * allocated, so smart reallocation really isnt possible. Need to
  * use the SetupFrame routines for that.
  * This is implemented using a copy-and-swap idiom. The copy-constructor
  * will make a copy of rhs, which is then safely switched with *this. 
  * The original data in *this ends up in the copy of rhs, which is freed
  * when the function returns.
  */
Frame &Frame::operator=(const Frame &rhs) {
  //swap(*this, rhs);

  // Check for self assignment
  if ( this == &rhs ) return *this;

  // Copy members that require no allocation
  natom_ = rhs.natom_;
  if (rhs.maxnatom_ > natom_)
    maxnatom_ = rhs.maxnatom_;
  else
    maxnatom_ = natom_;
  // rhs.maxnatom may not be equal to natom. Use maxNcoord for alloc.
  int maxNcoord = maxnatom_ * 3;
  Ncoord_ = rhs.Ncoord_;
  memcpy(box_, rhs.box_, BOXSIZE_);
  T_ = rhs.T_;

  // Deallocate
  if (X_!=NULL) {delete[] X_; X_=NULL;}
  if (V_!=NULL) {delete[] V_; V_=NULL;}
  if (Mass_!=NULL) {delete[] Mass_; Mass_=NULL;}

  // Allocate and copy 
  // Coordinates
  if (rhs.X_ != NULL) {
    X_ = new double[ maxNcoord ];
    memcpy(X_, rhs.X_, Ncoord_*sizeof(double));
  }
  // Velocity
  if (rhs.V_ != NULL) {
    V_ = new double[ maxNcoord ];
    memcpy(V_, rhs.V_, Ncoord_*sizeof(double));
  }
  // Mass
  if (rhs.Mass_ != NULL) {
    Mass_ = new double[ maxnatom_ ];
    memcpy(Mass_, rhs.Mass_, natom_*sizeof(double));
  }

  return *this;
}

// Float array assignment
Frame &Frame::operator=(const std::vector<float> &farray) {
  int f_ncoord = (int) farray.size();
  if (f_ncoord > maxnatom_*3) {
    mprinterr("Error: Frame: Float array size %i > max #coords in frame %i\n",
              f_ncoord, maxnatom_*3);
    return *this;
  }
  natom_ = f_ncoord / 3;
  Ncoord_ = f_ncoord;
  for (int coord = 0; coord < Ncoord_; coord++)
    X_[coord] = (double)farray[coord];
  return *this;
}

// Frame::ConvertToFloat()
std::vector<float> Frame::ConvertToFloat(AtomMask &maskIn) {
  std::vector<float> farray;

  farray.reserve( maskIn.Nselected() * 3 );
  for (AtomMask::const_iterator atom = maskIn.begin();
                                atom != maskIn.end();
                                atom++)
  {
    int atomidx = (*atom) * 3;
    float fval = (float) X_[atomidx];
    farray.push_back( fval );
    fval = (float) X_[atomidx+1];
    farray.push_back( fval );
    fval = (float) X_[atomidx+2];
    farray.push_back( fval );
  }
  return farray;
}

// Frame::ConvertToDouble()
std::vector<double> Frame::ConvertToDouble() {
  return ( std::vector<double> (X_, X_ + Ncoord_ ) );
}

// Frame::DoubleArray()
double *Frame::DoubleArray() {
  if (natom_ < 1) return NULL;
  double *darray = new double[ Ncoord_ ];
  memcpy(darray, X_, Ncoord_*sizeof(double));
  return darray;
} 

// Frame::CoordPtr()
// NOTE: Right now only used for interfacing with older ptraj routines:
//       1) setting up mask in Topology, since the mask parser still uses 
//          double* for coordinates in order to remain compatible with ptraj 
//          routines.
//       2) Setting coordinates for use with ptraj routines (SetReferenceInfo).
// TODO: Get rid of the need for this routine.  
const double *Frame::CoordPtr() {
  return (const double*)X_;
}

// Frame::ConvertToPtrajXYZ()
// NOTE: As with CoordPtr() only used to interface with older ptraj routines.
void Frame::ConvertToPtrajXYZ(double *x_coord, double *y_coord, double *z_coord,
                              double *ptraj_box)
{
  double *Xptr = X_;
  for (int atom = 0; atom < natom_; atom++) {
    x_coord[atom] = *Xptr;
    ++Xptr;
    y_coord[atom] = *Xptr;
    ++Xptr;
    z_coord[atom] = *Xptr;
    ++Xptr;
  }
  // Protect state box coords
  memcpy(ptraj_box, box_, BOXSIZE_);
}

// Frame::SetFromPtrajXYZ()
void Frame::SetFromPtrajXYZ(double *x_coord, double *y_coord, double *z_coord) {
  double *Xptr = X_;
  for (int atom = 0; atom < natom_; atom++) {
    *Xptr = x_coord[atom];
    ++Xptr;
    *Xptr = y_coord[atom];
    ++Xptr;
    *Xptr = z_coord[atom];
    ++Xptr;
  }
}

// Frame::GetAtomXYZ()
// NOTE: Currently only used by molsurf to interface with ATOM data structure
void Frame::GetAtomXYZ(double *Coord, int atom) {
  int i3 = atom * 3;
  double *Xptr = X_ + i3;
  Coord[0] = *Xptr;
  ++Xptr;
  Coord[1] = *Xptr;
  ++Xptr;
  Coord[2] = *Xptr;
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
  memcpy(X_ + Ncoord_, XYZin, COORDSIZE_);
  ++natom_;
  Ncoord_ += 3;
}

// Frame::Natom()
int Frame::Natom() {
  return natom_;
}

// Frame::empty()
bool Frame::empty() {
  return (Ncoord_==0);
}

// Frame::MaxImagedDistance()
// NOTE: Used only by closest. May not be necessary
double Frame::MaxImagedDistance() {
  double maxD = box_[0] + box_[1] + box_[2];
  maxD *= maxD;
  return maxD;
}

// Frame::SetupFrame()
/** Set up frame for given number of atoms, no mass or velocity information.
  * Only reallocate memory if natomIn > maxnatom.
  */
int Frame::SetupFrame(int natomIn) {
  return ReallocateFrame(natomIn, NULL, false);
}

// Frame::SetupFrame()
/** Set up frame for given number of atoms. Store mass information
  * if passed in. No velocity info. Only reallocate memory if 
  * natomIn > maxnatom.
  */
int Frame::SetupFrame(int natomIn, double *massIn) {
  bool hasMass = (massIn!=NULL);
  ReallocateFrame(natomIn, hasMass, false);
  if (hasMass)
    memcpy(Mass_, massIn, natom_*sizeof(double));
  return 0;
}

// Frame::SetupFrameV()
/** Set up frame for given number of atoms. Store mass information
  * if passed in. Set up for velocity info if hasVelocity = true.
  * Only reallocate memory if natomIn > maxnatom
  */
int Frame::SetupFrameV(int natomIn, double *massIn, bool hasVelocity) {
  bool hasMass = (massIn!=NULL);
  ReallocateFrame(natomIn, hasMass, hasVelocity);
  if (hasMass)
    memcpy(Mass_, massIn, natom_*sizeof(double));
  return 0;
}

// Frame::SetupFrameFromMask()
/** Set up frame to hold # selected atoms in given mask. If mass 
  * information is passed in store the masses corresponding to
  * selected atoms in the mask. No velocity info. Only reallocate
  * memory if Nselected > maxnatom.
  */
int Frame::SetupFrameFromMask(AtomMask &maskIn, double *massIn) {
  bool hasMass = (massIn!=NULL);
  ReallocateFrame(maskIn.Nselected(), hasMass, false);
  if (hasMass) {
    double *massptr = Mass_;
    for (AtomMask::const_iterator atom = maskIn.begin();
                                  atom != maskIn.end(); atom++)
    {
       *massptr = massIn[ *atom ];
       ++massptr;
    }
  }
  return 0;
}

// Frame::ReallocateFrame
/** Set up frame for given number of atoms. Set up space for storing mass 
  * information if hasMass = true. Set up space for storing velocity if 
  * hasVelocity = true. Only reallocate memory if natomIn > maxnatom.
  */ 
int Frame::ReallocateFrame(int natomIn, bool hasMass, bool hasVelocity) {
  bool Reallocate = false;
  // Set # atoms and coords
  natom_ = natomIn;
  Ncoord_ = natom_ * 3;
  // Since box/T info will not necessarily be read in, initialize
  box_[0]=0; box_[1]=0; box_[2]=0; box_[3]=0; box_[4]=0; box_[5]=0;
  T_ = 0.0;
  // Check if reallocation must occur, reallocate coords if necessary
  // ASSUMES COORDS ARE ALWAYS ALLOCD; should be safe
  if (natomIn > maxnatom_) {
    Reallocate = true;
    maxnatom_ = natomIn;
    if (X_!=NULL) delete[] X_;
    X_ = new double[ Ncoord_ ];
  }
  // Assume coords are always read in, dont initialize.
  bool m_allocd = (Mass_ != NULL);
  bool v_allocd = (V_ != NULL);
  // Mass
  if (!hasMass) {
    if (m_allocd) delete[] Mass_;
    Mass_ = NULL;
  } else {
    if (Reallocate || !m_allocd) {
      if (m_allocd) delete[] Mass_;
      Mass_ = new double[ natom_ ];
    }
    // if hasMass mass should be passed in right after this routine is called.
  }
  // Velocity
  if (!hasVelocity) {
    if (v_allocd) delete[] V_;
    V_ = NULL;
  } else {
    if (Reallocate || ! v_allocd) {
      if (v_allocd) delete[] V_;
      V_ = new double[ Ncoord_ ];
    }
    // Since V will not necessarily be read in, initialize it
    memset(V_, 0, Ncoord_ * sizeof(double));
  }
  return 0;
}

// Frame::SetCoordinates()
/** Copy only coordinates from input frame to this frame based
  * on selected atoms in mask.
  */
void Frame::SetCoordinates(Frame &frameIn, AtomMask &maskIn) {
  if (maskIn.Nselected() > maxnatom_) {
    mprinterr("Error: Frame::SetCoordinates: Mask [%s] selected (%i) > max natom\n",
              maskIn.MaskString(),maskIn.Nselected());
    mprinterr("       for this frame (%i)\n",maxnatom_);
    return;
  }
  natom_ = maskIn.Nselected();
  Ncoord_ = natom_ * 3;
  double *newX = X_;
  for (AtomMask::const_iterator atom = maskIn.begin();
                                atom != maskIn.end(); atom++)
  {
    int oldatom3 = (*atom) * 3;
    memcpy(newX, frameIn.X_ + oldatom3, COORDSIZE_);
    newX += 3;
  }
}

//------------------------------------------------------------------------------
// NOTE: Should these map routines be part of a child class used only by AtomMap?
// Frame::SetCoordinates()
/** Reorder atoms in the given target frame so that it matches the given
  * atom map with format (tgtAtom = Map[refAtom]). End result is that
  * this frame will be in the same order as the original reference. Assumes
  * that the map is complete (i.e. no unmapped atoms).
  */
// TODO: Error checking. Map should be something besides an int*
// NOTE: Should this reorder masses as well?
void Frame::SetCoordinates(Frame &tgtIn, std::vector<int>& mapIn) {
  if (tgtIn.natom_ > maxnatom_) {
    mprinterr("Error: Frame::SetCoordinates: Input map frame atoms (%i) > max natom\n",
              tgtIn.natom_);
    mprinterr("       for this frame (%i)\n",maxnatom_);
    return;
  }
  double *newX = X_;
  natom_ = tgtIn.natom_;
  Ncoord_ = tgtIn.Ncoord_;
  for (int refatom = 0; refatom < tgtIn.natom_; refatom++) {
    int tgtatom3 = mapIn[refatom] * 3;
    memcpy(newX, tgtIn.X_ + tgtatom3, COORDSIZE_);
    newX += 3;
  }
}

// Frame::SetReferenceByMap()
/** Set this frame to include only atoms from the given reference that are 
  * mapped.
  */
// TODO: Error Checking
void Frame::SetReferenceByMap(Frame &refIn, std::vector<int>& mapIn) {
  double *newX = X_;
  double *refX = refIn.X_;
  natom_ = 0;
  for (int refatom = 0; refatom < (int)mapIn.size(); refatom++) {
    if (mapIn[refatom] != -1) {
      memcpy(newX, refX, COORDSIZE_);
      newX += 3;
      ++natom_;
    }
    refX += 3;
  }
  Ncoord_ = natom_ * 3;
}

// Frame::SetTargetByMap()
/** Set this frame to include only atoms from the given target frame, remapped
  * according to the given atom map.
  */
void Frame::SetTargetByMap(Frame &tgtIn, std::vector<int>& mapIn) {
  double *newX = X_;
  natom_ = 0;
  for (int refatom = 0; refatom < (int)mapIn.size(); refatom++) {
    int tgtatom = mapIn[refatom];
    if (tgtatom!=-1) {
      int tgtatom3 = tgtatom * 3;
      memcpy(newX, tgtIn.X_ + tgtatom3, COORDSIZE_);
      newX += 3;
      ++natom_;
    }
  }
  Ncoord_ = natom_ * 3;
}

// -----------------------------------------------------------------------------
// Frame::SetCoordinates()
/** Copy only coordinates from input frame to this frame.
  */
void Frame::SetCoordinates(Frame &frameIn) {
  if (frameIn.natom_ > maxnatom_) {
    mprinterr("Error: Frame::SetCoordinates: Input frame atoms (%i) > max natom\n",
              frameIn.natom_);
    mprinterr("       for this frame (%i)\n",maxnatom_);
    return;
  }
  natom_ = frameIn.natom_;
  Ncoord_ = frameIn.Ncoord_;
  memcpy(X_, frameIn.X_, Ncoord_ * sizeof(double));
}

// Frame::SetFrame()
/** Copy entire input frame (including velocity if defined in both) to this 
  * frame based on selected atoms in mask. Assumes coords and mass exist.
  * Everything must be already be allocated. Frame can be resized up to 
  * maxnatom.
  */
void Frame::SetFrame(Frame &frameIn, AtomMask &maskIn) {
  if (maskIn.Nselected() > maxnatom_) {
    mprinterr("Error: Frame::SetFrame: Mask [%s] selected (%i) > max natom\n",
              maskIn.MaskString(),maskIn.Nselected());
    mprinterr("       for this frame (%i)\n",maxnatom_);
    return;
  }
  natom_ = maskIn.Nselected();
  Ncoord_ = natom_ * 3;
  // Copy coords and mass
  double *newX = X_;
  double *newM = Mass_;
  for (AtomMask::const_iterator atom = maskIn.begin();
                                atom != maskIn.end(); atom++)
  {
    int oldatom3 = (*atom) * 3;
    memcpy(newX, frameIn.X_ + oldatom3, COORDSIZE_);
    *newM = frameIn.Mass_[ *atom ];
    newX += 3;
    ++newM;
  }
  // Copy box/T
  memcpy(box_, frameIn.box_, BOXSIZE_);
  T_ = frameIn.T_;
  // Copy velocity
  double *newV = V_;
  if (frameIn.V_ != NULL && V_ != NULL) {
    for (AtomMask::const_iterator atom = maskIn.begin();
                                atom != maskIn.end(); atom++)
    {
      int oldatom3 = (*atom) * 3;
      memcpy(newV, frameIn.V_ + oldatom3, COORDSIZE_);
      newV += 3;
    }
  }
}

// Frame::FrameCopy()
/** Return a copy of the frame */
Frame *Frame::FrameCopy( ) {
  Frame *newFrame;

  newFrame=new Frame( );
  newFrame->natom_ = this->natom_;
  newFrame->maxnatom_ = this->maxnatom_;
  int maxNcoord = maxnatom_ * 3;
  newFrame->Ncoord_ = this->Ncoord_;
  newFrame->X_ = new double[ maxNcoord ];
  memcpy(newFrame->X_, this->X_, Ncoord_ * sizeof(double));
  memcpy(newFrame->box_, this->box_, BOXSIZE_);
  newFrame->T_ = this->T_;
  if (this->V_!=NULL) {
    newFrame->V_ = new double[ maxNcoord ];
    memcpy(newFrame->V_, this->V_, Ncoord_ * sizeof(double));
  }
  if (this->Mass_!=NULL) {
    newFrame->Mass_ = new double[ maxnatom_ ];
    memcpy(newFrame->Mass_, this->Mass_, natom_ * sizeof(double));
  }

  return newFrame;
}

/* ------------------- Coordinate Manipulation Routines --------------------- */
// Frame::ZeroCoords()
/** Set all coords to 0.0 */
void Frame::ZeroCoords( ) {
  memset(X_, 0, Ncoord_ * sizeof(double));
}

// Frame::operator+=()
Frame &Frame::operator+=(const Frame &rhs) {
  // For now ensure same natom
  if (natom_ != rhs.natom_) {
    mprinterr("Error: Frame::operator+=: Attempting to add 2 frames with different natom.\n");
    return *this;
  }
  for (int i = 0; i < Ncoord_; i++)
    X_[i] += rhs.X_[i];
  return *this;
}

// Frame::operator-=()
Frame &Frame::operator-=(const Frame &rhs) {
  // For now ensure same natom
  if (natom_ != rhs.natom_) {
    mprinterr("Error: Frame::operator-=: Attempting to add 2 frames with different natom.\n");
    return *this;
  }
  for (int i = 0; i < Ncoord_; i++)
    X_[i] -= rhs.X_[i];
  return *this;
}

// Frame::operator*=()
Frame &Frame::operator*=(const Frame &rhs) {
  // For now ensure same natom
  if (natom_ != rhs.natom_) {
    mprinterr("Error: Frame::operator*=: Attempting to add 2 frames with different natom.\n");
    return *this;
  }
  for (int i = 0; i < Ncoord_; i++)
    X_[i] *= rhs.X_[i];
  return *this;
}

// Frame::operator*()
// NOTE: return const instance and not const reference to disallow
//       nonsense statements like (a * b) = c
const Frame Frame::operator*(const Frame &rhs) const {
  return (Frame(*this) *= rhs);
}

// Frame::Divide()
/** Divide all coord values of this frame by divisor and store in dividend.
  */
int Frame::Divide(Frame &dividend, double divisor) {
  if (divisor < SMALL) {
    mprinterr("Error: Frame::Divide(Frame,divisor): Detected divide by 0.\n");
    return 1;
  }
  for (int i=0; i<Ncoord_; i++)
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
  for (int coord=0; coord < Ncoord_; coord++)
    X_[coord] /= divisor;
}

// Frame::AddByMask()
/** Increment atoms in this frame by the selected atoms in given frame.
  */
// TODO: Add error checking
void Frame::AddByMask(Frame &frameIn, AtomMask &maskIn) {
  // SHOULD CHECK BOUNDS!
  double *newX = X_;
  for (AtomMask::const_iterator atom = maskIn.begin();
                                atom != maskIn.end(); atom++)
  {
    double *Xin = frameIn.X_ + ((*atom) * 3);
    *newX += *Xin; // Xcoord
    ++newX;
    ++Xin;
    *newX += *Xin; // Ycoord
    ++newX;
    ++Xin;
    *newX += *Xin; // Zcoord
    ++newX;
  }
}

// Frame::Translate()
/** Translate all coords by Vec.  */
void Frame::Translate(double * Vec) {
  double Vec0, Vec1, Vec2;

  Vec0=Vec[0];
  Vec1=Vec[1];
  Vec2=Vec[2];
  for (int i=0; i<Ncoord_; i+=3) {
    X_[i  ]+=Vec0;
    X_[i+1]+=Vec1;
    X_[i+2]+=Vec2;
  }
}

// Frame::Translate()
/** Translate atoms in range by Vec. */
// NOTE: SHOULD CHECK BOUNDS! 
void Frame::Translate(double *Vec, int firstAtom, int lastAtom) {
  int startatom3 = firstAtom * 3;
  int lastatom3 = lastAtom * 3;
  double V0 = Vec[0];
  double V1 = Vec[1];
  double V2 = Vec[2];

  for (int atom3 = startatom3; atom3 < lastatom3; atom3 += 3) {
    X_[atom3  ] += V0;
    X_[atom3+1] += V1;
    X_[atom3+2] += V2;
  }
}

// Frame::Trans_Rot_Trans()
/** Given an array Vec of size 6 containing two translations:
  *   T0x T0y T0z T1x T1y T1z
  * and a rotation matrix T, apply the first translation, then
  * the rotation, then the second translation.
  */
void Frame::Trans_Rot_Trans(double *Vec, double *T) {
  double Vec0, Vec1, Vec2, Vec3, Vec4, Vec5;
  double x, y, z;
  double T0,T1,T2,T3,T4,T5,T6,T7,T8;

  Vec0=Vec[0];
  Vec1=Vec[1];
  Vec2=Vec[2];
  Vec3=Vec[3];
  Vec4=Vec[4];
  Vec5=Vec[5];

  T0=T[0]; 
  T1=T[1]; 
  T2=T[2]; 
  T3=T[3]; 
  T4=T[4]; 
  T5=T[5]; 
  T6=T[6]; 
  T7=T[7]; 
  T8=T[8]; 

  for (int i=0; i<Ncoord_; i+=3) {
    x = X_[i  ] + Vec0;
    y = X_[i+1] + Vec1;
    z = X_[i+2] + Vec2;

    X_[i  ]=(x*T0) + (y*T1) + (z*T2) + Vec3;
    X_[i+1]=(x*T3) + (y*T4) + (z*T5) + Vec4;
    X_[i+2]=(x*T6) + (y*T7) + (z*T8) + Vec5;
  }
}

// Frame::Rotate()
/** Multiply natomx3 matrix X by 3x3 matrix T. If T is a rotation matrix
  * this rotates the coords in X. 
  */
void Frame::Rotate(double *T) {
  double x,y,z;
  double T0,T1,T2,T3,T4,T5,T6,T7,T8;
 
  T0=T[0]; 
  T1=T[1]; 
  T2=T[2]; 
  T3=T[3]; 
  T4=T[4]; 
  T5=T[5]; 
  T6=T[6]; 
  T7=T[7]; 
  T8=T[8]; 
  for (int i=0; i<Ncoord_; i+=3) {
    x=X_[i]; y=X_[i+1]; z=X_[i+2];

    X_[i  ]=(x*T0) + (y*T1) + (z*T2);
    X_[i+1]=(x*T3) + (y*T4) + (z*T5);
    X_[i+2]=(x*T6) + (y*T7) + (z*T8);
  }
} 

// Frame::InverseRotate()
/** Multiply natomx3 matrix X by transpose of 3x3 matrix T. If T is a rotation
  * matrix this rotates the coords in X in the opposite direction.
  */
void Frame::InverseRotate(double *T) {
  double x,y,z;

  for (int i=0; i<Ncoord_; i+=3) {
    x=X_[i]; y=X_[i+1]; z=X_[i+2];

    X_[i  ]=(x*T[0]) + (y*T[3]) + (z*T[6]);
    X_[i+1]=(x*T[1]) + (y*T[4]) + (z*T[7]);
    X_[i+2]=(x*T[2]) + (y*T[5]) + (z*T[8]);
  }
}

// Frame::Center()
/** Center coordinates to center of coordinates in Mask w.r.t. given XYZ in
  * boxcoord. When called from Action_Center boxcoord will be either origin 
  * or box center. Use geometric center if mass is NULL, otherwise center 
  * of mass will be used.
  */
void Frame::Center(AtomMask *Mask, bool origin, bool useMassIn) {
  double center[3];

  if (useMassIn)
    this->CenterOfMass(Mask, center);
  else
    this->GeometricCenter(Mask, center);
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
void Frame::CenterReference(double *Trans, bool useMassIn) {
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

  this->GeometricCenter(frameCOM,0,natom_);
  //mprinterr("  FRAME COM: %lf %lf %lf\n",frameCOM[0],frameCOM[1],frameCOM[2]); //DEBUG
  
  // Shift to common COM
  frameCOM[0]=-frameCOM[0]; frameCOM[1]=-frameCOM[1]; frameCOM[2]=-frameCOM[2];
  this->Translate(frameCOM);
}

// Frame::ImageNonortho()
void Frame::ImageNonortho(bool origin, AtomMask *ComMask, bool truncoct, bool center,
                          bool useMass, std::vector<int> &AtomPairs)
{
  double ucell[9], recip[9], boxTrans[3], Coord[3];
  double fc[3], ffc[3];
  // fcom and ixyz only needed for truncoct
  double fcom[3];
  int ixyz[3];

  BoxToRecip(ucell, recip);
  // Set up centering if putting nonortho cell into familiar trunc. oct. shape
  if (truncoct) {
    if (ComMask!=NULL) {
      // Use center of atoms in mask
      if (useMass)
        CenterOfMass(ComMask, fcom);
      else
        GeometricCenter(ComMask, fcom);
    } else if (origin) {
      // Use origin
      fcom[0] = 0;
      fcom[1] = 0;
      fcom[2] = 0;
    } else {
      // Use box center
      fcom[0] = box_[0] / 2;
      fcom[1] = box_[1] / 2;
      fcom[2] = box_[2] / 2;
    }
    //fprintf(stdout,"DEBUG: fcom = %lf %lf %lf\n",fcom[0],fcom[1],fcom[2]);
  }

  // Loop over atom pairs
  for (std::vector<int>::iterator atom = AtomPairs.begin();
                                  atom != AtomPairs.end();
                                  atom++)
  {
    int firstAtom = *atom;
    ++atom;
    int lastAtom = *atom;
    //if (debug>2)
    //  mprintf( "  IMAGE processing atoms %i to %i\n", firstAtom+1, lastAtom);
    // boxTrans will hold calculated translation needed to move atoms back into box
    boxTrans[0] = 0;
    boxTrans[1] = 0;
    boxTrans[2] = 0;
    // Set up Coord with position to check for imaging based on first atom or 
    // center of mass of atoms first to last.
    if (center) {
      if (useMass)
        CenterOfMass(Coord,firstAtom,lastAtom);
      else
        GeometricCenter(Coord,firstAtom,lastAtom);
    } else {
      int atomidx = firstAtom * 3;
      Coord[0] = X_[atomidx];
      ++atomidx;
      Coord[1] = X_[atomidx];
      ++atomidx;
      Coord[2] = X_[atomidx];
    }
   
    fc[0]=(Coord[0]*recip[0]) + (Coord[1]*recip[1]) + (Coord[2]*recip[2]);
    fc[1]=(Coord[0]*recip[3]) + (Coord[1]*recip[4]) + (Coord[2]*recip[5]);
    fc[2]=(Coord[0]*recip[6]) + (Coord[1]*recip[7]) + (Coord[2]*recip[8]);

    if ( origin ) {
      fc[0] += 0.5;
      fc[1] += 0.5;
      fc[2] += 0.5;
    }

    ffc[0] = floor(fc[0]);
    ffc[1] = floor(fc[1]);
    ffc[2] = floor(fc[2]);

    boxTrans[0] -= (ffc[0]*ucell[0] + ffc[1]*ucell[3] + ffc[2]*ucell[6]);
    boxTrans[1] -= (ffc[0]*ucell[1] + ffc[1]*ucell[4] + ffc[2]*ucell[7]);
    boxTrans[2] -= (ffc[0]*ucell[2] + ffc[1]*ucell[5] + ffc[2]*ucell[8]);

    // Put into familiar trunc. oct. shape
    if (truncoct) {
      Coord[0] += boxTrans[0];
      Coord[1] += boxTrans[1];
      Coord[2] += boxTrans[2];
      MinImageNonOrtho2(Coord, fcom, box_, (int)origin, ixyz, ucell, recip);
      if (ixyz[0] != 0 || ixyz[1] != 0 || ixyz[2] != 0) {
        boxTrans[0] += (ixyz[0]*ucell[0] + ixyz[1]*ucell[3] + ixyz[2]*ucell[6]);
        boxTrans[1] += (ixyz[0]*ucell[1] + ixyz[1]*ucell[4] + ixyz[2]*ucell[7]);
        boxTrans[2] += (ixyz[0]*ucell[2] + ixyz[1]*ucell[5] + ixyz[2]*ucell[8]);

        //if (debug > 2)
        //  mprintf( "  IMAGING, FAMILIAR OFFSETS ARE %i %i %i\n", 
        //          ixyz[0], ixyz[1], ixyz[2]);
      }
    }

    Translate(boxTrans, firstAtom, lastAtom);

  } // END loop over atom pairs
}

// Frame::ImageOrtho()
void Frame::ImageOrtho(bool origin, bool center, bool useMass, std::vector<int> &AtomPairs) 
{
  double bp[3], bm[3], boxTrans[3], Coord[3];

  // Set up boundary information for orthorhombic cell
  if (origin) {
    bp[0] = box_[0] / 2;
    bp[1] = box_[1] / 2;
    bp[2] = box_[2] / 2;
    bm[0] = -bp[0];
    bm[1] = -bp[1];
    bm[2] = -bp[2];
  } else {
    bp[0] = box_[0];
    bp[1] = box_[1];
    bp[2] = box_[2];
    bm[0] = 0;
    bm[1] = 0;
    bm[2] = 0;
  }

  // Loop over atom pairs
  for (std::vector<int>::iterator atom = AtomPairs.begin();
                                  atom != AtomPairs.end();
                                  atom++)
  {
    int firstAtom = *atom;
    ++atom;
    int lastAtom = *atom;
    //if (debug>2)
    //  mprintf( "  IMAGE processing atoms %i to %i\n", firstAtom+1, lastAtom);
    // boxTrans will hold calculated translation needed to move atoms back into box
    boxTrans[0] = 0;
    boxTrans[1] = 0;
    boxTrans[2] = 0;
    // Set up Coord with position to check for imaging based on first atom or 
    // center of mass of atoms first to last.
    if (center) {
      if (useMass)
        CenterOfMass(Coord,firstAtom,lastAtom);
      else
        GeometricCenter(Coord,firstAtom,lastAtom);
    } else {
      int atomidx = firstAtom * 3;
      Coord[0] = X_[atomidx];
      ++atomidx;
      Coord[1] = X_[atomidx];
      ++atomidx;
      Coord[2] = X_[atomidx];
    }

    // Determine how far Coord is out of box
    for (int i=0; i < 3; i++) {
      while (Coord[i] < bm[i]) {
        Coord[i] += box_[i];
        boxTrans[i] += box_[i];
      }
      while (Coord[i] > bp[i]) {
        Coord[i] -= box_[i];
        boxTrans[i] -= box_[i];
      }
    }

    // Translate atoms according to Coord
    Translate(boxTrans,firstAtom,lastAtom);
  } // END loop over atom pairs
}

// Frame::printAtomCoord()
/** Print XYZ coords of given atom */
void Frame::printAtomCoord(int atom) {
  int atmidx = atom * 3;
  if (atmidx >= Ncoord_) return;
  mprintf("ATOM %i: %lf %lf %lf\n",atom,
          X_[atmidx],X_[atmidx+1],X_[atmidx+2]);
}

// ----------------- Center of Mass Calculation Routines -------------------- 
// Frame::CenterOfMass()
/** Given an AtomMask put center of mass of atoms in mask into Coord. Return 
  * sum of masses in Mask.
  */
double Frame::CenterOfMass(AtomMask *Mask, double *Coord) {
  double sumMass,mass,Coord0,Coord1,Coord2;
 
  Coord0=0.0;
  Coord1=0.0;
  Coord2=0.0;
  sumMass=0.0;

  for (AtomMask::const_iterator atom = Mask->begin();
                                atom != Mask->end();
                                atom++)
  {
      int atmidx = (*atom) * 3;
      mass = Mass_[*atom];
      sumMass += mass;
      Coord0 += (X_[atmidx  ] * mass);
      Coord1 += (X_[atmidx+1] * mass);
      Coord2 += (X_[atmidx+2] * mass);
  }

  // NOTE: Not using == since it is unreliable for floating point numbers.
  // Should NEVER have a mass smaller than SMALL (vectormath.h)
  if (sumMass < SMALL) return 0;

  Coord[0] = Coord0 / sumMass;
  Coord[1] = Coord1 / sumMass;
  Coord[2] = Coord2 / sumMass;
  return sumMass;
}

// Frame::GeometricCenter()
/** Given an AtomMask put geometric center of atoms in mask into Coord. Return 
  * #atoms in Mask.
  */
double Frame::GeometricCenter(AtomMask *Mask, double *Coord) {
  double sumMass,Coord0,Coord1,Coord2;
 
  Coord0=0.0;
  Coord1=0.0;
  Coord2=0.0;

  for (AtomMask::const_iterator atom = Mask->begin(); 
                                atom != Mask->end();
                                atom++)
  {
      int atmidx = (*atom) * 3;
      Coord0 += (X_[atmidx  ]);
      Coord1 += (X_[atmidx+1]);
      Coord2 += (X_[atmidx+2]);
  }

  sumMass=(double) Mask->Nselected();

  // NOTE: Not using == since it is unreliable for floating point numbers.
  // Should NEVER have a mass smaller than SMALL (vectormath.h)
  if (sumMass < SMALL) return 0;

  Coord[0] = Coord0 / sumMass;
  Coord[1] = Coord1 / sumMass;
  Coord[2] = Coord2 / sumMass;
  return sumMass;
}

// Frame::CenterOfMass()
/** Put center of mass of all atoms between start and stop in frame into 
  * Coord. Return sum of masses.
  */
double Frame::CenterOfMass(double *Coord, int startAtom, int stopAtom) {  
  int i,m;
  int startAtom3, stopAtom3;
  double sumMass,mass,Coord0,Coord1,Coord2;
 
  Coord[0]=0.0;
  Coord[1]=0.0;
  Coord[2]=0.0; 
  Coord0=0.0;
  Coord1=0.0;
  Coord2=0.0;
  sumMass=0.0;
  m=startAtom;
  startAtom3 = startAtom * 3;
  stopAtom3 = stopAtom * 3;
  
  for (i=startAtom3; i<stopAtom3; i+=3) {
    mass=Mass_[m++];
    sumMass+=mass;
    Coord0+=(X_[i  ] * mass);
    Coord1+=(X_[i+1] * mass);
    Coord2+=(X_[i+2] * mass);
  }

  // NOTE: Not using == since it is unreliable for floating point numbers.
  // Should NEVER have a mass smaller than SMALL (vectormath.h)
  if (sumMass < SMALL) return 0;

  Coord[0] = Coord0 / sumMass;
  Coord[1] = Coord1 / sumMass;
  Coord[2] = Coord2 / sumMass;

  return sumMass;
}

// Frame::GeometricCenter()
/** Put geometric center of all atoms between start and stop in frame into 
  * Coord. Return #atoms used in calc.
  */
double Frame::GeometricCenter(double *Coord, int startAtom, int stopAtom) {  
  int i;
  int startAtom3, stopAtom3;
  double sumMass,Coord0,Coord1,Coord2;
 
  Coord[0]=0.0;
  Coord[1]=0.0;
  Coord[2]=0.0; 
  Coord0=0.0;
  Coord1=0.0;
  Coord2=0.0;
  startAtom3 = startAtom * 3;
  stopAtom3 = stopAtom * 3;
  
  for (i=startAtom3; i<stopAtom3; i+=3) {
    Coord0+=(X_[i  ]);
    Coord1+=(X_[i+1]);
    Coord2+=(X_[i+2]);
  }

  i = stopAtom - startAtom;
  sumMass=(double) i;

  // NOTE: Not using == since it is unreliable for floating point numbers.
  // Should NEVER have a mass smaller than SMALL (vectormath.h)
  if (sumMass < SMALL) return 0;

  Coord[0] = Coord0 / sumMass;
  Coord[1] = Coord1 / sumMass;
  Coord[2] = Coord2 / sumMass;

  return sumMass;
}

// -------------------- Coordinate Calculation Routines --------------------- 
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

// Frame::DIST2()
/** Call the appropriate distance calc for atoms in Mask1 and Mask2 based on
  * given box type.
  *   0 = None
  *   1 = Orthorhombic
  *   2 = Non-orthorhombic
  * Based on useMassIn, calculate geometric center (false) or center of mass 
  * (true) of the atoms in each mask.
  */
double Frame::DIST2(AtomMask *Mask1, AtomMask *Mask2, bool useMassIn, int boxType,
                    double *ucell, double *recip) {
  double a1[3], a2[3];

  if (useMassIn) {
    CenterOfMass(Mask1, a1);
    CenterOfMass(Mask2, a2);
  } else {
    GeometricCenter(Mask1, a1);
    GeometricCenter(Mask2, a2);
  }

  if (boxType == 0) 
    return DIST2_NoImage(a1, a2);
  else if (boxType == 1) 
    return DIST2_ImageOrtho(a1, a2, this->box_);
  else if (boxType == 2) 
    return DIST2_ImageNonOrtho(a1, a2, ucell, recip);

  mprintf("    Error: Frame::DIST: Unrecognized box type (%i)\n.", boxType);

  return (-1.0);
}

// Frame::DIST2()
/** Return the distance between atoms A1 and A2 with optional imaging.
  *   0 = None
  *   1 = Orthorhombic
  *   2 = Non-orthorhombic
  */
double Frame::DIST2(int A1, int A2, int boxType, double *ucell, double *recip) {
  int atom3;
  double a1[3], a2[3];

  atom3 = A1 * 3;
  a1[0] = X_[atom3  ];
  a1[1] = X_[atom3+1];
  a1[2] = X_[atom3+2];
  atom3 = A2 * 3;
  a2[0] = X_[atom3  ];
  a2[1] = X_[atom3+1];
  a2[2] = X_[atom3+2];

  if (boxType == 0)
    return DIST2_NoImage(a1, a2);
  else if (boxType == 1)
    return DIST2_ImageOrtho(a1, a2, this->box_);
  else if (boxType == 2) 
    return DIST2_ImageNonOrtho(a1, a2, ucell, recip);

  mprintf("    Error: Frame::DIST: Unrecognized box type (%i)\n.", boxType);

  return (-1.0);
}

// Frame::DIST2()
/** Return the distance between a point and atom A2 with optional imaging.
  *   0 = None
  *   1 = Orthorhombic
  *   2 = Non-orthorhombic
  */
double Frame::DIST2(double *a1, int A2, int boxType, double *ucell, double *recip) {
  int atom3;
  double a2[3];
  
  atom3 = A2 * 3;
  a2[0] = X_[atom3  ];
  a2[1] = X_[atom3+1];
  a2[2] = X_[atom3+2];
  
  if (boxType == 0)
    return DIST2_NoImage(a1, a2);
  else if (boxType == 1)
    return DIST2_ImageOrtho(a1, a2, this->box_);
  else if (boxType == 2) 
    return DIST2_ImageNonOrtho(a1, a2, ucell, recip);

  mprintf("    Error: Frame::DIST: Unrecognized box type (%i)\n.", boxType);

  return (-1.0);
}

// Frame::DIST()
/** Return the distance between atoms A1 and A2, no imaging.
  */
double Frame::DIST(int A1, int A2) {
  int i, j; // Actual indices into X
  double x,y,z,D;

  i = A1 * 3;
  j = A2 * 3;

  x = X_[i  ] - X_[j  ];
  y = X_[i+1] - X_[j+1];
  z = X_[i+2] - X_[j+2];

  x=x*x;
  y=y*y;
  z=z*z;

  D=sqrt(x + y + z);
  return D;
}

// Frame::DIST2()
/** Return the distance squared between atoms A1 and A2, no imaging.
  */
double Frame::DIST2(int A1, int A2) {
  int i, j; // Actual indices into X
  double x,y,z,D2;

  i = A1 * 3;
  j = A2 * 3;

  x = X_[i  ] - X_[j  ];
  y = X_[i+1] - X_[j+1];
  z = X_[i+2] - X_[j+2];

  x=x*x;
  y=y*y;
  z=z*z;

  D2=(x + y + z);
  return D2;
}

// Frame::COORDDIST()
/** Return the distance between atoms i and j, no imaging. i and j
  * should be actual indices into the coord array (i.e. atom# * 3).
  */
double Frame::COORDDIST(int i, int j) {
  double x,y,z,D;

  x = X_[i  ] - X_[j  ];
  y = X_[i+1] - X_[j+1];
  z = X_[i+2] - X_[j+2];

  x=x*x;
  y=y*y;
  z=z*z;

  D=sqrt(x + y + z);
  return D;
}

// Frame::COORDDIST2()
/** Return the distance between atoms i and j, no imaging. i and j
  * should be actual indices into the coord array (i.e. atom# * 3).
  */
double Frame::COORDDIST2(int i, int j) {
  double x,y,z,D;

  x = X_[i  ] - X_[j  ];
  y = X_[i+1] - X_[j+1];
  z = X_[i+2] - X_[j+2];

  x=x*x;
  y=y*y;
  z=z*z;

  D = x + y + z;
  return D;
}

// Frame::COORDVECTOR()
/** Calculate the vector from atom j to atom i. i and j should be
  * actual indices into the coord array (i.e. atom# * 3).
  */
void Frame::COORDVECTOR(double *V, int i, int j) {
  double *U = X_ + i;
  double *W = X_ + j;
  V[0] = *U - *W;
  ++U;
  ++W; 
  V[1] = *U - *W;
  ++U;
  ++W;
  V[2] = *U - *W;
}

// Frame::ANGLE()
/** Return the angle (in radians) between atoms in M1, M2, M3.
  * Adapted from PTRAJ.
  */
double Frame::ANGLE(AtomMask *M1, AtomMask *M2, AtomMask *M3,bool useMass) {
  double a1[3],a2[3],a3[3];
  double angle, xij, yij, zij, xkj, ykj, zkj, rij, rkj;
  
  if (useMass) {
    CenterOfMass(M1,a1);
    CenterOfMass(M2,a2);
    CenterOfMass(M3,a3);
  } else {
    GeometricCenter(M1,a1);
    GeometricCenter(M2,a2);
    GeometricCenter(M3,a3);
  }

  xij = a1[0] - a2[0];
  yij = a1[1] - a2[1];
  zij = a1[2] - a2[2];

  xkj = a3[0] - a2[0];
  ykj = a3[1] - a2[1]; 
  zkj = a3[2] - a2[2];

  rij = xij*xij + yij*yij + zij*zij;
  rkj = xkj*xkj + ykj*ykj + zkj*zkj;

  if (rij > SMALL && rkj > SMALL) {
    angle = (xij*xkj + yij*ykj + zij*zkj) / sqrt(rij*rkj);
    if (angle > 1.0)
      angle = 1.0;
    else if (angle < -1.0)
      angle = -1.0;
    angle = acos(angle);
  } else
    angle = 0.0;

  return angle;
}

// Frame::ANGLE()
/** Return the angle (in radians) between atoms specified by A1, A2, A3.
  * Adapted from PTRAJ.
  */
double Frame::ANGLE(int A1, int A2, int A3) {
  double angle, xij, yij, zij, xkj, ykj, zkj, rij, rkj;
  int a1, a2, a3;

  a1 = A1 * 3;
  a2 = A2 * 3;
  a3 = A3 * 3;
  xij = X_[a1  ] - X_[a2  ];
  yij = X_[a1+1] - X_[a2+1];
  zij = X_[a1+2] - X_[a2+2];

  xkj = X_[a3  ] - X_[a2  ];
  ykj = X_[a3+1] - X_[a2+1];
  zkj = X_[a3+2] - X_[a2+2];

  rij = xij*xij + yij*yij + zij*zij;
  rkj = xkj*xkj + ykj*ykj + zkj*zkj;

  if (rij > SMALL && rkj > SMALL) {
    angle = (xij*xkj + yij*ykj + zij*zkj) / sqrt(rij*rkj);
    if (angle > 1.0)
      angle = 1.0;
    else if (angle < -1.0)
      angle = -1.0;
    angle = acos(angle);
  } else
    angle = 0.0;

  return angle;
}

// Frame::DIHEDRAL()
/** Return dihedral angle between COM of atoms in M1-M4.
  * NOTE: Torsion returns angles in radians.
  */
double Frame::DIHEDRAL(AtomMask *M1, AtomMask *M2, AtomMask *M3, AtomMask *M4,
                       bool useMass) 
{
  double a1[3],a2[3],a3[3],a4[3];

  if (useMass) {
    CenterOfMass(M1,a1); 
    CenterOfMass(M2,a2); 
    CenterOfMass(M3,a3); 
    CenterOfMass(M4,a4); 
  } else {
    GeometricCenter(M1,a1);
    GeometricCenter(M2,a2);
    GeometricCenter(M3,a3);
    GeometricCenter(M4,a4);
  }

  return Torsion(a1,a2,a3,a4);
}

// Frame::DIHEDRAL()
/** Return dihedral angle between atoms A1-A4.
  * NOTE: Torsion returns angles in radians.
  */
double Frame::DIHEDRAL(int A1, int A2, int A3, int A4) {
  int a1,a2,a3,a4;
  
  a1 = A1 * 3;
  a2 = A2 * 3;
  a3 = A3 * 3;
  a4 = A4 * 3;
  return Torsion(X_+a1,X_+a2,X_+a3,X_+a4);
}

// Frame::PUCKER()
/** Return the pseudorotation between atoms in masks M1-M5 for the given
  * puckerMethod:
  *   0: Use Altona & Sundaralingam method/conventions
  *   1: Use Cremer & Pople method
  * If amplitude is true, return amplitude instead of pseudorotation.
  * NOTE: Pucker routines return angles in radians.
  */
double Frame::PUCKER(AtomMask *M1, AtomMask *M2, AtomMask *M3, AtomMask *M4, AtomMask *M5,
                     int puckerMethod, bool amplitude, bool useMassIn) 
{
  double a1[3],a2[3],a3[3],a4[3],a5[3]; 
  double angle, amp;

  if (useMassIn) {
    CenterOfMass(M1,a1);
    CenterOfMass(M2,a2);
    CenterOfMass(M3,a3);
    CenterOfMass(M4,a4);
    CenterOfMass(M5,a5);
  } else {
    GeometricCenter(M1,a1);
    GeometricCenter(M2,a2);
    GeometricCenter(M3,a3);
    GeometricCenter(M4,a4);
    GeometricCenter(M5,a5);
  }

  angle = 0.0;
  amp = 0.0;
  switch (puckerMethod) {
    case 0 : angle = Pucker_AS(a1,a2,a3,a4,a5,&amp); break;
    case 1 : angle = Pucker_CP(a1,a2,a3,a4,a5,&amp); break;
  }

  if (amplitude) return amp;
  return angle;
}

// Frame::RADGYR()
/** Return the radius of gyration of atoms in mask. Also set the maximum 
  * distance from center. Use center of mass if useMassIn is true.
  */
double Frame::RADGYR(AtomMask *Mask, bool useMassIn, double *max) {
  double mid[3], Coord[3];
  double currentMass, total_mass, maxMass, dist2, sumDist2;

  total_mass=0.0;
  maxMass=1.0;
  currentMass=1.0;
  sumDist2=0.0;
  *max=0.0;

  if (useMassIn)
    total_mass = this->CenterOfMass(Mask,mid);
  else
    total_mass = this->GeometricCenter(Mask,mid);

  for (AtomMask::const_iterator atom = Mask->begin();
                                atom != Mask->end();
                                atom++)
  {
      int atmidx = (*atom) * 3;
      if (useMassIn) 
        currentMass = Mass_[*atom];
      Coord[0] = X_[atmidx  ] - mid[0];
      Coord[1] = X_[atmidx+1] - mid[1];
      Coord[2] = X_[atmidx+2] - mid[2];
      dist2 = (Coord[0]*Coord[0]) + (Coord[1]*Coord[1]) + (Coord[2]*Coord[2]);
      dist2 *= currentMass;
      if (dist2 > *max) {
        *max = dist2;
        maxMass = currentMass;
      }
      sumDist2 += dist2;
  }

  // NOTE: Not using == since it is unreliable for floating point numbers.
  // Should NEVER have a mass smaller than SMALL (vectormath.h)
  if (total_mass < SMALL) return 0;

  currentMass = sqrt(sumDist2 / total_mass); // Radius of Gyration
  *max = sqrt(*max / maxMass);

  return currentMass;
}

// Frame::RMSD()
/** Get the RMSD of this Frame to Ref Frame. Ref frame must contain the same
  * number of atoms as this Frame - should be checked for before this routine
  * is called. Put the best-fit rotation matrix in U and the COM translation 
  * vectors in Trans. The translation is composed of two XYZ vectors; the first
  * is the shift of the XYZ coords to origin, and the second is the shift to Ref 
  * origin. To reproduce the fit perform the first translation (Trans[0...2]), 
  * then rotate (U), then the second translation (Trans[3...5]).
  * Adapted from PTRAJ.
  */
double Frame::RMSD( Frame *Ref, double *U, double *Trans, bool useMassIn) {
  double frameCOM[3], refCOM[3], rms_return, total_mass;
  double mwss, rot[9], rtr[9];
  double xt,yt,zt,xr,yr,zr;
  double *Evector[3], Eigenvalue[3], Emat[9];
  double b[9];
  double cp[3], sig3;
  int i, m;
  int i3,k,k3,j;

  U[0]=0.0; U[1]=0.0; U[2]=0.0;
  U[3]=0.0; U[4]=0.0; U[5]=0.0;
  U[6]=0.0; U[7]=0.0; U[8]=0.0;
  // Trans is set at the end
 
  // Rotation will occur around geometric center/center of mass
  if (useMassIn) {
    total_mass = this->CenterOfMass(frameCOM,0,natom_);
    Ref->CenterOfMass(refCOM,0,natom_);
  } else {
    total_mass = this->GeometricCenter(frameCOM,0,natom_);
    Ref->GeometricCenter(refCOM,0,natom_);
  }
  if (total_mass<SMALL) {
    mprinterr("Error: Frame::RMSD: Divide by zero.\n");
    return -1;
  }
  //fprintf(stderr,"  FRAME COM: %lf %lf %lf\n",frameCOM[0],frameCOM[1],frameCOM[2]); //DEBUG
  //fprintf(stderr,"  REF   COM: %lf %lf %lf\n",refCOM[0],refCOM[1],refCOM[2]); //DEBUG

  // Shift to common COM
  frameCOM[0]=-frameCOM[0]; frameCOM[1]=-frameCOM[1]; frameCOM[2]=-frameCOM[2];
  refCOM[0]  =-refCOM[0];   refCOM[1]  =-refCOM[1];   refCOM[2]  =-refCOM[2];
  this->Translate(frameCOM);
  Ref->Translate(refCOM);
  //mprintf("  SHIFTED FRAME 0: %lf %lf %lf\n",X[0],X[1],X[2]); //DEBUG
  //mprintf("  SHIFTED REF 0  : %lf %lf %lf\n",Ref->X[0],Ref->X[1],Ref->X[2]); //DEBUG

  // Use Kabsch algorithm to calculate optimum rotation matrix.
  // U = [(RtR)^.5][R^-1]
  mwss=0.0;
  rot[0]=0.0; rot[1]=0.0; rot[2]=0.0;
  rot[3]=0.0; rot[4]=0.0; rot[5]=0.0;
  rot[6]=0.0; rot[7]=0.0; rot[8]=0.0;
  // rtr is set below
  // Calculate covariance matrix of Coords and Reference (R = Xt * Ref)
  m=0;
  for (i=0; i<Ncoord_; i+=3) {
    xt = X_[i];
    yt = X_[i+1];
    zt = X_[i+2];
    xr = Ref->X_[i];
    yr = Ref->X_[i+1];
    zr = Ref->X_[i+2];

    // Use rms_return to hold mass for this atom if specified
    rms_return = 1.0;
    if (useMassIn) 
      rms_return = Mass_[m++];
    //total_mass+=rms_return;

    mwss += rms_return * ( (xt*xt)+(yt*yt)+(zt*zt)+(xr*xr)+(yr*yr)+(zr*zr) );

    // Calculate the Kabsch matrix: R = (rij) = Sum(yni*xnj)
    rot[0] += rms_return*xt*xr;
    rot[1] += rms_return*xt*yr;
    rot[2] += rms_return*xt*zr;

    rot[3] += rms_return*yt*xr;
    rot[4] += rms_return*yt*yr;
    rot[5] += rms_return*yt*zr;

    rot[6] += rms_return*zt*xr;
    rot[7] += rms_return*zt*yr;
    rot[8] += rms_return*zt*zr;
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
  if (!diagEsort(rtr, Emat, Evector, Eigenvalue))
    return(0);

  // a3 = a1 x a2 
  CROSS_PRODUCT(Evector[2][0], Evector[2][1], Evector[2][2],
                Evector[0][0], Evector[0][1], Evector[0][2],
                Evector[1][0], Evector[1][1], Evector[1][2]);

  // Evector dot transpose rot: b = R . ak 
  b[0] = Evector[0][0]*rot[0] + Evector[0][1]*rot[3] + Evector[0][2]*rot[6];
  b[1] = Evector[0][0]*rot[1] + Evector[0][1]*rot[4] + Evector[0][2]*rot[7];
  b[2] = Evector[0][0]*rot[2] + Evector[0][1]*rot[5] + Evector[0][2]*rot[8];
  normalize(&b[0]);
  b[3] = Evector[1][0]*rot[0] + Evector[1][1]*rot[3] + Evector[1][2]*rot[6];
  b[4] = Evector[1][0]*rot[1] + Evector[1][1]*rot[4] + Evector[1][2]*rot[7];
  b[5] = Evector[1][0]*rot[2] + Evector[1][1]*rot[5] + Evector[1][2]*rot[8];
  normalize(&b[3]);
  b[6] = Evector[2][0]*rot[0] + Evector[2][1]*rot[3] + Evector[2][2]*rot[6];
  b[7] = Evector[2][0]*rot[1] + Evector[2][1]*rot[4] + Evector[2][2]*rot[7];
  b[8] = Evector[2][0]*rot[2] + Evector[2][1]*rot[5] + Evector[2][2]*rot[8];
  normalize(&b[6]);

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
  for (k=k3=0; k<3; k++,k3+=3)
    for (i=i3=0; i<3; i++,i3+=3)
      for (j=0; j<3; j++)
        U[i3+j] += Evector[k][j] * b[k3+i]; 

  // E=E0-sqrt(mu1)-sqrt(mu2)-sig3*sqrt(mu3) 
  rms_return = mwss;
  rms_return -= sqrt(fabs(Eigenvalue[0]));
  rms_return -= sqrt(fabs(Eigenvalue[1]));
  rms_return -= (sig3*sqrt(fabs(Eigenvalue[2])));

  if (rms_return<0) {
    //fprintf(stderr,"RMS returned is <0 before sqrt, setting to 0 (%lf)\n",rms_return);
    rms_return=0.0;
  }
  else
    rms_return = sqrt((2.0*rms_return)/total_mass);

  /* Translation matrix - coords are shifted to common CoM first (origin), then
   * to original reference location.
   * Remember frameCOM and refCOM were negated above to facilitate translation to COM.
   */
  Trans[0] = frameCOM[0];
  Trans[1] = frameCOM[1];
  Trans[2] = frameCOM[2];
  Trans[3] = -refCOM[0];
  Trans[4] = -refCOM[1];
  Trans[5] = -refCOM[2];

  //DEBUG
  //printRotTransInfo(U,Trans);
  //fprintf(stdout,"RMS is %lf\n",rms_return);

  return rms_return;
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
double Frame::RMSD_CenteredRef( Frame &Ref, double U[9], double Trans[6], bool useMassIn) 
{
  double frameCOM[3], rms_return, total_mass, atom_mass;
  double mwss, rot[9], rtr[9];
  double xt,yt,zt,xr,yr,zr;
  double *Evector[3], Eigenvalue[3], Emat[9];
  double b[9];
  double cp[3], sig3;
  int i, m;
  int i3,k,k3,j;

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
  frameCOM[0]=-frameCOM[0]; frameCOM[1]=-frameCOM[1]; frameCOM[2]=-frameCOM[2];
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
  m=0;
  for (i=0; i<Ncoord_; i+=3) {
    xt = X_[i];
    yt = X_[i+1];
    zt = X_[i+2];
    xr = Ref.X_[i];
    yr = Ref.X_[i+1];
    zr = Ref.X_[i+2];

    // Use atom_mass to hold mass for this atom if specified
    atom_mass = 1.0;
    if (useMassIn) 
      atom_mass = Mass_[m++];
    //total_mass+=atom_mass;

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
  if (!diagEsort(rtr, Emat, Evector, Eigenvalue))
    return(0);

  // a3 = a1 x a2 
  CROSS_PRODUCT(Evector[2][0], Evector[2][1], Evector[2][2],
                Evector[0][0], Evector[0][1], Evector[0][2],
                Evector[1][0], Evector[1][1], Evector[1][2]);

  // Evector dot transpose rot: b = R . ak 
  b[0] = Evector[0][0]*rot[0] + Evector[0][1]*rot[3] + Evector[0][2]*rot[6];
  b[1] = Evector[0][0]*rot[1] + Evector[0][1]*rot[4] + Evector[0][2]*rot[7];
  b[2] = Evector[0][0]*rot[2] + Evector[0][1]*rot[5] + Evector[0][2]*rot[8];
  normalize(&b[0]);
  b[3] = Evector[1][0]*rot[0] + Evector[1][1]*rot[3] + Evector[1][2]*rot[6];
  b[4] = Evector[1][0]*rot[1] + Evector[1][1]*rot[4] + Evector[1][2]*rot[7];
  b[5] = Evector[1][0]*rot[2] + Evector[1][1]*rot[5] + Evector[1][2]*rot[8];
  normalize(&b[3]);
  b[6] = Evector[2][0]*rot[0] + Evector[2][1]*rot[3] + Evector[2][2]*rot[6];
  b[7] = Evector[2][0]*rot[1] + Evector[2][1]*rot[4] + Evector[2][2]*rot[7];
  b[8] = Evector[2][0]*rot[2] + Evector[2][1]*rot[5] + Evector[2][2]*rot[8];
  normalize(&b[6]);

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
  for (k=k3=0; k<3; k++,k3+=3)
    for (i=i3=0; i<3; i++,i3+=3)
      for (j=0; j<3; j++)
        U[i3+j] += Evector[k][j] * b[k3+i]; 

  // E=E0-sqrt(mu1)-sqrt(mu2)-sig3*sqrt(mu3) 
  rms_return = mwss;
  rms_return -= sqrt(fabs(Eigenvalue[0]));
  rms_return -= sqrt(fabs(Eigenvalue[1]));
  rms_return -= (sig3*sqrt(fabs(Eigenvalue[2])));

  if (rms_return<0) {
    //fprintf(stderr,"RMS returned is <0 before sqrt, setting to 0 (%lf)\n",rms_return);
    rms_return=0.0;
  }
  else
    rms_return = sqrt((2.0*rms_return)/total_mass);

  /* Translation matrix - coords are shifted to common CoM first (origin), then
   * to original reference location.
   * Remember frameCOM was negated above to facilitate translation to COM.
   * Reference translation should already be set
   */
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
double Frame::RMSD( Frame *Ref, bool useMass ) {
  double rms_return, total_mass, xx,yy,zz, currentMass;
  int i, m;

  rms_return = 0.0;
  total_mass = 0.0;
  currentMass = 1.0;
  m=0;

  for (i=0; i < Ncoord_; i+=3) {
    if (useMass) 
      currentMass = Mass_[m++];
    total_mass += currentMass;
    xx = Ref->X_[i]   - X_[i];
    yy = Ref->X_[i+1] - X_[i+1];
    zz = Ref->X_[i+2] - X_[i+2];
    rms_return += currentMass * (xx*xx + yy*yy + zz*zz);
  }
  if (total_mass<SMALL) {
    mprinterr("Error: Frame::RMSD: Divide by zero.\n");
    return -1;
  }
  if (rms_return < 0) {
    //mprinterr("Error: Frame::RMSD: Negative RMS. Coordinates may be corrupted.\n");
    //return -1;
    //mprinterr("RMS returned is <0 before sqrt, setting to 0 (%lf)\n",rms_return);
    return 0;
  }
  rms_return = sqrt(rms_return / total_mass);

  return rms_return;
}

// Frame::DISTRMSD()
/** Calcuate the distance RMSD of Frame to Ref. Frames must contain
  * same # of atoms. Should not be called for 0 atoms.
  */
double Frame::DISTRMSD( Frame *Ref ) {
  double TgtDist, RefDist;
  double diff, rms_return;
  double x,y,z;
  int a10,a11,a12;
  int a20,a21,a22; 
  double Ndistances = ((natom_ * natom_) - natom_) / 2;

  rms_return = 0;
  a10 = 0;
  for (int atom1 = 0; atom1 < natom_-1; atom1++) {
    a20 = a10 + 3;
    for (int atom2 = atom1+1; atom2 < natom_; atom2++) {
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
      x = Ref->X_[a10] - Ref->X_[a20];
      x = x * x;
      y = Ref->X_[a11] - Ref->X_[a21];
      y = y * y;
      z = Ref->X_[a12] - Ref->X_[a22];
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
void Frame::RotateAroundAxis(double *T, double theta, AtomMask &Rmask) {
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

