#include "CompactFrameArray.h"
#include "CpptrajStdio.h"
#include <algorithm> // std::fill

/** CONSTRUCTOR */
CompactFrameArray::CompactFrameArray() :
  currentIdx_(-1),
  maxIdx_(0)
{
  std::fill(componentIdx_, componentIdx_+CoordinateInfo::NCOMPONENTS, -1);
}

/** COPY CONSTRUCTOR */
CompactFrameArray::CompactFrameArray(CompactFrameArray const& rhs) :
  compactFrames_(rhs.compactFrames_),
  components_(rhs.components_),
  offsets_(rhs.offsets_),
  currentIdx_(rhs.currentIdx_),
  maxIdx_(rhs.maxIdx_)
{
  std::copy(rhs.componentIdx_, rhs.componentIdx_+CoordinateInfo::NCOMPONENTS, componentIdx_);
}

/** ASSIGNMENT OPERATOR */
CompactFrameArray& CompactFrameArray::operator=(CompactFrameArray const& rhs) {
  if (this == &rhs) return *this;
  compactFrames_ = rhs.compactFrames_;
  std::copy(rhs.componentIdx_, rhs.componentIdx_+CoordinateInfo::NCOMPONENTS, componentIdx_);
  components_ = rhs.components_;
  offsets_ = rhs.offsets_;
  currentIdx_ = rhs.currentIdx_;
  maxIdx_ = rhs.maxIdx_;
  return *this;
}

/** Allocate for specified number of frames. */
void CompactFrameArray::Resize(int nframes) {
  //rprintf("DEBUG: Calling CompactFrameArray::Resize for %i frames.\n", nframes);
  if (nframes > 0 && !offsets_.empty()) {
    compactFrames_.resize( offsets_.back() * nframes );
    maxIdx_ = nframes;
    if (currentIdx_ >= maxIdx_)
      currentIdx_ = maxIdx_ - 1;
  }
}

/** \return Size of a single frame in elements. */
unsigned int CompactFrameArray::FrameSize() const {
  if (offsets_.empty())
    return 0;
  else
    return (unsigned int)offsets_.back();
}

/** \return Total size of the array in bytes. */
unsigned int CompactFrameArray::SizeInBytes() const {
  return ( compactFrames_.size()       * sizeof(float) +
           CoordinateInfo::NCOMPONENTS * sizeof(int)   +
           components_.size()          * sizeof(CoordinateInfo::Component) +
           offsets_.size()             * sizeof(long int) +
           2                           * sizeof(int) // currentIdx_ and maxIdx_
         );
}

/** Compare components and offsets */
bool CompactFrameArray::operator!=(CompactFrameArray const& rhs) const {
  if (components_.size() != rhs.components_.size() ||
      offsets_.size()    != rhs.offsets_.size())
    return true;
  for (unsigned int cidx = 0; cidx != components_.size(); cidx++)
  {
    if (components_[cidx] != rhs.components_[cidx]) return true;
    if (offsets_[cidx]    != rhs.offsets_[cidx]   ) return true;
  }
  if (offsets_.size() > 0 && offsets_.back() != rhs.offsets_.back())
    return true;
  return false;
}

/** Add component of specified size to array, update the offset. */
void CompactFrameArray::addComponent(long int& currentOffset, CoordinateInfo::Component compIn,
                                     long int sizeIn)
{
  componentIdx_[compIn] = components_.size();
  components_.push_back( compIn );
  offsets_.push_back( currentOffset );
  currentOffset += sizeIn;
}

/** \return Size of a frame containing given coordinate info in bytes. */
unsigned int CompactFrameArray::EstimateFrameSizeInBytes(CoordinateInfo const& cinfoIn, unsigned int natoms)
{
  unsigned int frameSize = 0;
  if (cinfoIn.HasCrd())         frameSize += natoms * 3;
  if (cinfoIn.HasVel())         frameSize += natoms * 3;
  if (cinfoIn.HasForce())       frameSize += natoms * 3;
  if (cinfoIn.HasBox())         frameSize += 9;
  if (cinfoIn.HasTemp())        frameSize += 1;
  if (cinfoIn.Has_pH())         frameSize += 1;
  if (cinfoIn.HasRedOx())       frameSize += 1;
  if (cinfoIn.HasTime())        frameSize += 1;
  if (cinfoIn.HasStep())        frameSize += 1;
  if (cinfoIn.HasReplicaDims()) frameSize += cinfoIn.ReplicaDimensions().Ndims();
  if (cinfoIn.HasRepIdx())      frameSize += 1;
  if (cinfoIn.HasCrdIdx())      frameSize += 1;

  return frameSize * sizeof(float);
}

/** Set up frame array. */
int CompactFrameArray::SetupFrameArray(CoordinateInfo const& cinfoIn, unsigned int natoms, int nframes)
{
  std::fill(componentIdx_, componentIdx_+CoordinateInfo::NCOMPONENTS, -1);
  components_.clear();
  offsets_.clear();
  long int currentOffset = 0;

  if (cinfoIn.HasCrd())
    addComponent(currentOffset, CoordinateInfo::POSITION, (long int)natoms * 3);
  if (cinfoIn.HasVel())
    addComponent(currentOffset, CoordinateInfo::VELOCITY, (long int)natoms * 3);
  if (cinfoIn.HasForce())
    addComponent(currentOffset, CoordinateInfo::FORCE, (long int)natoms * 3);
  if (cinfoIn.HasBox())
    addComponent(currentOffset, CoordinateInfo::BOX, 9);
  if (cinfoIn.HasTemp())
    addComponent(currentOffset, CoordinateInfo::TEMPERATURE, 1);
  if (cinfoIn.Has_pH())
    addComponent(currentOffset, CoordinateInfo::PH, 1);
  if (cinfoIn.HasRedOx())
    addComponent(currentOffset, CoordinateInfo::REDOX, 1);
  if (cinfoIn.HasTime())
    addComponent(currentOffset, CoordinateInfo::TIME, 1);
  if (cinfoIn.HasStep())
    addComponent(currentOffset, CoordinateInfo::STEP, 1);
  if (cinfoIn.HasReplicaDims())
    addComponent(currentOffset, CoordinateInfo::REMD_INDICES, cinfoIn.ReplicaDimensions().Ndims());
  if (cinfoIn.HasRepIdx())
    addComponent(currentOffset, CoordinateInfo::REPIDX, 1);
  if (cinfoIn.HasCrdIdx())
    addComponent(currentOffset, CoordinateInfo::CRDIDX, 1);

  // Final "offset" is the total frame size
  offsets_.push_back( currentOffset );

  // Allocate for specified number of frames
  Resize(nframes);
  currentIdx_ = -1;
/*
  rprintf("DEBUG: CompactFrameArray (%u frames) has the following components:\n", maxIdx_);
  for (unsigned int cidx = 0; cidx != components_.size(); cidx++)
    rprintf("DEBUG:\t%20s : %li - %li\n", CoordinateInfo::ComponentStr(components_[cidx]),
            offsets_[cidx], offsets_[cidx+1]);
*/
  return 0;
}

// -----------------------------------------------------------------------------
/// \return 1 and print a error message that the given component is not in the array.
static inline int ComponentNotFoundErr(CoordinateInfo::Component cmpt) {
  mprinterr("Error: Component '%s' not present.\n", CoordinateInfo::ComponentStr(cmpt));
  return 1;
}

/** Seek to frame, allocate if necessary. */
void CompactFrameArray::SeekAndAllocate(unsigned int idx) {
  if ((int)idx >= maxIdx_) {
    compactFrames_.resize( (idx+1) * offsets_.back() );
    maxIdx_ = idx+1;
  }
  currentIdx_ = idx;
  //rprintf("DBG: Called SeekAndAllocate(%u), currentIdx=%i maxIdx=%i\n", idx, currentIdx_, maxIdx_);
  //unsigned int frameBeginIdx = idx * offsets_.back();
  //if (frameBeginIdx >= compactFrames_.size()) {
  // }
  //frameBegin_ = &compactFrames_[0] + frameBeginIdx;
}

/** Advance to next frame, allocte if necessary. */
void CompactFrameArray::NextAndAllocate() {
  SeekAndAllocate( currentIdx_+1 );
}

/** Set component at current frame from given double pointer. */
int CompactFrameArray::SetFromDblPtr(const double* ptrIn, CoordinateInfo::Component cmpt)
{
  int cidx = componentIdx_[cmpt];
  if (cidx < 0) return ComponentNotFoundErr(cmpt);
  float* frameBegin = (&compactFrames_[0]) + (currentIdx_ * offsets_.back());
  const double* ptr = ptrIn;
  for (long int i = offsets_[cidx]; i != offsets_[cidx+1]; ++i, ++ptr)
    frameBegin[i] = (float)(*ptr);
  return 0;
}

/** Set component at current frame from given integer pointer. */
int CompactFrameArray::SetFromIntPtr(const int* ptrIn, CoordinateInfo::Component cmpt)
{
  int cidx = componentIdx_[cmpt];
  if (cidx < 0) return ComponentNotFoundErr(cmpt);
  float* frameBegin = (&compactFrames_[0]) + (currentIdx_ * offsets_.back());
  const int* ptr = ptrIn;
  for (long int i = offsets_[cidx]; i != offsets_[cidx+1]; ++i, ++ptr)
    frameBegin[i] = (float)(*ptr);
  return 0;
}

/** Set single value component at current frame. */
int CompactFrameArray::SetFromDblVal(double dval, CoordinateInfo::Component cmpt)
{
  int cidx = componentIdx_[cmpt];
  if (cidx < 0) return ComponentNotFoundErr(cmpt);
  float* frameBegin = (&compactFrames_[0]) + (currentIdx_ * offsets_.back());
  frameBegin[offsets_[cidx]] = (float)dval;
  return 0;
}

/** Set single value component at current frame. */
int CompactFrameArray::SetFromIntVal(int ival, CoordinateInfo::Component cmpt)
{
  int cidx = componentIdx_[cmpt];
  if (cidx < 0) return ComponentNotFoundErr(cmpt);
  float* frameBegin = (&compactFrames_[0]) + (currentIdx_ * offsets_.back());
  frameBegin[offsets_[cidx]] = (float)ival;
  return 0;
}
// -----------------------------------------------------------------------------
/** Copy component to given double pointer. */
int CompactFrameArray::GetToDblPtr(double* ptrOut, unsigned int idx, CoordinateInfo::Component cmpt)
const
{
  int cidx = componentIdx_[cmpt];
  if (cidx < 0) return ComponentNotFoundErr(cmpt);
  const float* frameBegin = (&compactFrames_[0]) + (idx * offsets_.back());
  double* ptr = ptrOut;
  for (long int i = offsets_[cidx]; i != offsets_[cidx+1]; ++i, ++ptr)
    *ptr = (double)(frameBegin[i]);
  return 0;
}

/** Copy component to given double pointer with mask. For use with coords,
  * velocities, forces.
  */
int CompactFrameArray::GetToMaskDblPtr(double* ptrOut, std::vector<int> const& selected, unsigned int idx,
                                       CoordinateInfo::Component cmpt)
const
{
  int cidx = componentIdx_[cmpt];
  if (cidx < 0) return ComponentNotFoundErr(cmpt);
  const float* frameBegin = (&compactFrames_[0]) + (idx * offsets_.back());
  // Loop over atom indices in the mask.
  for (unsigned int j = 0; j != selected.size(); j++)
  {
    // Used for coordinate index in compactFrames_
    int          frmIdx = 3*selected[j];
    // Used for coordinate index in ptrOut (Frame)
    unsigned int ptrIdx = 3*j;
    ptrOut[ptrIdx  ] = (double)(frameBegin[offsets_[cidx] + frmIdx  ]);
    ptrOut[ptrIdx+1] = (double)(frameBegin[offsets_[cidx] + frmIdx+1]);
    ptrOut[ptrIdx+2] = (double)(frameBegin[offsets_[cidx] + frmIdx+2]);
  }
  //  ptrOut[j] = (double)(frameBegin[offsets_[cidx] + (3*selected[j])]);

  return 0;
}

/** Copy component to given integer pointer. */
int CompactFrameArray::GetToIntPtr(int* ptrOut, unsigned int idx, CoordinateInfo::Component cmpt)
const
{
  int cidx = componentIdx_[cmpt];
  if (cidx < 0) return ComponentNotFoundErr(cmpt);
  const float* frameBegin = (&compactFrames_[0]) + (idx * offsets_.back());
  int* ptr = ptrOut;
  for (long int i = offsets_[cidx]; i != offsets_[cidx+1]; ++i, ++ptr)
    *ptr = (int)(frameBegin[i]);
  return 0;
}

/** \return Single value component. */
float CompactFrameArray::GetVal(unsigned int idx, CoordinateInfo::Component cmpt)
const
{
  int cidx = componentIdx_[cmpt];
  if (cidx < 0) return ComponentNotFoundErr(cmpt);
  const float* frameBegin = (&compactFrames_[0]) + (idx * offsets_.back());
  return frameBegin[offsets_[cidx]];
}
