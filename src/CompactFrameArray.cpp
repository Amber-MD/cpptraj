#include "CompactFrameArray.h"
#include "CpptrajStdio.h"
#include <algorithm> // std::fill

/** CONSTRUCTOR */
CompactFrameArray::CompactFrameArray() //:
//  frameSize_(0)
{
  std::fill(componentIdx_, componentIdx_+CoordinateInfo::NCOMPONENTS, -1);
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

/** Allocate for specified number of frames. */
void CompactFrameArray::Resize(int nframes) {
  if (nframes > 0 && !offsets_.empty())
    compactFrames_.resize( offsets_.back() * nframes );
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

  mprintf("DEBUG: CompactFrameArray (%i frames) has the following components:\n", nframes);
  for (unsigned int cidx = 0; cidx != components_.size(); cidx++)
    mprintf("\t%20s : %li - %li\n", CoordinateInfo::ComponentStr(components_[cidx]),
            offsets_[cidx], offsets_[cidx+1]);

  return 0;
}

static inline int ComponentNotFoundErr(CoordinateInfo::Component cmpt) {
  mprinterr("Error: Component '%s' not present.\n", CoordinateInfo::ComponentStr(cmpt));
  return 1;
}

int CompactFrameArray::SetFromDblPtr(unsigned int idx, const double* ptrIn, CoordinateInfo::Component cmpt)
{
  int cidx = componentIdx_[cmpt];
  if (cidx < 0) return ComponentNotFoundErr(cmpt);
  float* frameBegin = (&compactFrames_[0]) + (idx * offsets_.back());
  const double* ptr = ptrIn;
  for (long int i = offsets_[cidx]; i != offsets_[cidx+1]; ++i, ++ptr)
    frameBegin[i] = (float)(*ptr);
  return 0;
}

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

int CompactFrameArray::SetFromIntPtr(unsigned int idx, const int* ptrIn, CoordinateInfo::Component cmpt)
{
  int cidx = componentIdx_[cmpt];
  if (cidx < 0) return ComponentNotFoundErr(cmpt);
  float* frameBegin = (&compactFrames_[0]) + (idx * offsets_.back());
  const int* ptr = ptrIn;
  for (long int i = offsets_[cidx]; i != offsets_[cidx+1]; ++i, ++ptr)
    frameBegin[i] = (float)(*ptr);
  return 0;
}

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

int CompactFrameArray::SetFromDblVal(unsigned int idx, double dval, CoordinateInfo::Component cmpt)
{
  int cidx = componentIdx_[cmpt];
  if (cidx < 0) return ComponentNotFoundErr(cmpt);
  float* frameBegin = (&compactFrames_[0]) + (idx * offsets_.back());
  frameBegin[offsets_[cidx]] = (float)dval;
  return 0;
}

float CompactFrameArray::GetVal(unsigned int idx, CoordinateInfo::Component cmpt)
const
{
  int cidx = componentIdx_[cmpt];
  if (cidx < 0) return ComponentNotFoundErr(cmpt);
  const float* frameBegin = (&compactFrames_[0]) + (idx * offsets_.back());
  return frameBegin[offsets_[cidx]];
}
