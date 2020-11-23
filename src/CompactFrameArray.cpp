#include "CompactFrameArray.h"
#include "CpptrajStdio.h"
#include <algorithm> // std::fill

/** CONSTRUCTOR */
CompactFrameArray::CompactFrameArray() :
  frameSize_(0)
{
  std::fill(hasComponent_, hasComponent_+CoordinateInfo::NCOMPONENTS, false);
}

/** Add component of specified size to array, update the offset. */
void CompactFrameArray::addComponent(long int& currentOffset, CoordinateInfo::Component compIn,
                                     long int sizeIn)
{
  hasComponent_[compIn] = true;
  components_.push_back( compIn );
  offsets_.push_back( currentOffset );
  currentOffset += sizeIn;
}

/** Set up frame array. */
int CompactFrameArray::SetupFrameArray(CoordinateInfo const& cinfoIn, unsigned int natoms, int nframes)
{
  std::fill(hasComponent_, hasComponent_+CoordinateInfo::NCOMPONENTS, false);
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
  if (nframes > 0)
    compactFrames_.resize( offsets_.back() * nframes );

  return 0;
}

/** \return index of specified component, or -1 if component not present. */
int CompactFrameArray::ComponentIndex(CoordinateInfo::Component cmpt) const {
  for (int i = 0; i < (int)components_.size(); ++i)
    if (components_[i] == cmpt) return i;
  mprinterr("Internal Error: Component not present.\n");
  return -1;
}

int CompactFrameArray::SetFromDblPtr(unsigned int idx, const double* ptrIn, CoordinateInfo::Component cmpt)
{
  int cidx = ComponentIndex(cmpt);
  if (cidx < 0) return 1;
  float* frameBegin = (&compactFrames_[0]) + (idx * offsets_.back());
  const double* ptr = ptrIn;
  for (long int i = offsets_[cidx]; i != offsets_[cidx+1]; ++i, ++ptr)
    frameBegin[i] = (float)(*ptr);
  return 0;
}

int CompactFrameArray::SetFromDblVal(unsigned int idx, double dval, CoordinateInfo::Component cmpt)
{
  int cidx = ComponentIndex(cmpt);
  if (cidx < 0) return 1;
  float* frameBegin = (&compactFrames_[0]) + (idx * offsets_.back());
  frameBegin[offsets_[cidx]] = (float)dval;
  return 0;
}
