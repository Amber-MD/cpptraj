#include "CompactFrameArray.h"
#include <algorithm> // std::fill

/** CONSTRUCTOR */
CompactFrameArray::CompactFrameArray() :
  frameSize_(0)
{
  std::fill(hasComponent_, hasComponent_+CoordinateInfo::NCOMPONENTS, false);
}

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

  // Final "offset" is the total frame size
  offsets_.push_back( currentOffset );

  return 0;
}
