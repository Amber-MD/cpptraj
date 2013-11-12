#include "DataSet_Coords.h"
#include "CpptrajStdio.h"

DataSet_Coords::DataSet_Coords() :
  DataSet_1D(COORDS, 8, 3),
  numBoxCrd_(0),
  numVel_(0)
{}

int DataSet_Coords::Allocate1D( size_t sizeIn ) {
  coords_.reserve( sizeIn );
  return 0;
}

void DataSet_Coords::Info() const {
  size_t sze = (((coords_.size() * (size_t)top_.Natom() * 3UL) + numBoxCrd_) * sizeof(float))
               / 1048576UL;
  if (sze == 0)
    mprintf(" (<1 MB)");
  else
    mprintf(" (%zu MB)", sze);
  // DEBUG
  mprintf(" box=%zu vel=%zu", numBoxCrd_, numVel_); 
  top_.Brief();
}

void DataSet_Coords::SetTopology(Topology const& topIn) {
  top_ = topIn;
  if (top_.ParmBox().HasBox())
    numBoxCrd_ = 6;
  else
    numBoxCrd_ = 0;
  if (top_.HasVelInfo())
    numVel_ = top_.Natom() * 3;
  else
    numVel_ = 0;
}
