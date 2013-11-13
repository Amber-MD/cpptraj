#include "DataSet_Coords.h"

Frame DataSet_Coords::AllocateFrame() const {
  Frame f;
  f.SetupFrameV( top_.Atoms(), (numVel_ > 0), 0 );
  return f;
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
