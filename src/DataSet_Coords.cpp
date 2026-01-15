#include "DataSet_Coords.h"
#include "CpptrajStdio.h"

Frame DataSet_Coords::AllocateFrame() const {
  Frame f;
  f.SetupFrameV( top_.Atoms(), cInfo_ );
  return f;
}

void DataSet_Coords::CommonInfo() const {
  if (cInfo_.HasBox()) mprintf(" Box Coords,");
  if (cInfo_.HasVel()) mprintf(" Velocities,");
  mprintf(" %i atoms", top_.Natom());
}

/** Update coords info box */
void DataSet_Coords::UpdateCoordsInfoBox( Box const& boxIn ) {
  cInfo_.SetBox( boxIn );
}

/** Copy the incoming topology if different from current topology. */
void DataSet_Coords::set_topology(Topology const& topIn)
{
  if (&topIn != &top_) {
    top_ = topIn;
  }
/*
  if (&topIn == &top_) {
    mprintf("DEBUG: set_topology called with existing topology (%s)\n", top_.c_str());
  } else {
    if (top_.Natom() > 0)
      mprintf("DEBUG: set_topology: overwriting top %s with top %s\n", top_.c_str(), topIn.c_str());
    top_ = topIn;
  }
*/
}
