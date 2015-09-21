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
