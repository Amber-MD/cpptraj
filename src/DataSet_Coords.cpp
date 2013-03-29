#include "DataSet_Coords.h"
#include "CpptrajStdio.h"

DataSet_Coords::DataSet_Coords() :
  DataSet_1D(COORDS, 8, 3),
  numBoxCrd_(0)
{}

int DataSet_Coords::Allocate1D( size_t sizeIn ) {
  coords_.reserve( sizeIn );
  return 0;
}

void DataSet_Coords::Info() const {
  top_.ParmInfo();
}

void DataSet_Coords::SetTopology(Topology const& topIn) {
  top_ = topIn;
  if (top_.ParmBox().HasBox())
    numBoxCrd_ = 6;
  else
    numBoxCrd_ = 0;
}
