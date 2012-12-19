#include "DataSet_Coords.h"
#include "CpptrajStdio.h"

DataSet_Coords::DataSet_Coords() :
  DataSet(COORDS, 8, 3, 4),
  numBoxCrd_(0)
{}

int DataSet_Coords::Allocate( int sizeIn ) {
  coords_.reserve( sizeIn );
  return 0;
}

void DataSet_Coords::Info() {
  top_.ParmInfo();
}

void DataSet_Coords::SetTopology(Topology const& topIn) {
  top_ = topIn;
  if (top_.ParmBox().HasBox())
    numBoxCrd_ = 6;
  else
    numBoxCrd_ = 0;
}
