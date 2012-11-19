#include "DataSet_Coords.h"
#include "CpptrajStdio.h"

DataSet_Coords::DataSet_Coords() :
  DataSet(COORDS, 8, 3, 4),
  top_(0)
{}

DataSet_Coords::~DataSet_Coords() {
  if (top_ != 0) delete top_;
}

void DataSet_Coords::AddFrame(Frame const& frameIn) {
  coords_.push_back( frameIn.ConvertToFloat(mask_) );
}

void DataSet_Coords::SetFrame(int idx, Frame const& frameIn) {
  coords_[idx] = frameIn.ConvertToFloat();
}

int DataSet_Coords::Allocate( int sizeIn ) {
  coords_.reserve( sizeIn );
  return 0;
}

void DataSet_Coords::Info() {
  //mprintf(", mask [%s] (%i atoms)", mask_.MaskString(), Natom());
  if (top_ != 0) top_->ParmInfo();
}

/** \return 0 on success, 1 on error, -1 if topology/mask already set up. */
int DataSet_Coords::SetupTopMask(std::string const& maskexpr, Topology& topIn) {
  if (top_ != 0) return -1;
  mask_.SetMaskString( maskexpr );
  if (topIn.SetupIntegerMask( mask_ )) return 1;
  top_ = topIn.modifyStateByMask( mask_ );
  if (top_ == 0) return 1;
  return 0;
}

void DataSet_Coords::SetTop(Topology* topIn) {
  if (top_ == 0) top_ = new Topology();
  *top_ = *topIn;
}
