#include <algorithm> // std::min, std::max
#include "DataSet_Cmatrix_NOMEM.h"

DataSet_Cmatrix_NOMEM::~DataSet_Cmatrix_NOMEM() {
  if (cdist_ != 0) delete cdist_;
}

void DataSet_Cmatrix_NOMEM::WriteBuffer(CpptrajFile& outfile, SizeArray const& pIn) const {
  int xidx = (int)std::min(pIn[0], pIn[1]);
  int yidx = (int)std::max(pIn[0], pIn[1]);
  if ( xidx >= sievedFrames_.ActualNframes() || yidx >= sievedFrames_.ActualNframes() )
    outfile.Printf(format_.fmt(), 0.0);
  else
    outfile.Printf(format_.fmt(), cdist_->FrameDist( sievedFrames_.IdxToFrame(yidx),
                                                     sievedFrames_.IdxToFrame(xidx) ));
}

int DataSet_Cmatrix_NOMEM::SetCdist(ClusterDist* cdistIn)
{
  if (cdistIn == 0) return 1;
  cdist_ = cdistIn->Copy();
  return 0;
}
