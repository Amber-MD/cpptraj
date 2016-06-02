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
    outfile.Printf(format_.fmt(), cdist_->FrameDist( idxToFrame_[yidx], idxToFrame_[xidx] ));
}

int DataSet_Cmatrix_NOMEM::SetupSieveAndCdist(size_t sizeIn, size_t sieveIn, int iseed,
                                              ClusterDist* cdistIn)
{
  if (cdistIn == 0) return 1;
  if (SetupWithSieve(sizeIn, sieveIn, iseed)) return 1;
  cdist_ = cdistIn->Copy();
  idxToFrame_ = sievedFrames_.Frames();
  return 0;
}
