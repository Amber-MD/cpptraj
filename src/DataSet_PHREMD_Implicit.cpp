#include "DataSet_PHREMD_Implicit.h"
#include "CpptrajStdio.h"

/// CONSTRUCTOR
DataSet_PHREMD_Implicit::DataSet_PHREMD_Implicit() :
  DataSet_PHREMD(PH_IMPL, TextFormat(TextFormat::DOUBLE, 10, 4))
{}

size_t DataSet_PHREMD_Implicit::MemUsageInBytes() const {
  size_t mySize = 0;
  for (RecArray::const_iterator it = records_.begin(); it != records_.end(); ++it)
    mySize += it->DataSize();
  return mySize;
}

int DataSet_PHREMD_Implicit::Allocate(SizeArray const& sizeIn) {
  if (!sizeIn.empty()) records_.reserve( sizeIn[0] ); 
  return 0;
}

void DataSet_PHREMD_Implicit::Info() const {
  mprintf(" (%zu residues, %zu frames)", residues_.size(), records_.size());
}
