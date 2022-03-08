#ifndef INC_DATAIO_CMATRIX_BINARY_H
#define INC_DATAIO_CMATRIX_BINARY_H
#include "DataIO.h"
#include "DataSet_PairwiseCache_MEM.h" // TODO just pairwisecache
/// Read/write cpptraj-format binary cluster pairwise distance matrix files. 
class DataIO_Cmatrix_Binary : public DataIO {
  public:
    DataIO_Cmatrix_Binary();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_Cmatrix_Binary(); }
    static void ReadHelp();
    static void WriteHelp();
    int processReadArgs(ArgList&);
    int ReadData(FileName const&,DataSetList&,std::string const&);
    int processWriteArgs(ArgList&);
    int WriteData(FileName const&,DataSetList const&);
    bool ID_DataFormat(CpptrajFile&);
    // -------------------------------------------
    int ReadCmatrix(FileName const&, DataSet_PairwiseCache_MEM&);
    int WriteCmatrix(FileName const&, DataSet_PairwiseCache_MEM const&); 
  private:
};
#endif
