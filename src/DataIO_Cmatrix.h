#ifndef INC_DATAIO_CMATRIX_H
#define INC_DATAIO_CMATRIX_H
#include "DataIO.h"
#include "DataSet_Cmatrix_MEM.h"
/// Read/write cpptraj-format binary cluster pairwise distance matrix files. 
class DataIO_Cmatrix : public DataIO {
  public:
    DataIO_Cmatrix();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_Cmatrix(); }
    static void ReadHelp();
    static void WriteHelp();
    int processReadArgs(ArgList&);
    int ReadData(FileName const&,DataSetList&,std::string const&);
    int processWriteArgs(ArgList&);
    int WriteData(FileName const&,DataSetList const&);
    bool ID_DataFormat(CpptrajFile&);
    // -------------------------------------------
    int ReadCmatrix(FileName const&, DataSet_Cmatrix_MEM&);
    int WriteCmatrix(FileName const&, DataSet_Cmatrix_MEM const&); 
  private:
    static const unsigned char Magic_[];
    /// For reading/writing 8 byte unsigned integers
    typedef unsigned long long int uint_8;
    /// For reading/writing 8 byte signed integers
    typedef long long int sint_8;
};
#endif
