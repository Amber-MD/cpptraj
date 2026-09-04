#ifndef INC_DATAIO_GBNSR6_H
#define INC_DATAIO_GBNSR6_H
#include "DataIO.h"
/// Read energies from Amber GBNSR6 output files. 
class DataIO_GBNSR6 : public DataIO {
  public:
    DataIO_GBNSR6();
    static void ReadHelp();
    static void WriteHelp();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_GBNSR6(); }
    int processReadArgs(ArgList&);
    int ReadData(FileName const&, DataSetList&, std::string const&);
    int processWriteArgs(ArgList&);
    int WriteData(FileName const&, DataSetList const&);
    bool ID_DataFormat(CpptrajFile&);
};
#endif
