#ifndef INC_DATAIO_AMBERENE_H
#define INC_DATAIO_AMBERENE_H
#include "DataIO.h"
/// Read Amber ASCII energy file 
class DataIO_AmberEne : public DataIO {
  public:
    DataIO_AmberEne();
    static void ReadHelp();
    static void WriteHelp();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_AmberEne(); }
    int processReadArgs(ArgList&);
    int ReadData(FileName const&, DataSetList&, std::string const&);
    int processWriteArgs(ArgList&);
    int WriteData(FileName const&, DataSetList const&);
    bool ID_DataFormat(CpptrajFile&);
};
#endif
