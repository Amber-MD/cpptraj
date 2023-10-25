#ifndef INC_DATAIO_AMBERLIB_H
#define INC_DATAIO_AMBERLIB_H
#include "DataIO.h"
/// <Enter description of DataIO_AmberLib here>
class DataIO_AmberLib : public DataIO {
  public:
    DataIO_AmberLib();
    static void ReadHelp();
    static void WriteHelp();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_AmberLib(); }
    int processReadArgs(ArgList&);
    int ReadData(FileName const&, DataSetList&, std::string const&);
    int processWriteArgs(ArgList&);
    int WriteData(FileName const&, DataSetList const&);
    bool ID_DataFormat(CpptrajFile&);
};
#endif
