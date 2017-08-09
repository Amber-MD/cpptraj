#ifndef INC_DATAIO_CHARMMOUTPUT_H
#define INC_DATAIO_CHARMMOUTPUT_H
#include "DataIO.h"
/// <Enter description of DataIO_CharmmOutput here>
class DataIO_CharmmOutput : public DataIO {
  public:
    DataIO_CharmmOutput();
    static void ReadHelp();
    static void WriteHelp();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_CharmmOutput(); }
    int processReadArgs(ArgList&);
    int ReadData(FileName const&, DataSetList&, std::string const&);
    int processWriteArgs(ArgList&);
    int WriteData(FileName const&, DataSetList const&);
    bool ID_DataFormat(CpptrajFile&);
};
#endif
