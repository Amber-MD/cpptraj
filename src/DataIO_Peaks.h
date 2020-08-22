#ifndef INC_DATAIO_PEAKS_H
#define INC_DATAIO_PEAKS_H
#include "DataIO.h"
/// <Enter description of DataIO_Peaks here>
class DataIO_Peaks : public DataIO {
  public:
    DataIO_Peaks();
    static void ReadHelp();
    static void WriteHelp();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_Peaks(); }
    int processReadArgs(ArgList&);
    int ReadData(FileName const&, DataSetList&, std::string const&);
    int processWriteArgs(ArgList&);
    int WriteData(FileName const&, DataSetList const&);
    bool ID_DataFormat(CpptrajFile&);
};
#endif
