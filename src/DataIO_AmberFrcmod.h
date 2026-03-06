#ifndef INC_DATAIO_AMBERFRCMOD_H
#define INC_DATAIO_AMBERFRCMOD_H
#include "DataIO.h"
class DataSet_LeapOpts;
/// <Enter description of DataIO_AmberFrcmod here>
class DataIO_AmberFrcmod : public DataIO {
  public:
    DataIO_AmberFrcmod();
    static void ReadHelp();
    static void WriteHelp();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_AmberFrcmod(); }
    int processReadArgs(ArgList&);
    int ReadData(FileName const&, DataSetList&, std::string const&);
    int processWriteArgs(ArgList&);
    int WriteData(FileName const&, DataSetList const&);
    bool ID_DataFormat(CpptrajFile&);
    // Read in data with leap options
    int ReadData(FileName const&, DataSetList&, std::string const&, DataSet_LeapOpts*);
};
#endif
