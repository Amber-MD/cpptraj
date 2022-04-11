#ifndef INC_DATAIO_EVECS_H
#define INC_DATAIO_EVECS_H
#include "DataIO.h"
/// Read write evecs (eigenmodes) data file.
class DataIO_Evecs : public DataIO {
  public:
    DataIO_Evecs();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_Evecs(); }
    static void ReadHelp();
    int processReadArgs(ArgList &);
    int ReadData(FileName const&,DataSetList&,std::string const&);
    int processWriteArgs(ArgList &)                     { return 0; }
    int WriteData(FileName const&,DataSetList const&);
    bool ID_DataFormat(CpptrajFile&);
  private:
    static const char* MatrixOutputString(MetaData::scalarType);
    int ibeg_, iend_;
};
#endif
