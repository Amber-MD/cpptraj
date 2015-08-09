#ifndef INC_DATAIO_OPENDX_H
#define INC_DATAIO_OPENDX_H
#include "DataIO.h"
/// Write OpenDx format data files.
class DataIO_OpenDx : public DataIO {
  public:
    DataIO_OpenDx() : DataIO(false, false, true) {} // Valid for 3D only
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_OpenDx(); }
    int processReadArgs(ArgList&) { return 0; }
    int ReadData(FileName const&, DataSetList&, std::string const&);
    int processWriteArgs(ArgList&)                               { return 0; }
    int WriteData(FileName const&,DataSetList const&)         { return 1; }
    bool ID_DataFormat(CpptrajFile&);
  private:
    int LoadGrid(const char*, DataSet&);
    int WriteSet3D( DataSet const&, CpptrajFile&);
};
#endif
