#ifndef INC_DATAIO_CCP4_H
#define INC_DATAIO_CCP4_H
#include "DataIO.h"
/// Write CCP4 format data files.
class DataIO_CCP4 : public DataIO {
  public:
    DataIO_CCP4() : DataIO(false, false, true) {} // Valid for 3D only
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_CCP4(); }

    bool ID_DataFormat(CpptrajFile&);
    int processReadArgs(ArgList&) { return 0; }
    int ReadData(FileName const&, DataSetList&, std::string const&);
    static void WriteHelp();
    int processWriteArgs(ArgList&);
    int WriteData(FileName const&,DataSetList const&);
  private:
};
#endif
