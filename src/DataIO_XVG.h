#ifndef INC_DATAIO_XVG_H
#define INC_DATAIO_XVG_H
#include "DataIO.h"
/// Gromacs XVG data file
class DataIO_XVG : public DataIO {
  public:
    DataIO_XVG() {}
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_XVG(); }
    int processReadArgs(ArgList&) { return 0; }
    int ReadData(FileName const&,DataSetList&,std::string const&);
    int processWriteArgs(ArgList&) { return 0; }
    int WriteData(FileName const&,DataSetList const&) { return 1; }
    bool ID_DataFormat(CpptrajFile&);
};
#endif
