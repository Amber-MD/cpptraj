#ifndef INC_DATAIO_CHARMMFASTREP_H
#define INC_DATAIO_CHARMMFASTREP_H
#include "DataIO.h"
/// Read CHARMM fast temperature replica exchange log data.
class DataIO_CharmmFastRep : public DataIO {
  public:
    DataIO_CharmmFastRep();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_CharmmFastRep(); }
    static void ReadHelp();
    int processReadArgs(ArgList&);
    int ReadData(FileName const&,DataSetList&,std::string const&);
    int processWriteArgs(ArgList&) { return 0; }
    int WriteData(FileName const&, DataSetList const&) { return 1; }
    bool ID_DataFormat(CpptrajFile&);
};
#endif
