#ifndef INC_DATAIO_CHARMMRTFPRM_H
#define INC_DATAIO_CHARMMRTFPRM_H
#include "DataIO.h"
#include "BufferedLine.h"
#include "DataSet_Parameters.h"
/// Read in CHARMM topology / parameters. TODO topology read 
class DataIO_CharmmRtfPrm : public DataIO {
  public:
    DataIO_CharmmRtfPrm();
    static void ReadHelp();
    static void WriteHelp();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_CharmmRtfPrm(); }
    int processReadArgs(ArgList&);
    int ReadData(FileName const&, DataSetList&, std::string const&);
    int processWriteArgs(ArgList&);
    int WriteData(FileName const&, DataSetList const&);
    bool ID_DataFormat(CpptrajFile&);
};
#endif
