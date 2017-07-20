#ifndef INC_DATAIO_CHARMMREPLOG_H
#define INC_DATAIO_CHARMMREPLOG_H
#include "DataIO.h"
/// Read CHARMM replica exchange log data.
class DataIO_CharmmRepLog : public DataIO {
  public:
    DataIO_CharmmRepLog();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_CharmmRepLog(); }
    static void ReadHelp();
    int processReadArgs(ArgList&);
    int ReadData(FileName const&,DataSetList&,std::string const&);
    int processWriteArgs(ArgList&) { return 0; }
    int WriteData(FileName const&, DataSetList const&) { return 1; }
    bool ID_DataFormat(CpptrajFile&);
  private:
    int ReadReplogArray(FileName const&, DataSetList&, std::string const&);

    std::string crdidx_; ///< Hold user-specified starting coord indices arg
    int nrep_;           ///< Number of replicas.
};
#endif
