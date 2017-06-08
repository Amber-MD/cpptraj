#ifndef INC_DATAIO_CHARMMREPLOG_H
#define INC_DATAIO_CHARMMREPLOG_H
#include <map>
#include "DataIO.h"
#include "BufferedLine.h"
#include "DataSet_RemLog.h"
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
    // NOTE: Must match LogDescription TODO Combine with DataIO_RemLog?
    enum LogType { UNKNOWN = 0, TREMD, HREMD, MREMD, RXSGLD, PHREMD };
    static const char* LogDescription[];
};
#endif
