#ifndef INC_DATAIO_XPLOR_H
#define INC_DATAIO_XPLOR_H
#include "DataIO.h"
/// Write Xplor format data files.
class DataIO_Xplor : public DataIO {
  public:
    DataIO_Xplor() {}
    static DataIO* Alloc() { return (DataIO*)new DataIO_Xplor(); }
    int ReadData(std::string const&, ArgList&, DataSetList&, std::string const&);
    int processWriteArgs(ArgList&)                 { return 0; }
    int WriteData(std::string const&,DataSetList const&)         { return 1; }
    int WriteDataInverted(std::string const&,DataSetList const&) { return 1; }
    int WriteData2D(std::string const&, DataSet const&)         { return 1; }
    int WriteData3D(std::string const&, DataSet const&);
    bool ID_DataFormat(CpptrajFile&) { return false; }
  private:
    std::string title_;
    std::string remark_;
};
#endif
