#ifndef INC_DATAIO_NC_CMATRIX_H
#define INC_DATAIO_NC_CMATRIX_H
#include "DataIO.h"
#include "NC_Cmatrix.h"
/// Read/write NetCDF pairwise matrix files.
class DataIO_NC_Cmatrix : public DataIO {
  public:
    DataIO_NC_Cmatrix() {}
    ~DataIO_NC_Cmatrix() { file_.CloseCmatrix(); }
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_NC_Cmatrix(); }
    static void ReadHelp() {}
    static void WriteHelp() {}
    int processReadArgs(ArgList&) { return 0; }
    int ReadData(FileName const&,DataSetList&,std::string const&);
    int processWriteArgs(ArgList&) { return 0; }
    int WriteData(FileName const&,DataSetList const&);
    bool ID_DataFormat(CpptrajFile& f) { return NC_Cmatrix::ID_Cmatrix(f.Filename()); }
  private:
    NC_Cmatrix file_;
};
#endif
