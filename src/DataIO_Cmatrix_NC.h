#ifndef INC_DATAIO_CMATRIX_NC_H
#define INC_DATAIO_CMATRIX_NC_H
#include "DataIO.h"
#include "Cluster/Cmatrix_NC.h"
/// Read/write NetCDF pairwise matrix files.
class DataIO_Cmatrix_NC : public DataIO {
  public:
    DataIO_Cmatrix_NC();
    ~DataIO_Cmatrix_NC() { file_.CloseCmatrix(); }
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_Cmatrix_NC(); }
    static void ReadHelp() {}
    static void WriteHelp() {}
    int processReadArgs(ArgList&) { return 0; }
    int ReadData(FileName const&,DataSetList&,std::string const&);
    int processWriteArgs(ArgList&) { return 0; }
    int WriteData(FileName const&,DataSetList const&);
    bool ID_DataFormat(CpptrajFile& f) {
      return Cpptraj::Cluster::Cmatrix_NC::ID_Cmatrix(f.Filename());
    }
  private:
    Cpptraj::Cluster::Cmatrix_NC file_;
};
#endif
