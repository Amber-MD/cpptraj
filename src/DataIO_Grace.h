#ifndef INC_DATAIO_GRACE_H
#define INC_DATAIO_GRACE_H
#include "DataIO.h"
// Class: DataIO_Grace
/// Read/write Grace data files.
class DataIO_Grace : public DataIO {
  public:
    DataIO_Grace() {}
    static DataIO* Alloc() { return (DataIO*)new DataIO_Grace(); } 
    int ReadData(std::string const&,ArgList&,DataSetList&,std::string const&);

    int processWriteArgs(ArgList &);
    int WriteData(std::string const&,DataSetList const&);
    int WriteDataInverted(std::string const&,DataSetList const&);
    int WriteData2D(std::string const&, DataSet const&) { return 1; }
    int WriteData3D(std::string const&, DataSet const&) { return 1; }
    bool ID_DataFormat(CpptrajFile&) { return false; }
};
#endif
