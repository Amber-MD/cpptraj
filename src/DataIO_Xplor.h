#ifndef INC_DATAIO_XPLOR_H
#define INC_DATAIO_XPLOR_H
#include "DataIO.h"
/// Write Xplor format data files.
class DataIO_Xplor : public DataIO {
  public:
    DataIO_Xplor() : DataIO(false,false,true) {} // Valid for 3D only
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_Xplor(); }
    int processReadArgs(ArgList&) { return 0; }
    int ReadData(FileName const&, DataSetList&, std::string const&);
    int processWriteArgs(ArgList&);
    int WriteData(FileName const&,DataSetList const&);
    bool ID_DataFormat(CpptrajFile&) { return false; }
  private:
    int WriteSet3D(DataSet const&, CpptrajFile&) const;
    void WriteXplorHeader(CpptrajFile&, std::string const&,
                          int, int, int, int, int, int, int, int, int,
                          Matrix_3x3 const&) const;
    std::string title_;
    std::string remark_;
};
#endif
