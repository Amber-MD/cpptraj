#ifndef INC_DATAIO_OPENDX_H
#define INC_DATAIO_OPENDX_H
#include "DataIO.h"
/// Write OpenDx format data files.
class DataIO_OpenDx : public DataIO {
  public:
    DataIO_OpenDx() : DataIO(false, false, true), gridWriteMode_(BIN_CORNER) {} // Valid for 3D only
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_OpenDx(); }
    int processReadArgs(ArgList&) { return 0; }
    int ReadData(FileName const&, DataSetList&, std::string const&);
    static void WriteHelp();
    int processWriteArgs(ArgList&);
    int WriteData(FileName const&,DataSetList const&);
    bool ID_DataFormat(CpptrajFile&);
  private:
    enum GridWriteType { BIN_CORNER = 0, BIN_CENTER, WRAP, EXTENDED };
    int LoadGrid(const char*, DataSet&);
    int WriteSet3D(DataSet const&, CpptrajFile&) const;
    void WriteDxHeader(CpptrajFile&, size_t, size_t, size_t, double, double, double,
                       Matrix_3x3 const&, Vec3 const&) const;
    int WriteGridWrap(DataSet const&, CpptrajFile&) const;
    int WriteGrid(DataSet const&, CpptrajFile&) const;

    GridWriteType gridWriteMode_;
};
#endif
