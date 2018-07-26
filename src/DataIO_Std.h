#ifndef INC_DATAIO_STD_H
#define INC_DATAIO_STD_H
#include "DataIO.h"
/// Read/write standard data files.
class DataIO_Std : public DataIO {
  public:
    DataIO_Std();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_Std(); }
    static void ReadHelp();
    static void WriteHelp();
    int processReadArgs(ArgList&);
    int ReadData(FileName const&,DataSetList&,std::string const&);
    int processWriteArgs(ArgList&);
    int WriteData(FileName const&,DataSetList const&);
    bool ID_DataFormat(CpptrajFile&) { return false; }
  private:
    static const char* SEPARATORS;
    static int Get3Double(std::string const&, Vec3&);
    int Read_1D(std::string const&,DataSetList&,std::string const&);
    int Read_2D(std::string const&,DataSetList&,std::string const&);
    int Read_3D(std::string const&,DataSetList&,std::string const&);
    int Read_Vector(std::string const&,DataSetList&,std::string const&);
    int Read_Mat3x3(std::string const&,DataSetList&,std::string const&);
    static void WriteNameToBuffer(CpptrajFile&, std::string const&, int,  bool);
    int WriteCmatrix(CpptrajFile&, DataSetList const&);
    int WriteDataNormal(CpptrajFile&,DataSetList const&);
    int WriteDataInverted(CpptrajFile&,DataSetList const&);
    int WriteData2D(CpptrajFile&, DataSetList const&);
    int WriteData3D(CpptrajFile&, DataSetList const&);
    int WriteSet2D(DataSet const&, CpptrajFile&);
    int WriteSet3D(DataSet const&, CpptrajFile&);
    enum modeType {READ1D=0, READ2D, READ3D, READVEC, READMAT3X3};
    enum precType {UNSPEC, FLOAT, DOUBLE};
    modeType mode_;    ///< Read mode
    precType prec_;
    int indexcol_;     ///< Read: column containing index (X) values
    bool isInverted_;  ///< For 1D writes invert X/Y.
    bool hasXcolumn_;
    bool writeHeader_;
    bool square2d_;
    Vec3 origin_; ///< 3d reads, grid origin
    Vec3 delta_; ///< 3d reads, grid delta
    size_t dims_[3]; ///< 3d reads, grid dims
};
#endif
