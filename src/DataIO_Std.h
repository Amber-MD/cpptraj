#ifndef INC_DATAIO_STD_H
#define INC_DATAIO_STD_H
#include "DataIO.h"
// Class: DataIO_Std
/// Read/write standard data files.
class DataIO_Std : public DataIO {
  public:
    DataIO_Std();
    static DataIO* Alloc() { return (DataIO*)new DataIO_Std(); }
    static void ReadHelp();
    static void WriteHelp();
    int ReadData(std::string const&,ArgList&,DataSetList&,std::string const&);
    int processWriteArgs(ArgList&);
    int WriteData(std::string const&,DataSetList const&);
    int WriteDataInverted(std::string const&,DataSetList const&);
    int WriteData2D(std::string const&, DataSet const&);
    int WriteData3D(std::string const&, DataSet const&);
    bool ID_DataFormat(CpptrajFile&) { return false; }
  private:
    static void WriteNameToBuffer(CpptrajFile&, std::string const&, int,  bool);
    bool hasXcolumn_;
    bool writeHeader_;
    bool square2d_;
};
#endif
