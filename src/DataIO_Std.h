#ifndef INC_DATAIO_STD_H
#define INC_DATAIO_STD_H
#include "DataIO.h"
// Class: DataIO_Std
/// Read/write standard data files.
class DataIO_Std : public DataIO {
  public:
    DataIO_Std();

    int ReadData(std::string const&,DataSetList&);
    int processWriteArgs(ArgList &);
    int WriteData(std::string const&,DataSetList&);
    int WriteDataInverted(std::string const&,DataSetList&);
    int WriteData2D(std::string const&, DataSet& set );
  private:
    bool writeHeader_;
    bool square2d_;

    // FIXME: Use dimension for this in DataIO
    std::string y_label_;
    double ymin_;
    double ystep_;

    void WriteNameToBuffer(CpptrajFile&, DataSet*, bool);
};
#endif
