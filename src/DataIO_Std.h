#ifndef INC_DATAIO_STD_H
#define INC_DATAIO_STD_H
#include "DataIO.h"
// Class: DataIO_Std
/// Read/write standard data files.
class DataIO_Std : public DataIO {
  public:
    DataIO_Std();

    int ReadData(DataSetList&);
    int processWriteArgs(ArgList &);
    int WriteData(DataSetList&);
    int WriteDataInverted(DataSetList&);
    int WriteData2D( DataSet& set );
  private:
    bool writeHeader_;
    bool square2d_;

    // FIXME: Use dimension for this in DataIO
    std::string y_label_;
    double ymin_;
    double ystep_;

    void WriteNameToBuffer(DataSet*, bool);
};
#endif
