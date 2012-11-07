#ifndef INC_DATAIO_GRACE_H
#define INC_DATAIO_GRACE_H
#include "DataIO.h"
// Class: DataIO_Grace
/// Read/write Grace data files.
class DataIO_Grace : public DataIO {
  public:
    DataIO_Grace();

    int ReadData(std::string const&,DataSetList&);
    int processWriteArgs(ArgList &);
    int WriteData(std::string const&,DataSetList&);
    int WriteDataInverted(std::string const&,DataSetList&);
  private:
    std::string y_label_;
    double ymin_;
    double ystep_;
};
#endif
