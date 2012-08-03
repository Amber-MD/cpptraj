#ifndef INC_DATAIO_GRACE_H
#define INC_DATAIO_GRACE_H
#include "DataIO.h"
// Class: DataIO_Grace
/// Read/write Grace data files.
class DataIO_Grace : public DataIO {
  public:
    DataIO_Grace();
    ~DataIO_Grace();

    int ReadData(DataSetList&);
    int processWriteArgs(ArgList &);
    int WriteData(DataSetList&);
    int WriteDataInverted(DataSetList&);
  private:
    std::string y_label_;
    double ymin_;
    double ystep_;

    char* readbuffer_; 
};
#endif
