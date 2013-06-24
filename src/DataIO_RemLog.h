#ifndef INC_DATAIO_REMLOG_H
#define INC_DATAIO_REMLOG_H
#include "DataIO.h"
// Class: DataIO_RemLog
/// Read replica exchange log data.
class DataIO_RemLog : public DataIO {
  public:
    DataIO_RemLog();

    int ReadData(std::string const&,DataSetList&);
};
#endif
