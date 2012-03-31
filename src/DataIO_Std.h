#ifndef INC_DATAIO_STD_H
#define INC_DATAIO_STD_H
#include "DataIO.h"
// Class: StdDataFile
/// Read/write standard data files.
class StdDataFile : public DataIO {
  public:
    StdDataFile();

    int processWriteArgs(ArgList &);
    int WriteData(DataSetList&);
    int WriteDataInverted(DataSetList&);
  private:
    bool writeHeader_;
};
#endif
