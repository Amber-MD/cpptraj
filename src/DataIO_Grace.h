#ifndef INC_DATAIO_GRACE_H
#define INC_DATAIO_GRACE_H
#include "DataIO.h"
// Class: GraceDataFile
/// Read/write Grace data files.
class GraceDataFile : public DataIO {
  public:
    GraceDataFile();

    int processWriteArgs(ArgList &);
    int WriteData(DataSetList&);
    int WriteDataInverted(DataSetList&);
  private:
    std::string y_label_;
    double ymin_;
    double ystep_; 
};
#endif
