#ifndef INC_DATAIO_GNUPLOT_H
#define INC_DATAIO_GNUPLOT_H
#include "DataIO.h"
// Class: GnuplotDataFile
/// Read/write Grace data files.
class GnuplotDataFile : public DataIO {
  public:
    GnuplotDataFile();

    int processWriteArgs(ArgList &);
    int WriteData(DataSetList&);
    //int WriteDataInverted(DataSetList&);
  private:
    std::string y_label_;
    double ymin_;
    double ystep_;

    bool useMap_;
    bool printLabels_; 
};
#endif
