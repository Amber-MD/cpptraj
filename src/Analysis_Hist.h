#ifndef INC_ANALYSIS_HIST_H
#define INC_ANALYSIS_HIST_H
// Class: Hist
/// Create an N-dimensional histogram from N input datasets
#include "Analysis.h"
#include "Histogram.h"
#include <vector>
class Hist : public Analysis {
    Histogram hist;
    std::vector<DataSet*> histdata;
    ArgList dimensionArgs;

    bool calcFreeE;
    double Temp;
    bool normalize;
    bool gnuplot;
    bool circular;
    int Ndata;
    char *outfilename;

    double min;
    bool defaultMinSet;
    double max;
    bool defaultMaxSet;
    double step;
    int bins;
    DataSetList histout;

    int CheckDimension(char *,DataSetList *);
    int setupDimension(char *,DataSet*);
  public :
    Hist();
    ~Hist();

    int Setup(DataSetList*);
    int Analyze();
    void Print(DataFileList*);
};
#endif
