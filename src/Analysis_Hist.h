#ifndef INC_ANALYSIS_HIST_H
#define INC_ANALYSIS_HIST_H
#include "Analysis.h"
#include "Histogram.h"
// Class: Analysis_Hist
/// Create an N-dimensional histogram from N input datasets
class Analysis_Hist : public Analysis {
  public :
    Analysis_Hist();

    int Setup(DataSetList*);
    int Analyze();
    void Print(DataFileList*);
  private:
    Histogram hist_;
    std::vector<DataSet*> histdata_;
    std::vector<ArgList> dimensionArgs_;

    bool calcFreeE_;
    double Temp_;
    bool normalize_;
    bool gnuplot_;
    bool circular_;
    char *outfilename_;

    Dimension default_dim_;
    bool minArgSet_;
    bool maxArgSet_;

    DataSetList histout_;

    int CheckDimension(char *, DataSetList *);
    int setupDimension(ArgList&, DataSet*);
};
#endif
