#ifndef INC_ANALYSIS_HIST_H
#define INC_ANALYSIS_HIST_H
#include "Analysis.h"
#include "Histogram.h"
// Class: Analysis_Hist
/// Create an N-dimensional histogram from N input datasets
class Analysis_Hist : public Analysis {
  public :
    Analysis_Hist();

    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_Hist(); }
    static void Help();

    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    DataFile* outfile_;                  ///< Output DataFile for 1 or 2 dims.
    Histogram* hist_;
    std::vector<DataSet*> histdata_;
    std::vector<ArgList> dimensionArgs_;

    int debug_;
    bool calcFreeE_;
    double Temp_;
    bool normalize_;
    bool gnuplot_;
    bool circular_;
    std::string outfilename_;

    Dimension default_dim_;
    bool minArgSet_;
    bool maxArgSet_;

    int CheckDimension(std::string const&, DataSetList*);
    int setupDimension(ArgList&, DataSet*);
};
#endif
