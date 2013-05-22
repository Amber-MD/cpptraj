#ifndef INC_ANALYSIS_HIST_H
#define INC_ANALYSIS_HIST_H
#include "Analysis.h"
#include "DataSet_1D.h"
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
    int CheckDimension(std::string const&, DataSetList*);
    int setupDimension(ArgList&, DataSet_1D const&, size_t&);
    int CalcFreeE();
    int Normalize();
    // ---------------------------------
    int BinIndicesToIndex(std::vector<int> const&);
    bool IncrementBinIndices(std::vector<int>&, int, bool&);
    void PrintBins();

    DataFile* outfile_;                  ///< Output DataFile.
    DataSet* hist_;                      ///< Histogram data set.
    std::vector<double> Bins_;           ///< Histogram data - double in case free E calculated
    std::vector<DataSet_1D*> histdata_;  ///< Array of data sets to be binned.
    std::vector<ArgList> dimensionArgs_; ///< Array of args defining histogram dims
    std::vector<Dimension> dimensions_;  ///< Histogram dimensions.

    int debug_;                          ///< Debug level
    bool calcFreeE_;                     ///< If true, calc free E from hist populations.
    double Temp_;                        ///< temperature to calc free E at.
    bool normalize_;                     ///< if true, normalize histogram.
    bool gnuplot_;                       ///< For internal write only
    bool circular_;                      ///< If true, wrap histogram dimensions.
    bool nativeOut_;                     ///< If true, use built in output routine.
    std::string outfilename_;            ///< Stored in case internal write used (DIM > 3)
    size_t N_dimensions_;                ///< # of histogram dimensions.

    Dimension default_dim_;              ///< Hold default dimensions.
    bool minArgSet_;                     ///< True if global min arg has been specified.
    bool maxArgSet_;                     ///< True if global max arg has been specified.
};
#endif
