#ifndef INC_ANALYSIS_AVERAGE_H
#define INC_ANALYSIS_AVERAGE_H
#include "Analysis.h"
#include "Array1D.h"
class Analysis_Average : public Analysis {
  public:
    Analysis_Average();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_Average(); }
    void Help() const;
  
    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    Array1D input_dsets_;   ///< Input sets to average.
    DataSet* avgOfSets_;    ///< Series holding average over all input sets at each point.
    DataSet* sdOfSets_;     ///< Series holding stdev over all input sets at each point.
    DataSet* data_avg_;     ///< Overall average of each set.
    DataSet* data_sd_;      ///< Overall stdev of each set.
    DataSet* data_ymin_;    ///< Minimum Y value of each set.
    DataSet* data_ymax_;    ///< Maximum Y value of each set.
    DataSet* data_yminIdx_; ///< Index at which minimum Y value found in each set.
    DataSet* data_ymaxIdx_; ///< Index at which maximum Y value found in each set.
    DataSet* data_names_;   ///< Legend of each input data set.
    bool calcAvgOverSets_;  ///< If true calculate avg/stdev over all input sets at each point.
    bool toStdout_;         ///< If true write avg results to STDOUT when no outfile specified
};
#endif
