#ifndef INC_ANALYSIS_EXPCURVEFIT_H
#define INC_ANALYSIS_EXPCURVEFIT_H
#include "Analysis.h"
/// Perform curve fit to sum of exponentials.
class Analysis_ExpCurveFit : public Analysis {
  public:
    Analysis_ExpCurveFit() : dset_(0), finalY_(0), tolerance_(0.0), maxIt_(0),
                             nexp_(0) {}
    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_ExpCurveFit(); }
    static void Help();

    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    DataSet* dset_;
    DataSet* finalY_;
    double tolerance_;
    int maxIt_;
    int nexp_; ///< Number of exponentials
};
#endif
