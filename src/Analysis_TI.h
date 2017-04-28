#ifndef INC_ANALYSIS_TI_H
#define INC_ANALYSIS_TI_H
#include "Analysis.h"
#include "Array1D.h"
class Analysis_TI : public Analysis {
    enum ModeType { GAUSSIAN_QUAD = 0, TRAPEZOID };
  public:
    Analysis_TI() : Analysis(HIDDEN), nskip_(0), dAout_(0), mode_(GAUSSIAN_QUAD) {}
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_TI(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    int SetQuadAndWeights(int);

    typedef std::vector<int> Iarray;
    typedef std::vector<double> Darray;
    typedef std::vector<DataSet*> DSarray;

    Array1D input_dsets_; ///< Input DV/DL data sets
    Iarray nskip_;        ///< Numbers of data points to skip in calculating <DV/DL>
    DataSet* dAout_;      ///< Free energy data set
    DSarray curve_;       ///< TI curve data set for each skip value
    Darray xval_;         ///< Hold abscissas corresponding to data sets.
    Darray wgt_;          ///< Hold Gaussian quadrature weights
    ModeType mode_;       ///< Integration mode
    int debug_;
};
#endif
