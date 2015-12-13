#ifndef INC_ANALYSIS_TI_H
#define INC_ANALYSIS_TI_H
#include "Analysis.h"
#include "Array1D.h"
class Analysis_TI : public Analysis {
  public:
    Analysis_TI() : nskip_(0), dAout_(0) {}
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_TI(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    int SetQuadAndWeights(int);

    Array1D input_dsets_; ///< Input DV/DL data sets
    typedef std::vector<int> Iarray;
    Iarray nskip_;        ///< Numbers of data points to skip in calculating <DV/DL>
    DataSet* dAout_;      ///< Free energy data set
    typedef std::vector<DataSet*> DSarray;
    DSarray curve_;       ///< TI curve data set for each skip value
    typedef std::vector<double> Darray;
    Darray quad_;         ///< Hold Gaussian quadrature abscissas
    Darray wgt_;          ///< Hold Gaussian quadrature weights
};
#endif
