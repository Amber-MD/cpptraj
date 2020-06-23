#ifndef INC_ANALYSIS_EVALPLATEAU_H
#define INC_ANALYSIS_EVALPLATEAU_H
#include "Analysis.h"
#include "Array1D.h"
class DataSet_Mesh;
/// Can be used to evaluate if a time series has reached a plateau (i.e. near zero slope). 
class Analysis_EvalPlateau : public Analysis {
  public:
    Analysis_EvalPlateau();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_EvalPlateau(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    /// Enumeration for output data
    enum OdataType {
      A0 = 0, ///< The initial density
      A1,     ///< The rate constant
      A2,     ///< The final density
      ONEA1,  ///< 1/A1, the decay time
      //FVAL,   ///< F value from linear regression
      CORR,   ///< Correlation coefficient of exp. fit
      VALA,
      //LCHISQ, ///< Chi^2 of linear fit
      CHISQ,  ///< Hold chi^2 of exp. fit for each input set
      PLTIME, ///< Time at which slope cutoff was satisfied
      FSLOPE, ///< Final slope of the single exponential
      NAME,   ///< Hold name of each input set
      RESULT,
      NDATA
    };
    /// Hold aspect for each output data
    static const char* OdataStr_[];
    /// Hold data type for each output data
    static DataSet::DataType OdataType_[];
    /// Evaluate whether slope cutoff is satisfied.
    bool CheckSlope(DataSet_Mesh const&, double&, double&);
    /// Used to add zero data when error occurs evaluating a set
    void BlankResult(long int, const char*);

    Array1D inputSets_;                ///< Will hold data to evaluate
    std::vector<DataSet*> outputSets_; ///< Will hold final fit curves
    std::string dsname_;               ///< Output set(s) base name
    CpptrajFile* statsout_;            ///< File to write stats to.
    std::vector<DataSet*> data_;       ///< Will hold output data
    double tolerance_;                 ///< Tolerance for non-linear curve fit
    double initpct_;                   ///< Percent of initial data to use as guess for A0
    double finalpct_;                  ///< Percent of final data to use as guess for A2
    double valaCut_;                   ///< Cutoff for long-term estimate from last half of data
    double chisqCut_;                  ///< Cutoff for non-linear fit chi^2
    double slopeCut_;                  ///< Cutoff for non-linear fit slope
    int maxIt_;                        ///< Max iterations to perform non-linear fit
    int debug_;
};
#endif
