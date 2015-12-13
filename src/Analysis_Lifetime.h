#ifndef INC_ANALYSIS_LIFETIME_H
#define INC_ANALYSIS_LIFETIME_H
#include "Analysis.h"
#include "Array1D.h"
class Analysis_Lifetime : public Analysis {
  public:
    Analysis_Lifetime();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_Lifetime(); }
    void Help() const;
    Analysis::RetType ExternalSetup(Array1D const&, CpptrajFile*);
    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    Array1D inputDsets_;
    Array1D outputDsets_;
    Array1D curveSets_;
    Array1D maxDsets_;
    Array1D avgDsets_;
    CpptrajFile* standalone_;
    int windowSize_;
    int fuzzCut_;
    double cut_;
    bool averageonly_;
    bool cumulative_;
    bool deltaAvg_;
    bool normalizeCurves_; ///< If true normalize lifetime curves
    typedef bool (*CompareFxn)(double,double);
    CompareFxn Compare_;
    static bool Compare_GreaterThan(double l, double r) { return (l > r); }
    static bool Compare_LessThan(double l, double r)    { return (l < r); }
};
#endif
