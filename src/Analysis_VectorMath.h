#ifndef INC_ANALYSIS_VECTORMATH_H
#define INC_ANALYSIS_VECTORMATH_H
#include "Analysis.h"
#include "DataSet_Vector.h"
class Analysis_VectorMath : public Analysis {
  public:
    Analysis_VectorMath();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_VectorMath(); }
    void Help() const;
    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    enum ModeType {  DOTPRODUCT = 0, DOTANGLE, CROSSPRODUCT };
    static const char* ModeString[];

    int DotProduct();
    int CrossProduct();

    ModeType mode_;
    DataSet_Vector* vinfo1_;
    DataSet_Vector* vinfo2_;
    DataSet* DataOut_;       ///< Output data set
    bool norm_;
};
#endif
