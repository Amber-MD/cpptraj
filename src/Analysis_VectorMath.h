#ifndef INC_ANALYSIS_VECTORMATH_H
#define INC_ANALYSIS_VECTORMATH_H
#include "Analysis.h"
// Forward declarations
class DataSet_Vector;
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

    int DotProduct(DataSet*, DataSet_Vector&, DataSet_Vector&,
                   unsigned int, unsigned int, unsigned int) const;
    int CrossProduct(DataSet*, DataSet_Vector&, DataSet_Vector&,
                     unsigned int, unsigned int, unsigned int) const;
    int DoMath(DataSet*, DataSet_Vector&, DataSet_Vector&) const;

    ModeType mode_;
    typedef std::vector<DataSet_Vector*> DVarray;
    DVarray vinfo1_;
    DVarray vinfo2_;
    typedef std::vector<DataSet*> DSarray;
    DSarray DataOut_;                      ///< Output data set(s)
    bool norm_;
};
#endif
