#ifndef INC_ANALYSIS_MODES_H
#define INC_ANALYSIS_MODES_H
#include "Analysis.h"
#include "ModesInfo.h"
class Analysis_Modes : public Analysis {
  public:
    Analysis_Modes();
    ~Analysis_Modes();

    int Setup(DataSetList*);
    int Analyze();

  private:
    enum modeAnalysisType { FLUCT=0, DISPLACE, CORR };
    static const char analysisTypeString[][22];
    modeAnalysisType type_; // iarg1
    int beg_;
    int end_;
    bool bose_;
    double factor_;
    ModesInfo* modinfo_;
    ModesInfo::modesSource source_;
    std::string filename_;
    ModesInfo::modestackType atompairStack_;
    double* results_;
};
#endif
