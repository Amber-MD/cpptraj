#ifndef INC_ANALYSIS_CRANKSHAFT_H
#define INC_ANALYSIS_CRANKSHAFT_H
#include "Analysis.h"
class Analysis_CrankShaft : public Analysis {
  public: 
    Analysis_CrankShaft();
    
    int Setup(DataSetList*);
    int Analyze();
  private:
    std::string filename_;
    int start_;
    int stop_;
    int offset_;
    enum CStype { ANGLE=0, DISTANCE };
    CStype type_;
    DataSet *scalar1_;
    DataSet *scalar2_;
    std::string info_;
};
#endif
