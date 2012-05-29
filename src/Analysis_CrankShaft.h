#ifndef INC_ANALYSIS_CRANKSHAFT_H
#define INC_ANALYSIS_CRANKSHAFT_H
#include "Analysis.h"
class Analysis_CrankShaft : public Analysis {
  public: 
    Analysis_CrankShaft();
    
    int Setup(DataSetList*);
    int Analyze();
  private:
    enum CSangleType { NOTYPE = 0, EPSILON_ZETA, ALPHA_GAMMA };
    enum CStype { ANGLE=0, DISTANCE };
    static const char CSstring[][9];
    static const char distance_ss_2D[][6][9];
    static const char torsion_ss_2D[][6][6];

    std::string filename_;
    int start_;
    int stop_;
    int offset_;
    CStype type_;
    CSangleType angletype_;
    DataSet *scalar1_;
    DataSet *scalar2_;
    std::string info_;
};
#endif
