#ifndef INC_ANALYSIS_MATRIX_H
#define INC_ANALYSIS_MATRIX_H
#include "Analysis.h"
#include "MatrixType.h"
class Analysis_Matrix : public Analysis {
  public:
    Analysis_Matrix();

    int Setup(DataSetList*);
  private:
    enum ThermoFlag { OFF=0, THERMO, ORDER };

    MatrixType* minfo_;
    std::string outfilename_; // carg2
    std::string orderparamfile_; // carg4
    ThermoFlag thermo_; // iarg1
    int nevec_; // iarg2
    bool reduce_; // iarg3
    ModesInfo* modinfo_; // carg3
};
#endif
