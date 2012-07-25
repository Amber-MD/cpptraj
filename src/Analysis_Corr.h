#ifndef INC_ANALYSIS_CORR_H
#define INC_ANALYSIS_CORR_H
#include "Analysis.h"
#include "DataSet_double.h"
// Class: Analysis_Corr
/// Calculate autocorrelation or correlation of dataset(s).
class Analysis_Corr : public Analysis {
  public:
    Analysis_Corr();

    int Setup(DataSetList*);
    int Analyze();
    void Print(DataFileList*);
  private:
    DataSet *D1_;
    DataSet *D2_;
    int lagmax_;
    int Nelements_;
    DataSet_double Ct_;
    char *outfilename_;
};
#endif
