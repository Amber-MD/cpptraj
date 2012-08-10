#ifndef INC_ANALYSIS_AUTOCORR_H
#define INC_ANALYSIS_AUTOCORR_H
#include "Analysis.h"
class Analysis_AutoCorr : public Analysis {
  public:
    Analysis_AutoCorr();

    int Setup( DataSetList* );
    int Analyze();
    void Print( DataFileList* );

  private:
    DataSetList dsets_;
    std::string outfilename_;
    std::string setname_;
    std::vector<DataSet*> outputData_;
    int lagmax_;
};
#endif
