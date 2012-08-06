#ifndef INC_ANALYSIS_CROSSCORR_H
#define INC_ANALYSIS_CROSSCORR_H
#include "Analysis.h"
class Analysis_CrossCorr : public Analysis {
  public:
    Analysis_CrossCorr();

    int Setup( DataSetList* );
    int Analyze();
    void Print( DataFileList* );

  private:
    DataSetList dsets_;
    std::string outfilename_;
    std::string setname_;
    DataSet* matrix_;
    //std::string Xlabels_;
    std::string Ylabels_;
};
#endif
