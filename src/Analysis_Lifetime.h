#ifndef INC_ANALYSIS_LIFETIME_H
#define INC_ANALYSIS_LIFETIME_H
#include "Analysis.h"
class Analysis_Lifetime : public Analysis {
  public:
    Analysis_Lifetime();

    int Setup(DataSetList*);
    int Analyze();
    void Print(DataFileList*);
  private:
    DataSetList inputDsets_;
    std::string outfilename_;
    std::vector<DataSet*> outputDsets_;
    std::vector<DataSet*> maxDsets_;
    std::vector<DataSet*> avgDsets_;
    int windowSize_;
    double cut_;
    bool averageonly_;
    bool cumulative_;

    void PrintListToFile(DataFileList *, std::vector<DataSet*>&, std::string const&);
};
#endif
