#ifndef INC_ANALYSIS_AVERAGE_H
#define INC_ANALYSIS_AVERAGE_H
#include "Analysis.h"
#include "Array1D.h"
class Analysis_Average : public Analysis {
  public:
    Analysis_Average() : outfile_(0), avgOfSets_(0), sdOfSets_(0),
                         writeHeader_(true), calcAvgOverSets_(false) {}
    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_Average(); }
    static void Help();
  
    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    Array1D input_dsets_;
    CpptrajFile* outfile_;
    DataSet* avgOfSets_;
    DataSet* sdOfSets_;
    bool writeHeader_;
    bool calcAvgOverSets_;
};
#endif
