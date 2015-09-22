#ifndef INC_ANALYSIS_REMLOG_H
#define INC_ANALYSIS_REMLOG_H
#include "Analysis.h"
#include "DataSet_RemLog.h"
class Analysis_RemLog : public Analysis {
  public:
    Analysis_RemLog();
    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_RemLog(); }
    static void Help();
    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    enum ModeType { NONE = 0, CRDIDX, REPIDX };

    /// Track stats for each dimension.
    typedef std::vector<int> Iarray;
    class RepStats {
      public:
        RepStats(int nreps) : acceptUp_(nreps, 0), acceptDown_(nreps, 0), attempts_(0) {}

        Iarray acceptUp_;
        Iarray acceptDown_;
        int attempts_;
    };

    int debug_;
    bool calculateStats_;
    bool calculateLifetimes_;
    bool printIndividualTrips_;
    DataSet_RemLog* remlog_;
    ModeType mode_;
    std::vector<DataSet*> outputDsets_;
    CpptrajFile* lifetimes_;
    CpptrajFile* statsout_;
    CpptrajFile* reptime_;
    CpptrajFile* acceptout_;
    int calcRepFracSlope_;
    CpptrajFile* repFracSlope_;
};
#endif
