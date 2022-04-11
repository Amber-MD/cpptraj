#ifndef INC_ANALYSIS_REMLOG_H
#define INC_ANALYSIS_REMLOG_H
#include "Analysis.h"
#include "DataSet_RemLog.h"
#include "DataSet_integer_mem.h"
class Analysis_RemLog : public Analysis {
  public:
    Analysis_RemLog();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_RemLog(); }
    void Help() const;
    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    enum ModeType { NONE = 0, CRDIDX, REPIDX };

    typedef std::vector<int> Iarray;
    typedef std::vector<DataSet_integer_mem> DSI_array;
    /// Track exchange stats for each dimension.
    class RepStats {
      public:
        RepStats(int nreps) :
          acceptUp_(nreps, 0), acceptDown_(nreps, 0), attempts_(0) {}

        Iarray acceptUp_;     ///< # exchanges to the right accepted for each replica.
        Iarray acceptDown_;   ///< # exchanges to the left accepted for each replica.
        int attempts_;        ///< Total # exchanges in dimension.
    };
    /// Trip status
    enum RepStatusType { UNKNOWN = 0, HIT_BOTTOM, HIT_TOP };
    /// Track trip stats for each dimension
    class TripStats {
      public:
        TripStats(int nreps) :
          status_(nreps, UNKNOWN), bottom_(nreps, 0), roundTrip_(nreps) {}

        Iarray status_;       ///< Current status of each crdidx.
        Iarray bottom_;       ///< Frame at which each crdidx hit bottom.
        DSI_array roundTrip_; ///< Array of round trip times for each crdidx.
    };

    int debug_;
    bool calculateStats_;
    bool calculateLifetimes_;
    bool printIndividualTrips_;
    DataSet_RemLog* remlog_;
    ModeType mode_;
    std::vector<DataSet*> outputDsets_;
    std::vector<DataSet*> eSets_;  ///< Hold energies extracted from replica logs
    DataFile* lifetimes_;
    CpptrajFile* statsout_;
    CpptrajFile* reptime_;
    CpptrajFile* acceptout_;
    int calcRepFracSlope_;
    CpptrajFile* repFracSlope_;
    std::string dsname_; ///< Output data set name
    AnalysisSetup Setup_; ///< Hold DSL for lifetime data sets
};
#endif
