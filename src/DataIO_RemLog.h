#ifndef INC_DATAIO_REMLOG_H
#define INC_DATAIO_REMLOG_H
#include <map>
#include "DataIO.h"
#include "BufferedLine.h"
#include "DataSet_RemLog.h"
/// Read replica exchange log data.
class DataIO_RemLog : public DataIO {
  public:
    DataIO_RemLog();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_RemLog(); }
    static void ReadHelp();
    int processReadArgs(ArgList&);
    int ReadData(FileName const&,DataSetList&,std::string const&);
    int processWriteArgs(ArgList&) { return 0; }
    int WriteData(FileName const&, DataSetList const&) { return 1; }
    bool ID_DataFormat(CpptrajFile&) { return false; }
  private:
    // NOTE: Must match LogDescription
    enum LogType { UNKNOWN = 0, TREMD, HREMD, MREMD, RXSGLD, PHREMD };
    static const char* LogDescription[];
    typedef std::vector<std::string> Sarray; // TODO FileName array?
    typedef std::map<double,int> TmapType; // FIXME: Use ReplicaMap

    int ReadRemlogHeader(BufferedLine&, LogType&, unsigned int) const;
    int ReadRemdDimFile(FileName const&, DataSet_RemLog::GdimArray&);
    TmapType SetupTemperatureMap(BufferedLine&,std::vector<int>&) const;
    TmapType Setup_pH_Map(BufferedLine&, std::vector<int>&) const;
    int CountHamiltonianReps(BufferedLine&) const;
    int OpenMremdDims(std::vector<BufferedLine>&, Sarray const&, unsigned int);
    void SetupDim1Group( int, DataSet_RemLog::GdimArray& );
    void PrintReplicaStats(DataSet_RemLog const&);

    Sarray logFilenames_; ///< Replica log file names.
    std::string dimfile_;
    std::string crdidx_;
    int n_mremd_replicas_;
    bool processMREMD_;
    bool searchForLogs_;

    //std::vector<ExchgType> DimTypes_;
    ReplicaDimArray DimTypes_;
    // Used for getting temps/coord indices from T-remlog
    struct TlogType {
      double t0;
      int crdidx;
    };
    struct TlogType_cmp {
      inline bool operator()(TlogType const& first, TlogType const& second) const {
        return (first.t0 < second.t0);
      }
    };
};
#endif
