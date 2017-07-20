#ifndef INC_DATAIO_REMLOG_H
#define INC_DATAIO_REMLOG_H
#include <map>
#include "DataIO.h"
#include "BufferedLine.h"
#include "DataSet_RemLog.h"
/// Read Amber replica exchange log data.
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
    typedef DataSet_RemLog::IdxArray IdxArray;

    /// Read remlog header
    int ReadRemlogHeader(BufferedLine&, LogType&, unsigned int) const;
    /// Set up groups and dimensions from replica dim file
    int ReadRemdDimFile(FileName const&, DataSet_RemLog::GdimArray&, ReplicaDimArray&);

    /// Set up replica temperature map
    TmapType SetupTemperatureMap(BufferedLine&, IdxArray&) const;
    /// Set up replica pH map
    TmapType Setup_pH_Map(BufferedLine&, IdxArray&) const;
    /// Count number of Hamiltonian replicas
    int CountHamiltonianReps(BufferedLine&) const;

    /// Open replica logs for all dimensions.
    int OpenMremdDims(std::vector<BufferedLine>&, Sarray const&, LogType);

    Sarray logFilenames_; ///< Replica log file names.
    std::string dimfile_; ///< remd.dim file name
    std::string crdidx_;  ///< Starting coordinate indices
    bool searchForLogs_;  ///< True if need to search for logs for individual dimensions

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
