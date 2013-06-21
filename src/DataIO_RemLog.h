#ifndef INC_DATAIO_REMLOG_H
#define INC_DATAIO_REMLOG_H
#include "DataIO.h"
// Class: DataIO_RemLog
/// Read replica exchange log data.
class DataIO_RemLog : public DataIO {
  public:
    DataIO_RemLog();

    int ReadData(std::string const&,DataSetList&);
  private:
    class ReplicaFrame {
      public:
        ReplicaFrame() : idx_(-1), temp0_(0.0), PE_x1_(0.0), PE_x2_(0.0) {}
      private:
        int idx_;      ///< Position in ensemble
        double temp0_; ///< Replica bath temperature.
        double PE_x1_; ///< (HREMD) Potential energy with coords 1.
        double PE_x2_; ///< (HREMD) Potential energy with coords 2.
    };
    class ReplicaArray {
      public:
        ReplicaArray() : originalIdx_(-1), currentIdx_(-1) {}
      private:
        int originalIdx_;
        int currentIdx_;
        std::vector<ReplicaFrame> frames_;
    };
    typedef std::vector<ReplicaArray> ReplicaEnsemble;
};
#endif
