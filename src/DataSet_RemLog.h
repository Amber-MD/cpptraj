#ifndef INC_DATASET_REMLOG_H
#define INC_DATASET_REMLOG_H
#include <map>
#include "DataSet.h"
class DataSet_RemLog : public DataSet {
  public:
    DataSet_RemLog();
    /// T-REMD temperature map.
    typedef std::map<double,int> TmapType;
    /// Hold info for a single replica at one exchange.
    class ReplicaFrame;
    /// Hold info for all exchanges of a single replica.
    typedef std::vector<ReplicaFrame> ReplicaArray;
    /// Allocate for given # of replicas
    void AllocateReplicas(int);
    /// \return replica in ensemble.
    ReplicaArray& Replica(int i)  { return ensemble_[i];         }
    /// Ensemble iterator
    typedef std::vector<ReplicaArray>::const_iterator ensemble_it;
    ensemble_it begin()     const { return ensemble_.begin();    }
    ensemble_it end()       const { return ensemble_.end();      }
    /// Iterator to replica in ensemble
    typedef ReplicaArray::const_iterator replica_it;
    replica_it begin(int i) const { return ensemble_[i].begin(); }
    replica_it end(int i)   const { return ensemble_[i].end();   }
    // ----- DataSet routines --------------------
    int Size() { return ensemble_.size(); }
  private:
    /// Hold info for all exchanges of all replicas.
    typedef std::vector<ReplicaArray> ReplicaEnsemble;
    ReplicaEnsemble ensemble_;
};
// ----- Public Class Definitions ----------------------------------------------
class DataSet_RemLog::ReplicaFrame {
  public:
    ReplicaFrame() : replicaIdx_(-1), partnerIdx_(-1), coordsIdx_(-1),
                     success_(false), temp0_(0.0),
                     PE_x1_(0.0), PE_x2_(0.0) {}
    int SetTremdFrame( const char*, TmapType const& );
    int SetHremdFrame( const char*, int );
    int ReplicaIdx() const { return replicaIdx_; }
    int PartnerIdx() const { return partnerIdx_; }
    int CoordsIdx()  const { return coordsIdx_;  }
    bool Success()   const { return success_;    }
    double Temp0()   const { return temp0_;      }
    double PE_X1()   const { return PE_x1_;      }
    double PE_X2()   const { return PE_x2_;      }
  private:
    int replicaIdx_; ///< Position in ensemble.
    int partnerIdx_; ///< Position to exchange to.
    int coordsIdx_;  ///< Coordinate index.
    bool success_;   ///< Successfully exchanged.
    double temp0_;   ///< Replica bath temperature.
    double PE_x1_;   ///< (HREMD) Potential energy with coords 1.
    double PE_x2_;   ///< (HREMD) Potential energy with coords 2.
};
#endif
