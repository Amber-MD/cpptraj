#ifndef INC_DATASET_REMLOG_H
#define INC_DATASET_REMLOG_H
#include <vector>
#include "DataSet.h"
#include "ReplicaDimArray.h"
/** Store data from Amber REMD log. For each exchange, the replica index and
  * coordinates index for each replica is stored. These indices start from 1.
  */
class DataSet_RemLog : public DataSet {
  public:
    DataSet_RemLog();
    static DataSet* Alloc() { return (DataSet*)new DataSet_RemLog();}
    /// Replica location in dimension
    enum LocationType { BOTTOM = 0, MIDDLE, TOP };
    // -------------------------------------------
    /// Used to hold replica partner information.
    class GroupReplica;
    typedef std::vector<GroupReplica> GroupArray; ///< Used to hold rep partner info in a group 
    typedef std::vector<GroupArray> GroupDimType; ///< Used to hold group info in a dimension
    typedef std::vector<GroupDimType> GdimArray;  ///< User to hold group info for all dimensions
    GdimArray const& GroupDims() const { return groupDims_; }
    //GdimArray& SetupGroupDims()        { return groupDims_; }
    // -------------------------------------------
    /// Hold topological info for replica in dimension.
    class RepInfo;
    typedef std::vector<RepInfo> RepInfoDimArray;      ///< Hold rep info in each dimension
    typedef std::vector<RepInfoDimArray> RepInfoArray; ///< Hold rep dim info for all replicas
    RepInfoArray const& ReplicaInfo() const { return repInfo_; }
    // -------------------------------------------
    /// Hold info for a single replica at one exchange.
    class ReplicaFrame;
    /// Add given replica frame to specified ensemble member.
    void AddRepFrame(int rep, ReplicaFrame const& frm) { ensemble_[rep].push_back(frm); }
    /// \return replica frame at exchange in specified ensemble member.
    ReplicaFrame const& RepFrame(int exch, int rep) const { return ensemble_[rep][exch];  }
    /// \return replica frame at last exchange in specified ensemble member.
    ReplicaFrame const& LastRepFrame(int rep)       const { return ensemble_[rep].back(); }

    ReplicaDimArray const& DimTypes() const { return repDims_; }
    /// Allocate for given # of replicas, 1 dimension
    void AllocateReplicas(int, ReplicaDimArray const&, int);
    /// Allocate for given # of replicas
    void AllocateReplicas(int, GdimArray const&, ReplicaDimArray const&, int);
    /// \return number of exchanges
    int NumExchange() const;
    /// \return true if ensemble is valid.
    bool ValidEnsemble() const;
    /// Trim last replica frame.
    void TrimLastExchange();
    /// Print data to stdout
    void PrintReplicaStats() const;
    
    // ----- DataSet routines --------------------
    size_t Size()                       const { return ensemble_.size(); }
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&) { return 1; }
#   endif
    void Info()                         const { return;                  }
    void WriteBuffer(CpptrajFile&, SizeArray const&) const { return; }
    int Allocate(SizeArray const&) { return 0; } // TODO implement?
    // TODO: Remove
    void Add( size_t, const void* ) { return;      }
    int Append(DataSet*) { return 1; }
  private:
    /// Hold info for all exchanges of a single replica.
    typedef std::vector<ReplicaFrame> ReplicaArray;
    /// Hold info for all exchanges of all replicas.
    typedef std::vector<ReplicaArray> ReplicaEnsemble;

    /// Setup a single 1D group.
    void SetupDim1Group(int);

    ReplicaEnsemble ensemble_; // [replica][exchange]
    GdimArray groupDims_;      // [dim][group][idx]
    RepInfoArray repInfo_;     // [replica][dim]
    ReplicaDimArray repDims_;  // [dim]
};
// ----- Public Class Definitions ----------------------------------------------
class DataSet_RemLog::ReplicaFrame {
  public:
    ReplicaFrame() :
      temp0_(0.0), PE_x1_(0.0), PE_x2_(0.0),
      replicaIdx_(-1), partnerIdx_(-1), coordsIdx_(-1), dim_(-1), success_(false) {}
    ReplicaFrame(int r, int p, int c, int d, bool s, double t0, double pe1, double pe2) :
      temp0_(t0), PE_x1_(pe1), PE_x2_(pe2),
      replicaIdx_(r), partnerIdx_(p), coordsIdx_(c), dim_(d), success_(s) {} 
    int ReplicaIdx() const { return replicaIdx_; }
    int PartnerIdx() const { return partnerIdx_; }
    int CoordsIdx()  const { return coordsIdx_;  }
    int Dim()        const { return dim_;        }
    bool Success()   const { return success_;    }
    double Temp0()   const { return temp0_;      }
    double PE_X1()   const { return PE_x1_;      }
    double PE_X2()   const { return PE_x2_;      }
  private:
    double temp0_;   ///< Replica bath temperature.
    double PE_x1_;   ///< (HREMD) Potential energy with coords 1.
    double PE_x2_;   ///< (HREMD) Potential energy with coords 2.
    int replicaIdx_; ///< Position in ensemble.
    int partnerIdx_; ///< Position to exchange to.
    int coordsIdx_;  ///< Coordinate index.
    int dim_;        ///< Dimension of the exchange
    bool success_;   ///< Successfully exchanged?
};

class DataSet_RemLog::GroupReplica {
  public:
    GroupReplica() : l_partner_(-1), me_(-1), r_partner_(-1) {}
    GroupReplica(const GroupReplica& rhs) :
      l_partner_(rhs.l_partner_), me_(rhs.me_), r_partner_(rhs.r_partner_) {}
    GroupReplica(int l, int m, int r) : l_partner_(l), me_(m), r_partner_(r) {}
    int L_partner() const { return l_partner_; }
    int Me()        const { return me_;        }
    int R_partner() const { return r_partner_; }
  private:
    int l_partner_, me_, r_partner_;
};

class DataSet_RemLog::RepInfo {
  public:
    RepInfo() : grp_(-1), lp_(-1), rp_(-1), loc_(BOTTOM) {}
    RepInfo(int g, int l, int r, LocationType loc) :
      grp_(g), lp_(l), rp_(r), loc_(loc) {}
    int GroupID() const { return grp_; }
    int LeftID()  const { return lp_;  }
    int RightID() const { return rp_;  }
    LocationType Location() const { return loc_; }
  private:
    int grp_; ///< Group this replica belongs to in this dim
    int lp_;  ///< Left partner of this replica in this dim
    int rp_;  ///< Right partner of this replica in this dim
    LocationType loc_; ///< Location of replica in this dim
};
#endif
