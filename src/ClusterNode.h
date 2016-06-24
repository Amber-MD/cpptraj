#ifndef INC_CLUSTERNODE_H
#define INC_CLUSTERNODE_H
#include "ClusterDist.h"
#include "DataSet_Cmatrix.h"
/// Hold information for a cluster in a ClusterList
class ClusterNode {
  public:
    ClusterNode();
    ~ClusterNode();
    ClusterNode(ClusterDist*,ClusterDist::Cframes const&, int);
    ClusterNode(const ClusterNode&);
    ClusterNode& operator=(const ClusterNode&);
    /// Used to sort clusters by # of frames in cluster
    inline bool operator<(const ClusterNode&) const;
    /// Merge frames from another cluster to this cluster
    inline void MergeFrames(ClusterNode const&);
    /// Find and set frame in the cluster that has lowest distance to all other frames.
    int FindBestRepFrame(ClusterDist* Cdist);
    /// Calculate eccentricity for frames in this cluster.
    void CalcEccentricity(DataSet_Cmatrix const&);
    /// Calculate centroid of members of this cluster.
    void CalculateCentroid(ClusterDist* Cdist) {
      // FIXME: Could potentially get rid of this branch.
      if (centroid_ == 0)
        centroid_ = Cdist->NewCentroid( frameList_ );
      else
        Cdist->CalculateCentroid( centroid_, frameList_ );
    }
    /// Calculate average distance of all members to centroid
    double CalcAvgToCentroid( ClusterDist*);
    // Iterator over frame numbers
    typedef ClusterDist::Cframes::const_iterator frame_iterator;
    frame_iterator beginframe() const { return frameList_.begin(); }
    frame_iterator endframe()   const { return frameList_.end();   }
    /// \return Frame number at given index.
    int ClusterFrame(int idx)   const { return frameList_[idx];    }
    // Return internal variables
    double Eccentricity()      const { return eccentricity_;          }
    int Num()                  const { return num_;                   }
    int Nframes()              const { return (int)frameList_.size(); }
    int BestRepFrame()         const { return repFrame_;              }
    Centroid* Cent()           const { return centroid_;              }
    std::string const& Cname() const { return name_;                  }
    double RefRms()            const { return refRms_;                }
    // Set internal variables 
    void AddFrameToCluster(int fnum)   { frameList_.push_back( fnum );  }
    void SetNum(int numIn)             { num_ = numIn;                  }
    inline void SetName(std::string const&, double);
    /// Sort internal frame list
    void SortFrameList();
    /// \return true if given frame is in this cluster.
    bool HasFrame(int) const;
    /// Remove specified frame from cluster if present.
    void RemoveFrameFromCluster(int);
    /// Remove specified frame from cluster and update centroid.
    void RemoveFrameUpdateCentroid(ClusterDist*, int);
    /// Add specified fram to cluster and update centroid.
    void AddFrameUpdateCentroid(ClusterDist*, int);
  private:
    double eccentricity_;             ///< Maximum distance between any 2 frames.
    double refRms_;                   ///< Cluster rms to reference (if assigned)
    int num_;                         ///< Cluster number.
    int repFrame_;                    ///< Frame number with lowest dist. to all other frames.
    ClusterDist::Cframes frameList_;  ///< List of frames belonging to this cluster.
    Centroid* centroid_;              ///< Centroid of all frames in this cluster.
    std::string name_;                ///< Cluster name assigned from reference.
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
/** Use > since we give higher priority to larger clusters. */
bool ClusterNode::operator<(const ClusterNode& rhs) const {
  return ( frameList_.size() > rhs.frameList_.size() );
}
/** Frames from rhs go to this cluster. */
void ClusterNode::MergeFrames( ClusterNode const& rhs) {
  frameList_.insert(frameList_.end(), rhs.frameList_.begin(), rhs.frameList_.end());
}

void ClusterNode::SetName(std::string const& nameIn, double rmsIn) {
  name_ = nameIn;
  refRms_ = rmsIn;
}
#endif
