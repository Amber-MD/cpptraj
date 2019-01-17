#ifndef INC_CLUSTER_NODE_H
#define INC_CLUSTER_NODE_H
#include "PairwiseMatrix.h"
namespace Cpptraj {
namespace Cluster {

/// Hold frame indices for a given cluster.
class Node {
  public:
    Node();
    ~Node();
    Node(Metric&, Cframes const&, int);
    Node(const Node&);
    Node& operator=(const Node&);
    /// Used to pair a representative frame number with a score.
    typedef std::pair<int,float> RepPair;
    /// Used to hold a list of representative frames/scores
    typedef std::vector<RepPair> RepPairArray;
    /// Used to sort clusters by # of frames in cluster
    inline bool operator<(const Node&) const;
    /// Merge frames from another cluster to this cluster
    inline void MergeFrames(Node const&);
    /// Find and set frame in the cluster that has lowest distance to all other frames.
//    int SetBestRep_CumulativeDist(DataSet_Cmatrix const&);
    /// Calculate eccentricity for frames in this cluster.
    void CalcEccentricity(PairwiseMatrix const&);
    /// Calculate centroid of members of this cluster.
    void CalculateCentroid(Metric& Cdist) {
      // FIXME: Could potentially get rid of this branch.
      if (centroid_ == 0)
        centroid_ = Cdist.NewCentroid( frameList_ );
      else
        Cdist.CalculateCentroid( centroid_, frameList_ );
    }
    /// Calculate average distance of all members to centroid
    double CalcAvgToCentroid( Metric&) const;
    // Iterator over frame numbers
    typedef Cframes::const_iterator frame_iterator;
    frame_iterator beginframe() const { return frameList_.begin(); }
    frame_iterator endframe()   const { return frameList_.end();   }
    /// \return Frame number at given index.
    int ClusterFrame(int idx)   const { return frameList_[idx];    }
    /// \return cluster eccentricity
    double Eccentricity()      const { return eccentricity_;          }
    /// \return internal cluster number
    int Num()                  const { return num_;                   }
    /// \return number of frames in cluster.
    int Nframes()              const { return (int)frameList_.size(); }
    /// \return best representative frame number, or -1 if no best rep set.
    int BestRepFrame()         const {
      if (bestReps_.empty())
        return -1;
      else
        return bestReps_.front().first;
    }
    Centroid* Cent()           const { return centroid_;              }
    std::string const& Cname() const { return name_;                  }
    double RefRms()            const { return refRms_;                }
    // Set internal variables 
    void AddFrameToCluster(int fnum) { frameList_.push_back( fnum );  }
    void SetNum(int numIn)           { num_ = numIn;                  }
    /// Access representative frame list
    RepPairArray const& BestReps() const { return bestReps_; }
    RepPairArray&       BestReps()       { return bestReps_; }
    inline void SetNameAndRms(std::string const&, double);
    /// Sort internal frame list
    void SortFrameList();
    /// \return true if given frame is in this cluster.
    bool HasFrame(int) const;
    /// Remove specified frame from cluster if present.
    void RemoveFrameFromCluster(int);
    /// Remove specified frame from cluster and update centroid.
    void RemoveFrameUpdateCentroid(Metric&, int);
    /// Add specified frame to cluster and update centroid.
    void AddFrameUpdateCentroid(Metric&, int);
  private:
    Cframes frameList_;     ///< List of frames belonging to this cluster.
    Centroid* centroid_;    ///< Centroid of all frames in this cluster.
    std::string name_;      ///< Cluster name assigned from reference.
    RepPairArray bestReps_;
    double eccentricity_;   ///< Maximum distance between any 2 frames.
    double refRms_;         ///< Cluster rms to reference (if assigned).
    int num_;               ///< Cluster number, used for bookkeeping.
    bool needsUpdate_;      ///< True if internal metrics need updating (e.g. after frames added).
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
/** Use > since we give higher priority to larger clusters. */
bool Node::operator<(const Node& rhs) const {
  return ( frameList_.size() > rhs.frameList_.size() );
}
/** Frames from rhs go to this cluster. */
void Node::MergeFrames( Node const& rhs) {
  frameList_.insert(frameList_.end(), rhs.frameList_.begin(), rhs.frameList_.end());
}
/** Set reference name and RMS. */
void Node::SetNameAndRms(std::string const& nameIn, double rmsIn) {
  name_ = nameIn;
  refRms_ = rmsIn;
}

} /* END namespace Cluster */
} /* END namespace Cpptraj */
#endif
