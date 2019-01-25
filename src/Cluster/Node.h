#ifndef INC_CLUSTER_NODE_H
#define INC_CLUSTER_NODE_H
#include "PairwiseMatrix.h"
namespace Cpptraj {
namespace Cluster {

// TODO implement needsUpdate_

/// Hold frame indices for a given cluster.
class Node {
  public:
    Node();
    ~Node();
    // NOTE: Taking a pointer to Metric here instead of a reference allows
    //       PairwiseMatrix to be passed in as const to routines while still
    //       allowing Metric to be used. Metric needs to be non-const because
    //       things like calculating RMSD modify Metric itself to avoid
    //       always reallocating Frames.
    /// CONSTRUCTOR - Take Metric for calculating centroid, frames, and cluster index.
    Node(Metric*, Cframes const&, int);
    /// COPY CONSTRUCTOR
    Node(const Node&);
    /// ASSIGNMENT
    Node& operator=(const Node&);

    /// Used to pair a representative frame number with a score.
    typedef std::pair<int,float> RepPair;
    /// Used to hold a list of representative frames/scores
    typedef std::vector<RepPair> RepPairArray;
    /// Used to pair frame numbers with silhouette values.
    typedef std::pair<int,float> SilPair;
    /// Used to hold list of frame numbers/silhouette values.
    typedef std::vector<SilPair> SilPairArray;

    /// Used to sort clusters by # of frames in cluster
    inline bool operator<(const Node&) const;
    /// Find and set frame in the cluster that has lowest distance to all other frames.
//    int SetBestRep_CumulativeDist(DataSet_Cmatrix const&);
    /// Calculate average distance of all members to centroid
    double CalcAvgToCentroid( Metric*) const;
    /// Const iterator over frame numbers
    typedef Cframes::const_iterator frame_iterator;
    /// Const iterator to beginning of frames
    frame_iterator beginframe()    const { return frameList_.begin();     }
    /// Const iteratir to end of frames.
    frame_iterator endframe()      const { return frameList_.end();       }
    /// \return Frame number at given index.
    int ClusterFrame(int idx)      const { return frameList_[idx];        }
    /// \return cluster eccentricity
    double Eccentricity()          const { return eccentricity_;          }
    /// \return internal cluster number
    int Num()                      const { return num_;                   }
    /// \return number of frames in cluster.
    int Nframes()                  const { return (int)frameList_.size(); }
    /// \return best representative frame number, or -1 if no best rep set.
    int BestRepFrame()             const {
      if (bestReps_.empty())
        return -1;
      else
        return bestReps_.front().first;
    }
    /// \return Cluster centroid.
    Centroid* Cent()               const { return centroid_;              }
    /// \return name assigned via reference
    std::string const& Cname()     const { return name_;                  }
    /// \return RMS to reference
    double RefRms()                const { return refRms_;                }
    /// \return true if given frame is in this cluster.
    bool HasFrame(int f)           const { return frameList_.HasFrame(f); }
    /// Access representative frame list, const
    RepPairArray const& BestReps() const { return bestReps_;              }
    /// Access frame silhouette list
    SilPairArray const& FrameSilhouettes() const { return frameSil_;      }
    /// \return cluster silhoueete vaule
    double Silhouette()            const { return avgSil_;                }

    /// Calculate centroid of members of this cluster.
    void CalculateCentroid(Metric* Cdist) {
      // FIXME: Could potentially get rid of this branch.
      if (centroid_ == 0)
        centroid_ = Cdist->NewCentroid( frameList_ );
      else
        Cdist->CalculateCentroid( centroid_, frameList_ );
    }
    /// Add frame to cluster
    void AddFrameToCluster(int fnum)   { frameList_.push_back( fnum );  }
    /// Set cluster number (for bookkeeping).
    void SetNum(int numIn)             { num_ = numIn;                  }
    /// Access representative frame list
    RepPairArray& BestReps()           { return bestReps_;              }
    /// Access frame silhouette list
    SilPairArray& FrameSilhouettes()   { return frameSil_;              }
    /// Set cluster silhouette value
    void SetSilhouette(double s)       { avgSil_ = s;                   }
    /// Sort internal frame list
    void SortFrameList()               { frameList_.Sort();             }
    /// Remove specified frame from cluster if present.
    void RemoveFrameFromCluster(int f) { frameList_.Remove(f);          }
    /// Merge frames from another cluster to this cluster
    inline void MergeFrames(Node const&);
    /// Set cluster name and RMS to reference
    inline void SetNameAndRms(std::string const&, double);
    /// Calculate eccentricity for frames in this cluster.
    void CalcEccentricity(PairwiseMatrix const&);
    /// Remove specified frame from cluster and update centroid.
    void RemoveFrameUpdateCentroid(Metric*, int);
    /// Add specified frame to cluster and update centroid.
    void AddFrameUpdateCentroid(Metric*, int);
  private:
    Cframes frameList_;     ///< List of frames belonging to this cluster.
    Centroid* centroid_;    ///< Centroid of all frames in this cluster.
    std::string name_;      ///< Cluster name assigned from reference.
    RepPairArray bestReps_; ///< Hold best representative frames and their score.
    SilPairArray frameSil_; ///< Frame silhouette values.
    double avgSil_;         ///< Average silhouette value for cluster TODO s.d. as well?
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
  frameList_.Insert( rhs.frameList_ );
  //frameList_.insert(frameList_.end(), rhs.frameList_.begin(), rhs.frameList_.end());
}
/** Set reference name and RMS. */
void Node::SetNameAndRms(std::string const& nameIn, double rmsIn) {
  name_ = nameIn;
  refRms_ = rmsIn;
}

} /* END namespace Cluster */
} /* END namespace Cpptraj */
#endif
