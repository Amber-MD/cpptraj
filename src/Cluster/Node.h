#ifndef INC_CLUSTER_NODE_H
#define INC_CLUSTER_NODE_H
#include <string>
#include <utility> // std::pair
#include <vector>
#include "CentroidArray.h"
#include "Cframes.h" // Cframes::const_iterator
class DataSet_integer;
class DataSet_float;
namespace Cpptraj {
namespace Cluster {
class MetricArray;
// TODO implement needsUpdate_

/// Hold frame indices for a given cluster.
class Node {
  public:
    Node();
    ~Node();
    // NOTE: Taking a non-const reference to Metric here allows
    //       MetricArray to be passed in as const to routines while still
    //       allowing Metric to be used. Metric needs to be non-const because
    //       things like calculating RMSD modify Metric itself to avoid
    //       always reallocating Frames.
    /// CONSTRUCTOR - Take Metric for calculating centroid, frames, and cluster index.
    Node(MetricArray&, Cframes const&, int);
    /// COPY CONSTRUCTOR
    Node(const Node&);
    /// ASSIGNMENT
    Node& operator=(const Node&);

    /// Types of normalization for cluster pop v time.
    enum CnormType { NONE = 0, CLUSTERPOP, FRAME };

    /// Used to pair a representative frame number with a score.
    typedef std::pair<int,double> RepPair;
    /// Used to hold a list of representative frames/scores
    typedef std::vector<RepPair> RepPairArray;
    /// Used to pair frame numbers with silhouette values.
    typedef std::pair<int,double> SilPair;
    /// Used to hold list of frame numbers/silhouette values.
    typedef std::vector<SilPair> SilPairArray;

    /// Used to sort clusters by # of frames in cluster
    inline bool operator<(const Node&) const;
    /// Find and set frame in the cluster that has lowest distance to all other frames.
//    int SetBestRep_CumulativeDist(DataSet_Cmatrix const&);
    /// Calculate average distance of all members to centroid
    double CalcAvgToCentroid(MetricArray&) const;
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
    /// \return Cluster centroid array.
    CentroidArray const& Cent()    const { return centroids_;             }
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
    void CalculateCentroid(MetricArray&);
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
    void CalcEccentricity(MetricArray&);
    /// Remove specified frame from cluster and update centroid.
    void RemoveFrameUpdateCentroid(MetricArray&, int);
    /// Add specified frame to cluster and update centroid.
    void AddFrameUpdateCentroid(MetricArray&, int);

    /// Calculate cluster population vs time.
    void CalcCpopVsTime(DataSet_float&, unsigned int, CnormType) const;
    /// Create cluster lifetime set.
    void CreateLifetimeSet(DataSet_integer&, unsigned int) const;
  private:
    Cframes frameList_;       ///< List of frames belonging to this cluster.
    CentroidArray centroids_; ///< Centroids (1 for each metric) for all frames in this cluster.
    std::string name_;        ///< Cluster name assigned from reference.
    RepPairArray bestReps_;   ///< Hold best representative frames and their score.
    SilPairArray frameSil_;   ///< Frame silhouette values.
    double avgSil_;           ///< Average silhouette value for cluster TODO s.d. as well?
    double eccentricity_;     ///< Maximum distance between any 2 frames.
    double refRms_;           ///< Cluster rms to reference (if assigned).
    int num_;                 ///< Cluster number, used for bookkeeping.
    bool needsUpdate_;        ///< True if internal metrics need updating (e.g. after frames added).
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
