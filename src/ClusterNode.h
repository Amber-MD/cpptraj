#ifndef INC_CLUSTERNODE_H
#define INC_CLUSTERNODE_H
#include "ClusterDist.h" 
/// Hold information for a cluster in a ClusterList
class ClusterNode {
  public:
    ClusterNode();
    ~ClusterNode();
    ClusterNode(ClusterDist*,ClusterDist::Cframes const&, int);
    ClusterNode(const ClusterNode&);
    ClusterNode& operator=(const ClusterNode&);
    /// Used to sort clusters by # of frames in cluster
    bool operator<(const ClusterNode&) const;
    /// Merge frames from another cluster to this cluster
    void MergeFrames(ClusterNode&);
    /// Determine which frame in the cluster is centroid.
    int FindCentroidFrame(ClusterMatrix const&);
    /// Calculate eccentricity for frames in this cluster.
    void CalcEccentricity(ClusterMatrix const&);
    /// Calculate centroid of members of this cluster.
    void CalculateCentroid(ClusterDist* Cdist) {
      Cdist->CalculateCentroid( centroid_, frameList_ );
    }
    /// Calculate average distance of all members to centroid
    double CalcAvgToCentroid( ClusterDist*);
    // Iterator over frame numbers
    typedef ClusterDist::Cframes::const_iterator frame_iterator;
    const frame_iterator beginframe() const { return frameList_.begin(); }
    const frame_iterator endframe()   const { return frameList_.end();   }
    // Return internal variables
    inline double AvgDist()      const { return avgClusterDist_;        }
    inline double Eccentricity() const { return eccentricity_;          }
    inline int Num()             const { return num_;                   }
    inline int Nframes()         const { return (int)frameList_.size(); }
    inline int CentroidFrame()   const { return centroidframe_;         }
    inline Centroid* Cent()            { return centroid_;              }
    // Set internal variables 
    void SetAvgDist(double avg)        { avgClusterDist_ = avg; }
    void AddFrameToCluster(int fnum)   { frameList_.push_back( fnum ); }
    void SetNum(int);
  private:
    double avgClusterDist_;           ///< Avg distance of this cluster to all other clusters.
    double eccentricity_;             ///< Maximum distance between any 2 frames.
    int num_;                         ///< Cluster number.
    int centroidframe_;               ///< Frame number with lowest dist. to all other frames.
    ClusterDist::Cframes frameList_;  ///< List of frames belonging to this cluster.
    Centroid* centroid_;              ///< Centroid of all frames in this cluster. 
};
#endif
