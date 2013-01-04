#ifndef INC_CLUSTERNODE_H
#define INC_CLUSTERNODE_H
#include <list>
#include "ClusterMatrix.h"
#include "Frame.h"
/// Hold information for a cluster in a ClusterList
class ClusterNode {
  public:
    // For holding RMS/DME calc options
    struct RMSoptions {
      bool useDME;
      bool nofit;
      bool useMass;
      AtomMask mask;
    };
    ClusterNode();
    ClusterNode(std::list<int> const&, int);
    ClusterNode(DataSet*, int, RMSoptions const&);
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
    /// Calculate average distance between all frames in the cluster
    void CalcAvgFrameDist(ClusterMatrix const&);
    /// Calculate centroid of members of this cluster.
    void CalculateCentroid(DataSet*, RMSoptions const&);
    /// Calculate average distance of all members to centroid
    double CalcAvgToCentroid( DataSet*, RMSoptions const&);
    /// Calculate distance from this centroid to another nodes centroid.
    double CentroidDist( ClusterNode const&, RMSoptions const& );
    // Iterator over frame numbers
    typedef std::list<int>::const_iterator frame_iterator;
    const frame_iterator beginframe() const { return frameList_.begin(); }
    const frame_iterator endframe()   const { return frameList_.end();   }
    // Return internal variables
    inline double AvgDist()      const { return avgClusterDist_;        }
    inline double InternalAvg()  const { return internalAvg_;           }
    inline double InternalSD()   const { return internalSD_;            }
    inline double Eccentricity() const { return eccentricity_;          }
    inline int Num()             const { return num_;                   }
    inline int Centroid()        const { return centroid_;              }
    inline int Nframes()         const { return (int)frameList_.size(); }
    // Set internal variables 
    void SetAvgDist(double avg)        { avgClusterDist_ = avg; }
    void AddFrameToCluster(int fnum)   { frameList_.push_back( fnum ); }
    void SetNum(int);
    void FrameSieveOffset(int);
  private:
    double avgClusterDist_;    ///< Avg distance of this cluster to all other clusters.
    double internalAvg_;
    double internalSD_;
    double eccentricity_;      ///< Maximum distance between any 2 frames.
    double cval_;              ///< Centroid value (avg) for non-coords DataSets.
    int num_;                  ///< Cluster number.
    int centroid_;             ///< Frame number with lowest distance to all other frames.
    std::list<int> frameList_; ///< List of frames belonging to this cluster.
    Frame cframe_;             ///< Centroid frame (avg) coords.
};
#endif
