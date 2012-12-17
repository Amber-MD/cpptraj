#ifndef INC_CLUSTERNODE_H
#define INC_CLUSTERNODE_H
#include <list>
#include "DataSet_Coords.h"
#include "TriangleMatrix.h"
/// Hold information for a cluster in a ClusterList
class ClusterNode {
  public:
    ClusterNode();
    ClusterNode(std::list<int> const&, int);
    ClusterNode(const ClusterNode&);
    ClusterNode& operator=(const ClusterNode&);
    /// Used in sorting
    bool operator<(const ClusterNode&) const;
    /// Merge frames from another cluster to this cluster
    void MergeFrames(ClusterNode&);
    /// Calculate centroid frame from frames in this cluster.
    void CalcCentroidFrame(DataSet_Coords const&, AtomMask const&);
    /// Determine which frame in the cluster is centroid.
    int FindCentroid(TriangleMatrix const&);
    /// Calculate eccentricity for frames in this cluster.
    void CalcEccentricity(TriangleMatrix const&);
    /// Calculate average distance between all frames in the cluster
    void CalcAvgFrameDist(TriangleMatrix const&);
    /// Calculate average distance of all frames to centroid
    double CalcAvgToCentroid( DataSet_Coords const&, AtomMask const& );
    /// Calculate distance from this centroid to another nodes centroid.
    double CentroidDist( ClusterNode const& );
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
    int num_;                  ///< Cluster number.
    int centroid_;             ///< Frame number with lowest distance to all other frames.
    std::list<int> frameList_; ///< List of frames belonging to this cluster.
    Frame cframe_;             ///< Centroid frame (avg) coords.
};
#endif
