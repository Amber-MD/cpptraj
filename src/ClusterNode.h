#ifndef INC_CLUSTERNODE_H
#define INC_CLUSTERNODE_H
#include <list>
/// Hold information for a cluster in a ClusterList
class ClusterNode {
  public:
    ClusterNode();
    ClusterNode(std::list<int>&, int);
    ClusterNode(const ClusterNode&);
    ClusterNode& operator=(const ClusterNode&);
    /// Used in sorting
    bool operator<(const ClusterNode&) const;
    /// Merge frames from another cluster to this cluster
    void MergeFrames(ClusterNode&);
    // Iterator over frames
    typedef std::list<int>::const_iterator frame_iterator;
    inline const frame_iterator beginframe() const { 
      return frameList_.begin(); 
    }
    inline const frame_iterator endframe() const { return frameList_.end(); }
    // Return internal variables
    inline double AvgDist() { return avgClusterDist_; }
    inline double Eccentricity() { return eccentricity_; }
    inline int Num() const { return num_; }
    inline int Centroid() const { return centroid_; }
    inline int Nframes() { return (int)frameList_.size(); }
    // Set internal variables 
    void SetAvgDist(double avg) { avgClusterDist_ = avg; }
    void SetEccentricity(double eccen) { eccentricity_ = eccen; }
    void SetNum(int);
    void SetCentroid(int cent) { centroid_ = cent; }

  private:
    double avgClusterDist_;
    double eccentricity_;
    int num_;
    int centroid_;
    std::list<int> frameList_;
};
#endif
