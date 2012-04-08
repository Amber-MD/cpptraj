#ifndef INC_CLUSTERNODE_H
#define INC_CLUSTERNODE_H
#include <list>
class ClusterNode {
  public:
    ClusterNode();
    ClusterNode(std::list<int>&, int);
    ClusterNode(const ClusterNode&);
    ClusterNode& operator=(const ClusterNode&);

    //bool operator==(const ClusterNode&);
    //bool operator!=(const ClusterNode&);
    bool operator<(const ClusterNode&);
    //bool operator>(const ClusterNode&);

    void MergeFrames(ClusterNode&);

    typedef std::list<int>::const_iterator frame_iterator;
    inline const frame_iterator beginframe() const { 
      return frameList_.begin(); 
    }
    inline const frame_iterator endframe() const { return frameList_.end(); }

    inline double AvgDist() { return avgClusterDist_; }
    inline double Eccentricity() { return eccentricity_; }
    inline const int Num() const { return num_; }
    inline const int Centroid() const { return centroid_; }
    inline int Nframes() { return (int)frameList_.size(); }

    void SetAvgDist(double avg) { avgClusterDist_ = avg; }
    void SetCentroid(int cent) { centroid_ = cent; }
    void SetEccentricity(double eccen) { eccentricity_ = eccen; }
    void SetNum(int);

  private:
    double avgClusterDist_;
    double eccentricity_;
    int num_;
    int centroid_;
    std::list<int> frameList_;
};
#endif
