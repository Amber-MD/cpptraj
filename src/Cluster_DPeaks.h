#ifndef INC_CLUSTER_DPEAKS_H
#define INC_CLUSTER_DPEAKS_H
#include "ClusterList.h"
class Cluster_DPeaks : public ClusterList {
  public:
    Cluster_DPeaks();
    static void Help();
    int SetupCluster(ArgList&);
    void ClusteringInfo();
    int Cluster();
    void AddSievedFrames();
    void ClusterResults(CpptrajFile&) const;
  private:
    double epsilon_;
    typedef std::vector<int> Iarray;
    typedef std::vector<double> Darray;
};
#endif
