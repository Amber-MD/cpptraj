//this cluster h file is to help define clusterDict object in helm.cpp
#ifndef INC_CPPTRAJ_MDANCE_CLUSTER_H
#define INC_CPPTRAJ_MDANCE_CLUSTER_H
#ifdef HAS_EIGEN
#include "types.h"
namespace Cpptraj {
namespace Mdance {
class Cluster{
    private:
        /*
        |   indices: cluster indices of merged clusters
        |   c_sum: feature array of the column-wise sum of data
        |   sq_sum: feature array of the column-wise of squared data
        |   n: number of elements in cluster

        */
        Veci indices;
        Vec cSum; 
        Vec sqSum;
        int n;
        Mat cluster;
        
    public:
        Cluster();
        Cluster(Veci indices, Vec cSum, Vec sqSum, int n);
        Cluster(Veci indices, Vec cSum, Vec sqSum, int n, Mat cluster);
        Veci getIndices();
        Vec getCsum();
        Vec getSQsum();
        int getN();
        Mat getCluster();
};
}
}
#endif
#endif
