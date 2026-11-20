#include "cluster.h"
#ifdef HAS_EIGEN
using namespace Cpptraj::Mdance;
Cluster::Cluster() : n(0) {
}
Cluster::Cluster(Veci indices, Vec cSum, Vec sqSum, int n){
    this->indices = indices;
    this->cSum = cSum;
    this->sqSum = sqSum;
    this->n = n;
}
Cluster::Cluster(Veci indices, Vec cSum, Vec sqSum, int n, Mat cluster){
    this->indices = indices;
    this->cSum = cSum;
    this->sqSum = sqSum;
    this->n = n;
    this->cluster = cluster;
}

//getter functions
Veci Cluster::getIndices(){
    return indices;
}
Vec Cluster::getCsum(){
    return cSum;
}
Vec Cluster::getSQsum(){
    return sqSum;
}
int Cluster::getN(){
    return n;
}
Mat Cluster::getCluster(){
    return cluster;
}
#endif
