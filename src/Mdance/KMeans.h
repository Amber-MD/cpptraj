#ifndef CPPTRAJ_MDANCE_KMEANS_H
#define CPPTRAJ_MDANCE_KMEANS_H
#ifdef HAS_EIGEN

//#include <iostream>
//#include <stdexcept>
//#include <limits>

#include "bts.h"
#include "types.h"
#include "scores.h"
#include "../Random.h"

namespace Cpptraj {
namespace Mdance {
class KmeansNANI{
    Mat data;
    Mat centers;
    Mat dist;
    Veci labels;
    MD::KinitType kinit;
    int seed;
    int kClusters;
    MD::Metric mt;
    int nAtoms;
    int percentage;
    int vectorizationThreshold;
    Random_Number RNG_;

    void set_seed();
    void set_vectorization_threshold(int threshold);
    int randint(int low, int high);
    int discrete_rand(Vec &p);
    void select_without_replacement(int N, int K, Vec &chosenIDs);
    void sampleRowsRandom();
    void sampleRowsPlusPlus();
    void reduced_init_Mu(bool isComp);
    void init_Mu();
    void pairwise_distance(Mat &X, Mat &Mu, Mat &Dist);
    double assignClosest();
    void calcMu();
    void run_lloyd(int Niter);

public: 
    KmeansNANI(ArrayXXd data, int kClusters, MD::Metric mt, MD::KinitType kinit = MD::KinitType::StratAll, int nAtoms = 1, int percentage = 10, int vectThreshold=16, int seedIn=0);
    KmeansNANI(ArrayXXd data, int kClusters, MD::Metric mt, Mat centers, int nAtoms = 1, int percentage = 10, int vectThreshold=16, int seedIn=0);
    map<int,vector<Index>> createClusterDict();
    pair<double, double> computeScores();
    Veci getLabels();
    Mat getCenters();
};
} /** END namespace Mdance */
} /* END namespace Cpptraj */
#endif /* HAS_EIGEN */
#endif
