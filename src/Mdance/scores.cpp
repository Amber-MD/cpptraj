#include "scores.h"
#ifdef HAS_EIGEN
using namespace Cpptraj::Mdance;
/* Returns the Calinski-Harabasz score for a given cluster array.
 * 
 * Parameters:
 *  - data: A feature array of size (nSamples, nFeatures)
 *  - labels: A list of labels.
 * 
 * Returns: the CH score
 * 
 * Reference: This code does not include any of the parameter checking: https://github.com/scikit-learn/scikit-learn/blob/c5497b7f7eacfaff061cf68e09bcd48aa93d4d6b/sklearn/metrics/cluster/_unsupervised.py#L325
*/
double Cpptraj::Mdance::calinskiHarabaszScore(const ArrayXXd& data, const VectorXi& labels) {
    std::map<int, vector<int>> clusters;

    for (int i=0; i<labels.size(); ++i) {
        if (clusters.find(labels[i]) == clusters.end()){
            vector<int> k = {i};
            clusters[labels[i]] = k;
        } else {
            clusters[labels[i]].emplace_back(i);
        }
    }

    // With a single cluster there is no between-cluster dispersion and the
    // ratio below is 0/0. sklearn rejects k < 2 outright; returning 0 keeps the
    // score a real number for callers that serialize it.
    if (clusters.size() < 2)
        return 0;

    double extraDisp = 0;
    double intraDisp = 0;
    ArrayXd mean = data.colwise().sum() / data.rows();
    
    for (auto k=clusters.begin(); k != clusters.end(); ++k) {
        ArrayXXd cluster = data(k->second,Eigen::all);
        ArrayXd clusterMean = cluster.colwise().sum() / cluster.rows();

        extraDisp += cluster.rows() * (clusterMean - mean).square().sum();

        for (int i=0; i<cluster.rows(); ++i) {
            intraDisp += (cluster.row(i) - clusterMean.transpose()).square().sum();
        }
    }

    if (intraDisp == 0)
        return 1;
    
    return extraDisp * (data.rows() - clusters.size()) / (intraDisp * (clusters.size() - 1.0));
}

/* Returns the Davies-Bouldin score for a given cluster array.
 * 
 * Parameters:
 *  - data: A feature array of size (nSamples, nFeatures)
 *  - labels: A list of labels.
 * 
 * Returns: the DB score
 * 
 * Reference: This code does not include any of the parameter checking: https://github.com/scikit-learn/scikit-learn/blob/c5497b7f7eacfaff061cf68e09bcd48aa93d4d6b/sklearn/metrics/cluster/_unsupervised.py#L396
*/
double Cpptraj::Mdance::daviesBouldinScore(const ArrayXXd& data, const VectorXi& labels) {
    std::map<int, vector<int>> clusters;

    for (int i=0; i<labels.size(); ++i) {
        if (clusters.find(labels[i]) == clusters.end()){
            vector<int> k = {i};
            clusters[labels[i]] = k;
        } else {
            clusters[labels[i]].emplace_back(i);
        }
    }

    ArrayXXd centroids(clusters.size(), data.row(0).size());
    ArrayXd intraDists(clusters.size());
    ArrayXXd scores(clusters.size(), clusters.size());
    ArrayXd maxScores(clusters.size());
    int c=0;

    for (auto k=clusters.begin(); k != clusters.end(); ++k) {
        ArrayXXd cluster = data(k->second,Eigen::all);
        centroids.row(c) = cluster.colwise().sum() / cluster.rows();
        intraDists[c] = 0;
        for (int i=0; i<cluster.rows(); ++i) {
            intraDists[c] += std::sqrt((cluster.row(i) - centroids.row(c)).square().sum());
        }
        intraDists[c] /= cluster.rows();

        scores.row(c)[c] = 0;
        maxScores[c] = 0;
        for (int i=0; i<c; ++i) {
            scores.row(c)[i] = (intraDists[i] + intraDists[c]) / std::sqrt((centroids.row(i) - centroids.row(c)).square().sum());
            scores.row(i)[c] = scores.row(c)[i];

            if (scores.row(i)[c] > maxScores[i])
                maxScores[i] = scores.row(i)[c];
            if (scores.row(c)[i] > maxScores[c])
                maxScores[c] = scores.row(c)[i];
        }

        ++c;
    }    
    
    return maxScores.sum() / maxScores.size();
}
#endif /* HAS_EIGEN */
