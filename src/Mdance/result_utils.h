#ifndef INC_CPPTRAJ_MDANCE_RESULT_UTILS_H
#define INC_CPPTRAJ_MDANCE_RESULT_UTILS_H
#ifdef HAS_EIGEN
#include <vector>
#include "types.h"
#include "bts.h"
#include "hc_utils.h"
namespace Cpptraj {
namespace Mdance {
std::vector<HCTree> buildClusterTree(const ArrayXXd& data, const Eigen::ArrayXi& labels);

// Map the clusters HELM converged on back to one label per frame. Throws when
// initialLabels and nFrames disagree, which means the coordinate and label
// files were not produced from the same trajectory.
std::vector<int> labelsFromHelmClusters(const std::vector<HCTree>& finalClusters,
                                        const Eigen::ArrayXi& initialLabels,
                                        int nFrames);

std::vector<int> computeRepresentatives(const ArrayXXd& data, const std::vector<int>& labels,
                                         int nClusters, int nAtoms, MD::Metric mt);

std::vector<int> computeClusterSizes(const std::vector<int>& labels, int nClusters);

std::vector<double> computeClusterMSD(const ArrayXXd& data, const std::vector<int>& labels,
                                       int nClusters, int nAtoms);
}
}
#endif /* HAS_EIGEN */
#endif
