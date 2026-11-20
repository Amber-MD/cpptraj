#ifdef HAS_EIGEN
#include "result_utils.h"
#include "scores.h"

#include <map>
#include <stdexcept>
#include <string>

using namespace Cpptraj::Mdance;

std::vector<HCTree> Cpptraj::Mdance::buildClusterTree(const ArrayXXd& data, const Eigen::ArrayXi& labels) {
    // This has to be checked here, before the accumulation loop below indexes
    // data.row(j) for every label: a label array longer than the coordinate
    // matrix reads past the end of it, and Eigen does not bounds-check in a
    // release build. A mismatch means the two inputs came from different
    // trajectories, so there is nothing sensible to do but say so.
    if (labels.size() != data.rows()) {
        throw std::runtime_error(
            "initial labels (" + std::to_string(labels.size()) +
            ") do not match the number of frames (" + std::to_string(data.rows()) + ").");
    }

    std::set<int> uniqueLabels;
    for (int i = 0; i < labels.size(); i++) {
        uniqueLabels.insert(labels(i));
    }

    std::vector<HCTree> clusters;
    int clusterIdx = 0;
    for (int label : uniqueLabels) {
        // The tree's index list holds ORIGINAL label values (see the contract on
        // HCTree::insertRoot), not the dense position -- callers map the final
        // clusters back to frames by matching these against the label column.
        // Storing clusterIdx here instead silently mislabels every frame unless
        // the label set happens to be exactly {0..K-1}.
        // The z_index stays dense: it numbers the Z-matrix leaves.
        std::list<int> indices = {label};
        Vec cSumi = Vec::Zero(data.cols());
        Vec sqSumi = Vec::Zero(data.cols());
        int ni = 0;
        for (int j = 0; j < labels.size(); j++) {
            if (labels(j) == label) {
                ni++;
                cSumi += data.row(j);
                sqSumi += data.row(j).square();
            }
        }
        HCTree tree = HCTree();
        tree.insertRoot(indices, cSumi, sqSumi, ni, clusterIdx);
        clusters.push_back(tree);
        clusterIdx++;
    }

    return clusters;
}

std::vector<int> Cpptraj::Mdance::labelsFromHelmClusters(const std::vector<HCTree>& finalClusters,
                                         const Eigen::ArrayXi& initialLabels,
                                         int nFrames) {
    // HELM is seeded from a per-frame partition, so a label file that does not
    // line up with the coordinate file is a mismatched pair of inputs, not
    // something to paper over: the old code wrote past the end of the label
    // vector when there were more labels than frames.
    if ((int)initialLabels.size() != nFrames) {
        throw std::runtime_error(
            "initial labels (" + std::to_string(initialLabels.size()) +
            ") do not match the number of frames (" + std::to_string(nFrames) + ").");
    }

    // Invert once: original initial-cluster label -> final HELM cluster. Beats
    // rescanning every frame for every (final cluster, initial cluster) pair.
    std::map<int, int> labelToFinal;
    for (int c = 0; c < (int)finalClusters.size(); ++c) {
        // getRootClusterIndices is non-const; the tree copy is cheap (it holds
        // pointers) and this keeps the caller's vector untouched.
        HCTree tree = finalClusters[c];
        for (int initialLabel : tree.getRootClusterIndices()) {
            labelToFinal[initialLabel] = c;
        }
    }

    std::vector<int> labels(nFrames, -1);
    for (int j = 0; j < nFrames; ++j) {
        auto it = labelToFinal.find(initialLabels(j));
        if (it != labelToFinal.end()) labels[j] = it->second;
    }
    return labels;
}

std::vector<int> Cpptraj::Mdance::computeRepresentatives(const ArrayXXd& data, const std::vector<int>& labels,
                                         int nClusters, int nAtoms, MD::Metric mt) {
    std::vector<int> reps;
    for (int c = 0; c < nClusters; ++c) {
        std::vector<int> memberIndices;
        for (size_t i = 0; i < labels.size(); ++i) {
            if (labels[i] == c) memberIndices.push_back(i);
        }
        if (memberIndices.empty()) {
            reps.push_back(-1);
            continue;
        }
        ArrayXXd subData(memberIndices.size(), data.cols());
        for (size_t i = 0; i < memberIndices.size(); ++i) {
            subData.row(i) = data.row(memberIndices[i]);
        }
        Index medoidLocal = calculateMedoid(subData, nAtoms, mt);
        reps.push_back(memberIndices[medoidLocal]);
    }
    return reps;
}

std::vector<int> Cpptraj::Mdance::computeClusterSizes(const std::vector<int>& labels, int nClusters) {
    std::vector<int> sizes(nClusters, 0);
    for (int l : labels) {
        if (l >= 0 && l < nClusters) sizes[l]++;
    }
    return sizes;
}

std::vector<double> Cpptraj::Mdance::computeClusterMSD(const ArrayXXd& data, const std::vector<int>& labels,
                                       int nClusters, int nAtoms) {
    std::vector<double> msds;
    for (int c = 0; c < nClusters; ++c) {
        std::vector<int> members;
        for (size_t i = 0; i < labels.size(); ++i) {
            if (labels[i] == c) members.push_back(i);
        }
        if (members.size() < 2) {
            msds.push_back(0.0);
            continue;
        }
        ArrayXXd subData(members.size(), data.cols());
        for (size_t i = 0; i < members.size(); ++i) {
            subData.row(i) = data.row(members[i]);
        }
        msds.push_back(meanSqDev(subData, nAtoms));
    }
    return msds;
}
#endif /* HAS_EIGEN */
