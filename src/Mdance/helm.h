#ifndef INC_CPPTRAJ_MDANCE_HELM_H
#define INC_CPPTRAJ_MDANCE_HELM_H
#ifdef HAS_EIGEN
#include "types.h"
#include "hc_utils.h"
namespace Cpptraj {
namespace Mdance {
class Helm{
    public:
        Helm(vector<HCTree> clusterTree, int nAtoms, MD::Metric mt = MD::Metric::MSD, 
                MD::MergeScheme mergeScheme = MD::MergeScheme::Inter, int nClusters = 0, float eps = -1, 
                bool trimStart = false,
                float minSamples = 0.01,
                float trimVal=0, float trimK=0,
                bool savePairwiseSum = false,
                string inputTop ="", string inputTraj ="");
        vector<HCTree> run();
        pair<double, double> computeScores(vector<HCTree> clusters, Mat data);
        Mat getZMatrix();
    private:
        vector<HCTree> clusterTree;
        int nAtoms;
        int nClusters;
        float eps; // -1 means None 
        bool trimStart;
        float trimVal;
        float trimK;
        float minSamples;
        int trimIncoming;
        MD::Metric mt;
        MD::MergeScheme mergeScheme;
        Mat clusterDists;
        int totalIncoming;
        bool savePairwiseSum;
        string inputTop;
        string inputTraj;
        int totalSum;
        Mat zMatrix;

        Mat makeDataByRow(Vec a, Vec b);
        vector<HCTree> trimClusters();
        vector<HCTree> genNewClusters(int ZIdx);
        float calcHelmSim(HCTree& firstTree, HCTree& secondTree);
        void genClusterDists(vector<HCTree>& previousClusters);
        void updateZMatrix(int idxA, int idxB, int mergedClusts);
};
} /* END namespace Mdance */
} /* END namespace Cpptraj */
#endif /* HAS_EIGEN */
#endif
