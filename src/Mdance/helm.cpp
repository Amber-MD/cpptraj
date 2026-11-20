#ifdef HAS_EIGEN
#include <stdexcept>

#include "helm.h"
#include "bts.h"
#include "scores.h"

using namespace Cpptraj::Mdance;

Mat Helm::makeDataByRow(Vec a, Vec b){
    //a and b need to be same length
    //a will be csum, and b will be sqsum
    Mat data(2, a.size());
    data.row(0) = a;
    data.row(1) = b;

    return data;
}

vector<HCTree> Helm::trimClusters(){
    /*
        Trims the intial clusters based on the trimVal or trimK

        returns
        ------------
        new cluster tree after trimming (vector<HCTree>)
    */
    vector<pair<double, int>> clusterMsds;
    for(size_t i=0; i<clusterTree.size(); i++){
        Vec cSum = clusterTree.at(i).getRootCSum();
        Vec sqSum = clusterTree.at(i).getRootSQSum();
        int Nik = clusterTree.at(i).getRootNObjects();

        if(Nik < minSamples){
            continue;
        }

        Mat data = makeDataByRow(cSum, sqSum);
        clusterMsds.emplace_back(pair<double, int>(extendedComparison(data, Nik, nAtoms, true, mt), i));
    }

    sort(clusterMsds.begin(), clusterMsds.end());

    //trim the clusters based on the trim_k or trim_val
    vector<HCTree> newClusterTree;
    if(trimK){
        this->trimIncoming = clusterMsds.size()-trimK;
        // compare in double so an empty clusterMsds does not wrap around
        if (trimK >= (double)clusterMsds.size()-1){
            throw std::runtime_error("trimK is too large!");
        }
        else if(trimK >= clusterMsds.size()/2.0){
            std::cerr<<"trimK is more than 50% of the clusters. This may lead to poor clustering"<<std::endl;
        }
        for(int i=0; i<(double)clusterMsds.size()-trimK; i++){
            int index = clusterMsds[i].second;
            if (index < 0 || index >= (int)clusterTree.size()) {
                throw std::runtime_error("Index out of bounds when trimming clusters.");
            }
            newClusterTree.emplace_back(clusterTree[index]);
        }
    }
    else if(trimVal){
        this->trimIncoming = 0;
        for(auto i:clusterMsds){
            if(i.first < trimVal){
                trimIncoming++;
            }
        }

        for(auto i:clusterMsds){
            if(i.first < trimVal){
                newClusterTree.emplace_back(clusterTree[i.second]);
            }
        }
    }

    return newClusterTree;
}

vector<HCTree> Helm::genNewClusters(int ZIdx){
    /*
        Generates new cluster by merging two most similar clusters.

        Parameters
        ------------
        zIdx: the Z index to assign to the new merged cluster. This is used for keeping track of the order of merges and for creating the Z matrix.
    */
    int previousSize = clusterTree.size();
    vector<HCTree> previousClusters = clusterTree;

    // update distance matrix with new cluster distances
    if(clusterDists.rows()==0 && clusterDists.cols()==0){
        genClusterDists(previousClusters);
    }
    else{
        Vec distsToNewCluster = Vec::Constant(previousSize-1, INFINITY);
        for(int i=0; i < previousSize - 1; i++){
            float helmSim = calcHelmSim(previousClusters[i], previousClusters[previousSize - 1]); // assumes new cluster is at the end of clusterTree!
            distsToNewCluster[i] = helmSim;
        }

        //Add new cluster to distance matrix
        //np hstack 
        Mat temp1(clusterDists.rows(), clusterDists.cols() + 1);
        temp1.leftCols(clusterDists.cols()) = clusterDists;
        temp1.col(clusterDists.cols()) = distsToNewCluster;
        
        //np vstack
        Mat temp2(temp1.rows() + 1, temp1.cols());
        temp2.topRows(temp1.rows()) = temp1;
        temp2.row(temp1.rows()).setConstant(INFINITY);
        clusterDists = temp2;
    }   

    //Find the two most similar clusters
    Index minRow, minCol;       //minRow and minCol refer to two different clusters
    
    float mergeDist = clusterDists.minCoeff(&minRow, &minCol);
    // add check for minRow and minCol being same. This can cause UB or crash.
    if (minRow == minCol){
        throw std::runtime_error("Got same cluster idx for the two most similar clusters.");
    }

    // eps termination: if the best merge is not similar enough, discard it
    // before it is recorded anywhere. Bailing out down here used to leave a
    // phantom row for the rejected merge in the Z matrix.
    if (eps != -1 && !(mergeDist < eps)) {
        return vector<HCTree>();
    }

    //Save the new clusters after merging (mergeTree computes the combined
    //cSum/sqSum/nObjects itself)
    vector<HCTree> newClusters;
    newClusters.reserve(previousClusters.size()+1);
    for(Index i=0; i<(Index)previousClusters.size(); i++){
        if(i==minRow || i==minCol){;}
        else{
            newClusters.emplace_back(previousClusters[i]);
        }
    }

    // update zMatrix with new merged cluster info
    int idxMinRow = previousClusters[minRow].getRootZIdx();
    int idxMinCol = previousClusters[minCol].getRootZIdx();
    int nClustsMerged = previousClusters[minRow].getRootClusterIndices().size() + previousClusters[minCol].getRootClusterIndices().size();
    updateZMatrix(idxMinRow, idxMinCol, nClustsMerged);

    // merge clusters and add to newClusters
    previousClusters[minRow].mergeTree(previousClusters[minCol], ZIdx);
    newClusters.emplace_back(previousClusters[minRow]);

    //Remove distances of merged clusters
    Veci clustersToKeep(clusterDists.rows()-2);
    int ind=0;
    for(int i=0; i<clusterDists.rows(); i++){
        if(i!= minRow && i!=minCol){
            clustersToKeep[ind] = i;
            ind++;
        }
        
    }
    Mat clusterDists_temp1 = clusterDists(Eigen::all, clustersToKeep);
    Mat clusterDists_temp2 = clusterDists_temp1(clustersToKeep, Eigen::all);
    clusterDists = clusterDists_temp2;

    return newClusters;
}

void Helm::updateZMatrix(int idxA, int idxB, int mergedClusts){
    /* Update the Z matrix with the new merged cluster info. This is used for keeping track of the order of merges and for creating the Z matrix.

     First two columns are the indices of the merged clusters, third column is the distance between the merged clusters, and fourth column is the number of clusters merged to create the new cluster. 

     In our case, distance between merged clusters corresponds to the "rank" of the merge. I.e. the first merge will have distance 1, the second merge will have distance 2, and so on.*/
    if (zMatrix.rows() == 0 && zMatrix.cols() == 0){
        // first merge, initialize Z matrix
        zMatrix = Mat::Zero(1,4);
        zMatrix(0, 0) = idxA;
        zMatrix(0, 1) = idxB;
        zMatrix(0, 2) = 1;
        zMatrix(0, 3) = mergedClusts;
    }
    else{
        Vec zMatrixRow(4);
        int distZMatrix = zMatrix.rows() + 1;
        zMatrixRow << idxA, idxB, distZMatrix, mergedClusts;
        zMatrix.conservativeResize(zMatrix.rows()+1, zMatrix.cols());
        zMatrix.row(zMatrix.rows()-1) = zMatrixRow;
    }
}

float Helm::calcHelmSim(HCTree& firstTree, HCTree& secondTree){
    /*
        Calculates the similarity between two clusters. These are the roots of the two trees that are passed. 

        Parameters
        ------------
        firstTree: HCTree object containing info about first cluster at root
        secondTree: HCTree object containing info about second cluster at root
    */

    Vec cSumA = firstTree.getRootCSum();
    Vec cSumB = secondTree.getRootCSum();
    Vec sqSumA = firstTree.getRootSQSum();
    Vec sqSumB = secondTree.getRootSQSum();
    int nA = firstTree.getRootNObjects();
    int nB = secondTree.getRootNObjects();

    Mat dataA = makeDataByRow(cSumA, sqSumA);
    double simA = extendedComparison(dataA, nA, nAtoms, true, mt);
    Mat dataB = makeDataByRow(cSumB, sqSumB);
    double simB = extendedComparison(dataB, nB, nAtoms, true, mt);

    Vec cSum;
    Vec sqSum;
    cSum = cSumA + cSumB;
    sqSum = sqSumA + sqSumB;

    int n = nA + nB;
    Mat data = makeDataByRow(cSum, sqSum);
    double sim = extendedComparison(data, n, nAtoms, true, mt);

    //Different merging schemes for determining which clusters to merge
    float helmSim = 0;
    if (mergeScheme == MD::MergeScheme::Intra){
        helmSim = sim;
    }
    else if(mergeScheme == MD::MergeScheme::Inter){
        helmSim = ((sim*pow(n,2)) - (simA*pow(nA,2)) - (simB*pow(nB,2)))/(nA*nB);
    }
    else if(mergeScheme == MD::MergeScheme::Half){
        helmSim = sim - ((simA + simB)/2);
    }

    return helmSim;
}

void Helm::genClusterDists(vector<HCTree>& previousClusters){
    /*
        Generates pairwise similairty matrix for initial clusters

        Parameters
        -------------
        previousClusters: contains info about clusters in kth iteration

        Returns
        -------------
        pairwise similarity matrix (void)
    */

    clusterDists = Mat::Zero(previousClusters.size(), previousClusters.size()) + INFINITY;
    for(size_t i=0; i<previousClusters.size(); i++){
        for(size_t j=i+1; j<previousClusters.size(); j++){
            float helmSim = calcHelmSim(previousClusters[i], previousClusters[j]);
            clusterDists(i, j) = helmSim;
        }
    }
}


Helm::Helm(vector<HCTree> clusterTree, int nAtoms, MD::Metric mt, 
        MD::MergeScheme mergeScheme, int nClusters, float eps, 
        bool trimStart, float minSamples,
        float trimVal, float trimK,
        bool savePairwiseSum,
        string inputTop, string inputTraj){
            /* Constructor for Helm class. Initializes the class variables and checks for end conditions.
        Parameters
        ----------
        clusterTree: 
            vector of HCTree objects containing info about initial clusters at root
        nAtoms: 
            number of atoms in the system. This is used for calculating the similarity between clusters
        mt: 
            metric to use for calculating similarity between clusters. This is used for calculating the similarity between clusters
        mergeScheme: 
            the scheme to use for determining which clusters to merge. This is used for calculating the similarity between clusters
        nClusters: 
            the number of clusters to return.
        eps:
            epsilon MSD value to terminate clustering process.
        trimStart:
            whether to trim the initial clusters before starting the clustering process. This is used for improving the clustering results and reducing the runtime. If true, then either trimVal or trimK must be provided. This is used for improving the clustering results and reducing the runtime.
        minSamples:
            the minimum number of samples a cluster must have to be considered for trimming. This is used for improving the clustering results and reducing the runtime.
        trimVal:
            the MSD value to use for trimming the initial clusters. This is used for improving the clustering results and reducing the runtime. 
        trimK:
            the number of clusters to keep after trimming the initial clusters. This is used for improving the clustering results and reducing the runtime. 
        savePairwiseSum:
            wether to save pairwise similarity matrix. Note: this is currently unused. 
        inputTop:
            topology file of the MD system
        inputTraj:
            trajectory file of the MD system
            */
    this->clusterTree = clusterTree;
    this->mt = mt;
    this->nAtoms = nAtoms;
    this->mergeScheme = mergeScheme;
    this->nClusters = nClusters;
    this->eps = eps;
    this->trimStart = trimStart;
    this->minSamples = minSamples;
    this->trimVal = trimVal;
    this->trimK = trimK;
    this->savePairwiseSum = savePairwiseSum;
    this->inputTop = inputTop;
    this->inputTraj = inputTraj;

    this->totalIncoming = clusterTree.size();;
    totalSum = 0;
    for(auto tree:clusterTree){
        this->totalSum += tree.getRootNObjects();
    }

    //check end conditions
    if((this->nClusters==0 && this->eps==-1) || 
        (this->nClusters > 0 && this->eps!=-1)){
            throw std::invalid_argument("You must provide either nClusters or eps, but not both.");
        }

    if(this->trimStart && !(this->trimVal || this->trimK)){
        throw std::invalid_argument("If trimStart is true, then either trimVal or trimK must be provided.");
    }
    if(this->trimVal && this->trimK){
        throw std::invalid_argument("You can only provided either trimVal or trimK, but not both.");
    }

    if(this->minSamples < 0){
        throw std::invalid_argument("minSamples must be greater than 0.");
    }
    else if(0 < this->minSamples && this->minSamples < 1){
        this->minSamples = int(this->minSamples * this->totalSum);
    }
    else if(this->minSamples >= 1){
        this->minSamples = int(this->minSamples);
    }

}

vector<HCTree> Helm::run(){
    /*
        Performs HELM clustering of initial clusters

        Returns
        -----------
        tree of clusters (vector<HCTree>) after performing HELM clustering
    */

    if(nClusters == 0){
        nClusters = 1;
    }

    if(trimStart){
        clusterTree = trimClusters();
    }

    //perform clustering
    int n = clusterTree.size();
    int currentZInd = clusterTree.size();
    if (nClusters == n) {
        std::cerr << "Number of clusters is already equal to nClusters. No clustering will be performed." << std::endl;
        return clusterTree;
    }
    while(n>1){
        vector<HCTree> newClusters = genNewClusters(currentZInd);
        if (!newClusters.empty()){ 
            // if new clusters are empty, clustering is stopped. In order to return the results, we cannot return the empty clusters, so we return the previous clusters.
            clusterTree = newClusters;
        }
        currentZInd +=1; // increment Z index for new merged cluster
        //termination conditions
        if(n==(nClusters+1) || newClusters.empty()){
            break;
        }
        n-=1;
    }
    return clusterTree;
}

pair<double, double> Helm::computeScores(vector<HCTree> clusters, Mat data){
    /*
        Computes Calinksi-Harabasz and Davies-Bouldin scores of clusters
        using random labeling

        Parameters
        -------------
        clusters: vector of HCTree objects containing info about clusters at root
        data: the data matrix used for clustering. This is used for calculating the scores. The row order has to be the same as the cluster and frame order in clusters. E.g. given the following clusters = [cluster1: frames 0, 2, 4; cluster2: frames 1, 3, 5], the data matrix should have the rows in the order of [frame0, frame2, frame4, frame1, frame3, frame5].

        Returns
        -------------
        pair: first element is Calinkski, second element Davies-Bouldin
    */

    vector<int> frameLabel;
    int count=0;

    for(auto c:clusters){
        //clusterIndices.emplace_back(c.getIndices());
        int Nik=c.getRootNObjects();
        frameLabel.insert(frameLabel.end(), Nik, count);
        count++;
    }
    
    set<int> frameLabelSet;
    for(int i:frameLabel){
        frameLabelSet.insert(i);
    }
    if(frameLabelSet.size()==1){// everything is in one cluster
        return pair<double, double>{-1,-1};
    }
    else{
        Veci labels = Eigen::Map<Veci>(frameLabel.data(), frameLabel.size());
        double chScore = calinskiHarabaszScore(data, labels);
        double dbScore = daviesBouldinScore(data, labels);

        return pair<double, double>{chScore, dbScore};
    }
}

Mat Helm::getZMatrix(){
    /* 
        getter for Z matrix. This is used for keeping track of the order of merges and for creating the Z matrix.

        Returns
        -------------
        Z matrix (Mat)
    */
    if (zMatrix.rows() == 0 && zMatrix.cols() == 0){
        std::cerr<<"Z matrix is empty. No clusters were merged."<<std::endl;
    }
    return zMatrix;
}
#endif /* HAS_EIGEN */
