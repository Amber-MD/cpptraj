#include "hc_utils.h"
#ifdef HAS_EIGEN
using namespace Cpptraj::Mdance;
/**
 * @brief Constructor for the HCTreeNode class.
 * @param clusterIdx list of cluster identifiers corresponding to the initial clustering (before running HELM). Each value represents an original cluster label contained in this tree node. These are cluster labels, not frame indices.
 * @param cSumIn A 1D Eigen Array representing the columnwise sum of the data for the molecules in this node. 
 * @param sqsum A 1D Eigen Array representing the columnwise sum of squares of the data for the molecules in this node.
 * @param nObjects An integer representing the number of molecules/frames in this node.
 * @param z_ind An integer representing the z_index for this node. This is the index of the node in the linkage matrix Z.
 */
HCTreeNode::HCTreeNode(std::list<int> clusterIdx, Vec cSumIn, Vec sqsum, int nObjects, int z_ind){
    setCSum(cSumIn);
    setSQSum(sqsum);
    setClusterIndices(clusterIdx);
    setNObjects(nObjects);
    setLeftPtr(nullptr);
    setRightPtr(nullptr);
    setZIdx(z_ind);
};

/**
 * @brief Get the columwise sum of data for this node.
 * @return Vec The cSum for this node.
 */
Vec HCTreeNode::getCSum(){
    return cSum;
}


void HCTreeNode::setCSum(Vec cSumIn){
    cSum = cSumIn;
}

Vec HCTreeNode::getSQSum(){
    return sqSum;
}

void HCTreeNode::setSQSum(Vec sqSumIn){
    sqSum = sqSumIn;
}

/**
 * @brief Get the z_index for this node.
 * @return int The z_index for this node.
 */
int HCTreeNode::getZIdx(){
    return zIndex;
}

/**
 * @brief Set the z_index for this node.
 * @param z_ind An integer representing the z_index for this node.
 */
void HCTreeNode::setZIdx(int z_ind){
    zIndex = z_ind;
}

/**
 * @brief Get the left pointer for this node.
 * @return HCTreeNode* A pointer to the left child node.
 */
std::shared_ptr<HCTreeNode> HCTreeNode::getLeftPtr(){
    return leftPtr;
}

/**
 * @brief Set the left pointer for this node.
 * @param input_left A pointer to the left child node.
 */
void HCTreeNode::setLeftPtr(std::shared_ptr<HCTreeNode> input_left){
    leftPtr = input_left;
}

/**
 * @brief Get the right pointer for this node.
 * @return std::shared_ptr<HCTreeNode> A pointer to the right child node.
 */
std::shared_ptr<HCTreeNode> HCTreeNode::getRightPtr(){
    return rightPtr;
}

/**
 * @brief Set the right pointer for this node.
 * @param input_right A pointer to the right child node.
 */
void HCTreeNode::setRightPtr(std::shared_ptr<HCTreeNode> input_right){
    rightPtr = input_right;
}

/**
 * @brief Get the number of molecules/frames in this node.
 * @return int The number of molecules/frames in this node.
 */
int HCTreeNode::getNObjects(){
    return nObjects;
}

/**
 * @brief Set the number of objects in this node.
 * @param n_mol An integer representing the number of objects in this node.
 */
void HCTreeNode::setNObjects(int n_mol){
    nObjects = n_mol;
}

/**
 * @brief Get the list of indices for the clusters in this node.
 * @return std::list<int> A list of integers representing the clusters in this node.
 */
std::list<int> HCTreeNode::getClusterIndices(){
    return clusterIndices;
}

/**
 * @brief Set the list of indices for the clusters in this node.
 * @param idx A list of integers representing the indices of the clusters in this node.
 */
void HCTreeNode::setClusterIndices(std::list<int> idx){
    clusterIndices=idx;
}



/**
 * @brief Constructor for the HCTree class.
 * 
 * @details
 * The HCTree class represents a hierarchical clustering tree. It contains a pointer to the root node
 * of the tree, which is an instance of the HCTreeNode class.
 * 
 */
HCTree::HCTree(){
    rootPtr = nullptr;
}

/**
 * @brief Get the root node of the hierarchical clustering tree.
 *
 * @return std::shared_ptr<HCTreeNode> A pointer to the root node of the tree.
 */
std::shared_ptr<HCTreeNode> HCTree::getRoot(){
    return rootPtr;
}

/**
 * @brief Set the root node of the hierarchical clustering tree.
 *
 * @param input_root A pointer to the new root node of the tree.
 */
void HCTree::setRoot(std::shared_ptr<HCTreeNode> input_root){
    rootPtr = input_root;
}

/**
 * @brief Get the columnwise sum of data for the root node of the hierarchical clustering tree.
 * 
 * @return Vec The columnwise sum of data for the root node of the tree.
 */
Vec HCTree::getRootCSum(){
    return rootPtr->getCSum();
}

Vec HCTree::getRootSQSum(){
    return rootPtr->getSQSum();
}

/**
 * @brief Get the number of objects (molecules/frames) in the root node of the hierarchical clustering tree.
 * 
 * @return int The number of objects in the root node of the tree.
 */
int HCTree::getRootNObjects(){
    return rootPtr->getNObjects();
}

/**
 * @brief Get the list of indices for the clusters in the root node of the hierarchical clustering tree. These indices are cluster identifiers corresponding to the initial clustering (before running HELM). Each value represents an original cluster label contained in this tree node. These are cluster labels, not frame indices.
 * 
 * @return std::list<int> A list of integers representing the indices of the clusters in the root node of the tree.
 */
std::list<int> HCTree::getRootClusterIndices(){
    return rootPtr->getClusterIndices();
}

/**
 * @brief Get the z_index for the root node of the hierarchical clustering tree.
 * 
 * @return int The z_index for the root node of the tree.
 */
int HCTree::getRootZIdx(){ //z_index
    return rootPtr->getZIdx();
}

/**
 * @brief Insert a new root node into the hierarchical clustering tree.
 * 
 * @param idx A list of integers representing the indices of the clusters in the new root node. These indices are cluster identifiers corresponding to the initial clustering (before running HELM). Each value represents an original cluster label contained in this tree node. These are cluster labels, not frame indices.
 * @param cSum A 1D Eigen Array representing the columnwise sum of the data for the molecules in the new root node.
 * @param sqsum A 1D Eigen Array representing the columnwise sum of squares of the data for the molecules in the new root node.
 * @param nObjects An integer representing the number of molecules/frames in the new root node.
 * @param z_ind An integer representing the z_index for the new root node.
 */
void HCTree::insertRoot(std::list<int> idx, Vec cSum, Vec sqsum, int nObjects, int z_ind){
    rootPtr = std::make_shared<HCTreeNode>(idx, cSum, sqsum, nObjects, z_ind);
}

/**
 * @brief Combine two hierarchical clustering trees into one.
 * 
 * @details
 * This method merges another tree into the current one according to the hierarchical clustering algorithm. It calculates
 * the new columnwise sum and indices for the combined tree, creates a new root node, and sets the
 * left and right pointers accordingly.
 * 
 * @param other_tree The other HCTree to be combined with this tree.
 * @param new_z_ind An integer representing the z_index for the new root node of the combined tree.
 */
void HCTree::mergeTree(HCTree other_tree, int new_z_ind){
    // This method combines two trees according, which is needed for hierarch. clust.
    
    // calculate new colsum for fingerprints, indices
    Vec csum_new = other_tree.getRootCSum() + getRootCSum();
    Vec sqsum_new = other_tree.getRootSQSum() + getRootSQSum();
    std::list<int> idx = other_tree.getRootClusterIndices();
    idx.splice(idx.end(), getRootClusterIndices());
    idx.sort();
    int nObjectsNew = other_tree.getRootNObjects() + getRootNObjects();
    //create new root node
    std::shared_ptr<HCTreeNode> new_root = std::make_shared<HCTreeNode>(idx, csum_new, sqsum_new, nObjectsNew, new_z_ind);

    // set left and right pointers so that z_left < z_right
    int z_other_tree = other_tree.getRootZIdx();
    int z_this_tree = getRootZIdx();
    if (z_this_tree < z_other_tree){
        new_root->setLeftPtr(getRoot());
        new_root->setRightPtr(other_tree.getRoot());
    }
    else{
        new_root->setLeftPtr(other_tree.getRoot());
        new_root->setRightPtr(getRoot());
    }
    // update root
    setRoot(new_root);
}
#endif /* HAS_EIGEN */
