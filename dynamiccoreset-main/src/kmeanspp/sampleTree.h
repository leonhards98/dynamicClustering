#include "../dynamic_coreset/coreset_point.h"
#include <memory>
#include <queue>
namespace fastKMeans
{
    /**
     * @brief A node of the sample tree. Contains one point and the weight of the edge to its parent
     * @author
     * @since Tue May 02 2023
     */
    struct treeNode
    {
        /**
         * The coordinates of the point stored in this node. Empty if it is not a leaf
         */
        blaze::DynamicVector<double> coordinates;

        /**
         * Weight of the point stored in this node
         */
        double pointWeight;

        /**
         * Index of the point stored in this node
         */
        int index;

        /**
         * Weight of the edge from this node to its parent
         */
        double KMeansWeight{0};

        /**
         * Shared pointer to the left child of this node
         */
        std::shared_ptr<treeNode> lChild;

        /**
         * Shared pointer to the right child of this node
         */
        std::shared_ptr<treeNode> rChild;

        /**
         * Pointer to the parent of this node
         */
        treeNode *parent;

        /**
         * @brief Initializes the weights of all edges (KMeansWeight). Only works if all points have a weight of 1
         * @param nrPoints Total number of points in the tree
         * @param M Constant used to caluclate the initial weight
         * @return (void)
         */
        void resetWeightWithoutPointWeight(int nrPoints, double M)
        {
            assert(nrPoints >= 1);
            this->KMeansWeight = nrPoints * M;

            if (nrPoints == 1)
                return;

            // At each point, this function is called again for the children of this node.
            // This procedure calculates those
            double rest = nrPoints - pow(2, blaze::floor(blaze::log2(nrPoints)));
            if (rest == 0)
            {
                lChild->resetWeightWithoutPointWeight(nrPoints / 2, M);
                rChild->resetWeightWithoutPointWeight(nrPoints / 2, M);
                return;
            }

            double halve = pow(2, blaze::floor(blaze::log2(nrPoints))) / 2;

            if (rest <= halve)
            {
                lChild->resetWeightWithoutPointWeight(halve + rest, M);
                rChild->resetWeightWithoutPointWeight(nrPoints - (halve + rest), M);
            }
            else
            {
                lChild->resetWeightWithoutPointWeight(halve * 2, M);
                rChild->resetWeightWithoutPointWeight(nrPoints - (halve * 2), M);
            }
            return;
        }
    };

    /**
     * @brief Represets a sample tree. This is a balanced binary tree, where each leaf corresponds to a data point.
     * Used to determine which point to chose as the next center. Each edge weight is the probability which which its
     * path (starting from  the root) should be taken when a new leaf is chose.
     */
    class sampleTree
    {
    private:
        /**
         * The root of the sample tree
         */
        std::shared_ptr<treeNode> root;

        /**
         * A queue containing leafs. When a new point needs to be inserted, it's new leaf should be contructed as a child of the
         * next leaf in this queue
         */
        std::queue<std::shared_ptr<treeNode>> insertionLeafs;

    public:
        /**
         * @brief Constructor
         */
        sampleTree()
        {
            this->root = std::make_shared<treeNode>();
            root->parent = nullptr;
            root->lChild = nullptr;
            root->rChild = nullptr;
            root->KMeansWeight = 0;
            root->pointWeight = 0;
            this->insertionLeafs.push(root);
        }

        /**
         * @brief Initializes the tree with a number of points. If not used dynamically, all points should be inserted this way.
         * Also initializes all weights of the tree. The index of each point is determined by its row in the data matrix
         * @param data A matrix containing all data points. Each row corresponds to the coordinates of one point.
         * @param weights A vector containing the weights of all points. The index in this vector should be the same as the
         * row in which the point that has the weight is stored
         * @param M A constant for determining the initial weights.
         * @return A vector containing the tree nodes of all points
         */
        std::vector<std::shared_ptr<treeNode>> initTree(const blaze::DynamicMatrix<double> &data, const blaze::DynamicVector<double> &weights, double M);

        /**
         * @brief Search for a point using its index. Uses a breadth first search
         * @param index Index of the point to search for
         * @return A pointer to the tree node containing the point with the given index
         */
        std::shared_ptr<treeNode> searchPoint(int index);

        /**
         * @brief Returns a vector containing pointer to all leafs (all nodes that store a point)
         * @return A vector containing pointers to all leafs
         */
        std::vector<std::shared_ptr<treeNode>> getAllLeafs();

        /**
         * @brief Searches for the leftmost leaf in the highest level that contains leafs
         * @return A pointer to the leftmost leaf in the highest level that contains leafs
         */
        std::shared_ptr<treeNode> searchLeftmostLeaf();

        /**
         * @brief Searches for the rightmost leaf in the lowest level
         * @return A pointer to the rightmost leaf in the lowest level
         */
        std::shared_ptr<treeNode> searchRightmostLeaf();

        /**
         * @brief Swaps the position of the two given leafs in the tree
         * @param leaf1
         * @param leaf2
         * @return
         */
        bool swapLeafs(std::shared_ptr<treeNode> leaf1, std::shared_ptr<treeNode> leaf2);

        /**
         * @brief Returns the depth of the tree
         * @return Depth of the tree
         */
        int getTreeDepth();

        const std::shared_ptr<treeNode> &getRoot();

        /**
         * @brief Updates all edge weights from a leaf where the weight canges on the path over its anchestors to the root.
         * Also updates the weight of the leaf
         * @param leaf Leaf for which the weight (weight of the edge from this leaf to its parent) changed
         * @param newDist New weight
         * @param oldDist Old weight
         * @return (void)
         */
        void updateWeights(treeNode *leaf, double newDist, double oldDist);

        /**
         * @brief Same as "update weights", but assumes the weight of the leaf has already been modified. Updates the edge weights
         * on the path from the parent of the given leaf node to the root.
         * @param node A leaf whose weight has changed
         * @param newDist New weight
         * @param oldDist Old weight
         * @return (void)
         */
        void updateInternalWeights(treeNode *node, double newDist, double oldDist);

        /**
         * @brief Initializes the weights of all edges (KMeansWeight). Also works if the points have weight that are not 1
         * @param nrPoints Total number of points in the tree
         * @param M Constant used to caluclate the initial weight
         * @return (void)
         */
        void resetWeightsWithPointWeights(double M);

        /**
         * @brief Inserts a point into the sample tree
         * @param cooridinates A vector containing the coordinates of the point
         * @param pointWeight The weight of the point
         * @param index Index of the point. Must be unique
         * @return  A shared pointer to the leaf that stores the point
         */
        std::shared_ptr<treeNode> insertPoint(blaze::DynamicVector<double> cooridinates, double pointWeight, double index);

        /**
         * @brief Removes a point from the sample tree
         * @param point a shared pointer to the leaf that contains the point that should be removed
         * @return (void)
         */
        void removePoint(std::shared_ptr<treeNode> point);
    };

}
