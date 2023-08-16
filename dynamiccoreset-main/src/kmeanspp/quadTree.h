#include "quadTreeNode.h"

namespace fastKMeans
{
    /**
     * @brief A quad tree used to approximate the distace between ever point and the closest center
     */
    class quadTree
    {
    private:
        
        /**
         * An approcimation for the maximal distance between any two points in the input data
        */
        double maxDist;

        /**
         * The root of the tree
        */
        std::shared_ptr<QuadtreeNode> root;

    public:
        /**
         * @brief Constructor
         * @param initialPoint The point used to initialize the quad tree
         * @param dimension The number of dimenstions of the points in the input data
         * @param maxDist At least maximal distance between any two points in the input data
         */
        quadTree(std::shared_ptr<treeNode> initialPoint, int dimension, double maxDist)
        {
            this->maxDist = maxDist;
            root = std::make_shared<QuadtreeNode>(dimension, blaze::sqrt(dimension) * 2 * maxDist);
            root->setAllBorders(initialPoint, maxDist);
        }

        /**
         * @brief Indert a point into the quad tree
         * @param point A pointer to the sample tree node that stores the point that should be inserted
         * @return (void)
         */
        void insertPoint(std::shared_ptr<treeNode> point);

        /**
         * @brief Remove a point from the quad tree
         * @param point Pointer to the sample tree node storing the poin to remove
         * @return (void)
         */
        void removePoint(std::shared_ptr<treeNode> point);

        /**
         * @brief Find the leaf containing a given point
         * @param point pointer to the sample tree node containing the point to search for
         * @return Pointer to the quad tree node containing the point
         */
        std::shared_ptr<QuadtreeNode> findLeaf(std::shared_ptr<treeNode> point);

        /**
         * @brief Marks all nodes on the path from the given tree node as marked, until the root or the first node
         * that is already marked is found.
         * @param point pointer to the quad tree node from which the marking should start
         * @return The first point that was not marked (either the root or the first point that was already marked)
         */
        QuadtreeNode *markPath(QuadtreeNode *point);

        /**
         * @brief Calculates the approximation of the tree for the distance between two points.
         * @param x Quad tree node containing the first point
         * @param y Sample tree node containing the second point
         * @param ancestor A common ancestor of the two points. At worst, the root of the tree
         * @return The calculated distance
         */
        double calcTreeDist(QuadtreeNode *x, treeNode *y, QuadtreeNode *ancestor);
    };

}
