#include "sampleTree.h"
#include <assert.h>
#include <boost/random.hpp>

namespace fastKMeans
{
    /**
     * @brief A node in the quad tree, responsible for a part of the space. If it is a leaf, containts exactly one point
     */
    class QuadtreeNode
    {
    public:
        /**
         * Pointer to the sample tree node that contains the point which is stored in this node
         */
        std::shared_ptr<treeNode> point;

        /**
         * Vector containing pointers to all the sample tree nodes that contain points
         * which are in leafs that are in the subtree of this node
         */
        std::vector<std::shared_ptr<treeNode>> pointInSubtree;

        /**
         * Vector containing shared pointers to the children of this node
         */
        std::vector<std::shared_ptr<QuadtreeNode>> children;

        /**
         * Pointer to the parent of this node
         */
        QuadtreeNode *parent;

        /**
         * A vector containing the left and right border of the part of space this node is representing for each dimenstion.
         * Each element in the vector deals with one dimension
         */
        std::vector<std::pair<double, double>> borders;

        /**
         * Weight of the edge to the parent node. Used to approximate the distance between points
         */
        double weightOfEdgeToParent;

        /**
         * Bool if the node is marked
         */
        bool isMarked{false};

        /**
         * @brief Constructor
         * @param dimension The dimenstion of the points that are stored in the tree
         * @param weightOfEdgeToParent weight of the edge to the parent
         */
        QuadtreeNode(int dimension, double weightOfEdgeToParent)
        {
            children.resize(std::pow(2, dimension));
            for (int i = 0; i < children.size(); i++)
                children.at(i) = nullptr;

            borders.resize(dimension * 2);
            parent = nullptr;
            this->weightOfEdgeToParent = weightOfEdgeToParent;
        }

        /**
         * @brief Destructor
         */
        ~QuadtreeNode()
        {
            point.reset();
        }

        /**
         * @brief Set the initial borders of the root of the tree with a random offset to make the quad
         * trees different. uses a single starting point to center the random shift around
         * @param point Initial point
         * @param border Original offest from the initial point to catach all remaining points in the space
         * @return (void)
         */
        void setAllBorders(std::shared_ptr<treeNode> point, double border);

        /**
         * @brief Inserts a Point into the tree. Checks if the node can be inserted into this node. If not, calls this function of
         * the child corresponding to the part of space the point is part of or creates a new node and inserts the point there.
         * Also updates the vector "pointsInSubtree"
         * @param point pointer to the sample tree leaf that contains the point that needs to be inserted
         * @return
         */
        bool insertPoint(std::shared_ptr<treeNode> point);

        /**
         * @brief Remove a point from the quad tree
         * @param point Pointer to the sample tree node that stores the point to be removed
         * @return
         */
        bool removePoint(std::shared_ptr<treeNode> point);

        
        /**
         * @brief Calculate the approximation of this tree for the point stored in this node to the closest center
         * @return The distance from the point in this node to the closest center
         */
        double calcMinDist();
    };
}
