#include "quadTree.h"


namespace fastKMeans
{
    /**
     * @brief Indert a point into the quad tree
     * @param point A pointer to the sample tree node that stores the point that should be inserted
     * @return (void)
     */
    void quadTree::insertPoint(std::shared_ptr<treeNode> point)
    {
        root->insertPoint(point);
    }

    /**
     * @brief Remove a point from the quad tree
     * @param point Pointer to the sample tree node storing the poin to remove
     * @return (void)
     */
    void quadTree::removePoint(std::shared_ptr<treeNode> point)
    {
        root->removePoint(point);
    }

    /**
     * @brief Find the leaf containing a given point
     * @param point pointer to the sample tree node containing the point to search for
     * @return Pointer to the quad tree node containing the point
     */
    std::shared_ptr<QuadtreeNode> quadTree::findLeaf(std::shared_ptr<treeNode> point)
    {
        int dimension = point->coordinates.size();
        std::shared_ptr<QuadtreeNode> currentNode = this->root;
        int childNr{0};

        while (currentNode->point.get() != point.get())
        {
            for (int i = 0; i < dimension; i++)
            {
                assert(point->coordinates.at(i) >= currentNode->borders.at(i).first);
                assert(point->coordinates.at(i) <= currentNode->borders.at(i).second);
            }

            childNr = 0;
            for (int i = 0; i < dimension; i++)
            {
                double cutOff = (currentNode->borders.at(i).first + currentNode->borders.at(i).second) / 2;

                if (point->coordinates.at(i) > cutOff)
                    childNr += pow(2, i);
            }
            assert(currentNode->children[childNr] != nullptr);
            currentNode = currentNode->children[childNr];
        }

        return currentNode;
    }

    /**
     * @brief Marks all nodes on the path from the given tree node as marked, until the root or the first node
     * that is already marked is found.
     * @param point pointer to the quad tree node from which the marking should start
     * @return The first point that was not marked (either the root or the first point that was already marked)
     */
    QuadtreeNode *quadTree::markPath(QuadtreeNode *point)
    {
        while (point->parent->isMarked == false)
        {
            point->isMarked = true;
            point = point->parent;
            if (point->parent == nullptr)
                break;
        }

        return point;
    }

    /**
     * @brief Calculates the approximation of the tree for the distance between two points.
     * @param x Quad tree node containing the first point
     * @param y Sample tree node containing the second point
     * @param ancestor A common ancestor of the two points. At worst, the root of the tree
     * @return The calculated distance
     */
    double quadTree::calcTreeDist(QuadtreeNode *x, treeNode *y, QuadtreeNode *ancestor)
    {
        //!!! Could also find the weight by moving up from x until y is in the subtree

        if (x->point.get() == y)
            return 0;
        int dimension = y->coordinates.size();
        int childNrX{0};
        int childNrY{0};
        double distance{0};

        // Starting from the given ancestor, move into the tree, as long as both points are part of the same child.
        //  This means "ancestor" points to the closest common ancestor of both points after this loop
        while (childNrX == childNrY)
        {
            childNrX = 0;
            childNrY = 0;
            for (int i = 0; i < dimension; i++)
            {
                double cutOff = (ancestor->borders.at(i).first + ancestor->borders.at(i).second) / 2;

                if (x->point->coordinates.at(i) > cutOff)
                    childNrX += pow(2, i);
                if (y->coordinates.at(i) > cutOff)
                    childNrY += pow(2, i);
            }
            assert(ancestor->children[childNrX] != nullptr);
            assert(ancestor->children[childNrY] != nullptr);

            if (childNrX == childNrY)
                ancestor = ancestor->children[childNrX].get();
            else
                break;
        }

        // Sum the weights of the edges from one point to the closest common ancestor.
        // This is the approximation of the distance
        while (x->parent != ancestor)
        {
            distance += x->weightOfEdgeToParent;
            x = x->parent;
        }
        distance += x->weightOfEdgeToParent;

        return distance;
    }
}
