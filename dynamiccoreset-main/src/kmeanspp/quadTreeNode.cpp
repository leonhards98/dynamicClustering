#include "quadTreeNode.h"
namespace fastKMeans
{
        void QuadtreeNode::setAllBorders(std::shared_ptr<treeNode> point, double border)
        {

            boost::random::mt19937 randEngine;
            boost::random::uniform_real_distribution<> r(0, 1);

            double rand{0};
            while (rand == 0)
                rand = r(randEngine) * (border / 2);

            if (r(randEngine) > 0.5)
                rand *= -1;

            assert(rand >= -border / 2);
            assert(rand <= border / 2);

            for (int i = 0; i < point->coordinates.size(); i++)
            {
                borders.at(i).first = point->coordinates.at(i) - border - border / 2 + rand;
                borders.at(i).second = point->coordinates.at(i) + border + border / 2 + rand;
            }
        }

	bool QuadtreeNode::insertPoint(std::shared_ptr<treeNode> point)
        {
            this->pointInSubtree.push_back(point);

            /**
             * A vector ob bools. Is set to 1 for an index, if the i'th dimesion of the point lies in the upper halve of the
             * range in space this node is responsible for.
             */
            std::vector<bool> childrenBits(std::pow(2, point->coordinates.size()), false);

            int dimension = point->coordinates.size();

            /*
            std::cout << point->coordinates << "\n";
            std::cout << borders.at(0).first;
            std::cout << borders.at(0).second << "\n";;
            std::cout << borders.at(1).first;
            std::cout << borders.at(1).second << "\n";
            double r = point->coordinates.at(0);
            double r2 = point->coordinates.at(1);
            */
            for (int i = 0; i < dimension; i++)
            {
                assert(point->coordinates.at(i) >= borders.at(i).first);
                assert(point->coordinates.at(i) <= borders.at(i).second);
            }

            /**
             * The index of the child in whose part of space the new point lies.
             * If the i'th dimenstion of the new point lies in the upper halve of the range tis node is
             * responsible for, 2^i is added to childNr
             */
            int childNr{0};

            // Split the space the current node is responsible for in halve for each dimension and
            // check in which halve the new point lies. update childNr accordingly
            for (int i = 0; i < dimension; i++)
            {
                double cutOff = (borders.at(i).first + borders.at(i).second) / 2;

                if (point->coordinates.at(i) > cutOff)
                {
                    childNr += pow(2, i);
                    childrenBits.at(i) = true;
                }
            }

            // If there is no child present for the range in which the new point lies, and the current node is not
            // a leaf (doesn't store a point), create a new child and insert the point.
            // If the current Node is a leaf (contains a point), first create a new child for the point in the current
            // node (call insertPoint()) and then call insertPoint for the original point that needs to be inserted
            if (children[childNr] == nullptr)
            {
                if (this->point != nullptr)
                {
                    int i{0};
                    for (int i = 0; i < dimension; i++)
                        if (point->coordinates.at(i) == this->point->coordinates.at(i))
                            i++;
                    assert(i < dimension);
                    std::shared_ptr<treeNode> p = this->point;
                    this->point.reset();
                    this->insertPoint(p);
                    this->insertPoint(point);
                    return 0;
                }
                children[childNr] = std::make_shared<QuadtreeNode>(dimension, this->weightOfEdgeToParent / 2);
                children[childNr]->point = point;
                children[childNr]->parent = this;

                children[childNr]->isMarked = this->isMarked;

                for (int i = 0; i < dimension; i++)
                {
                    double cutOff = (borders.at(i).first + borders.at(i).second) / 2;

                    if (childrenBits.at(i) == true)
                    {
                        children[childNr]->borders.at(i).first = cutOff;
                        children[childNr]->borders.at(i).second = this->borders.at(i).second;
                    }
                    else
                    {
                        children[childNr]->borders.at(i).first = this->borders.at(i).first;
                        children[childNr]->borders.at(i).second = cutOff;
                    }
                }
            }
            // If there is a child for the range in which the new point lies, call insertPoint for this child
            else
            {
                children[childNr]->insertPoint(point);
            }

            assert(this->point == nullptr);

            return 0;
        }

        bool QuadtreeNode::removePoint(std::shared_ptr<treeNode> point)
        {

            assert(this->point == nullptr);
            auto point_it = std::find(this->pointInSubtree.begin(), this->pointInSubtree.end(), point);
            assert(point_it != this->pointInSubtree.end());
            this->pointInSubtree.erase(point_it);

            std::vector<bool> childrenBits(std::pow(2, point->coordinates.size()), false);

            int dimension = point->coordinates.size();

            /*
            std::cout << point->coordinates << "\n";
            std::cout << borders.at(0).first;
            std::cout << borders.at(0).second << "\n";;
            std::cout << borders.at(1).first;
            std::cout << borders.at(1).second << "\n";
            double r = point->coordinates.at(0);
            double r2 = point->coordinates.at(1);
            */
            for (int i = 0; i < dimension; i++)
            {
                assert(point->coordinates.at(i) >= borders.at(i).first);
                assert(point->coordinates.at(i) <= borders.at(i).second);
            }

            /**
             * The index of the child in whose part of space the new point lies.
             * If the i'th dimenstion of the new point lies in the upper halve of the range tis node is
             * responsible for, 2^i is added to childNr
             */
            int childNr{0};

            // Get the index of the child responsible for the point that should be removed
            for (int i = 0; i < dimension; i++)
            {
                double cutOff = (borders.at(i).first + borders.at(i).second) / 2;

                if (point->coordinates.at(i) > cutOff)
                {
                    childNr += pow(2, i);
                    childrenBits.at(i) = true;
                }
            }
            assert(children[childNr] != nullptr);

            // If the point is stored directly in the responsible child, remove the point and then the child.
            // If only one other point remains in the subtree of this node, move the point to this node
            // and delete the remaining child of this node (this node becomes a leaf).
            if (children[childNr]->point == point)
            {
                children[childNr]->point.reset();
                assert(children[childNr]->pointInSubtree.size() == 0);
                children[childNr].reset();
                if (this->pointInSubtree.size() == 1)
                {
                    this->pointInSubtree.clear();
                    for (auto c = children.begin(); c != children.end(); c++)
                    {
                        if (*c == nullptr)
                            continue;

                        assert((*c)->point != nullptr);
                        this->point = (*c)->point;
                        (*c)->point.reset();
                        (*c).reset();
                        break;
                    }
                }
                return 0;
            }

            //If the point is not directly stored in the responsible child, call removePoint on that child
            children[childNr]->removePoint(point);


            // If only one other point remains in the subtree of this node, move the point to this node
            // and delete the remaining child of this node (this node becomes a leaf).
            if (this->pointInSubtree.size() == 1)
            {
                assert(this->point == nullptr);
                this->pointInSubtree.clear();
                for (auto c = children.begin(); c != children.end(); c++)
                {
                    if (*c == nullptr)
                        continue;

                    assert((*c)->point != nullptr);
                    this->point = (*c)->point;
                    (*c)->point.reset();
                    (*c).reset();
                    break;
                }
            }

            return 0;
        }

        double QuadtreeNode::calcMinDist()
        {
            double dist{0};
            QuadtreeNode *currentNode = this;

            while (currentNode->parent != nullptr)
            {
                dist += currentNode->weightOfEdgeToParent;
                currentNode = currentNode->parent;
                if (currentNode->isMarked)
                    break;
            }
            return dist;
        }
}
