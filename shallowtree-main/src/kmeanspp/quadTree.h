#include "sampleTree.h"
#include <assert.h>
#include <boost/random.hpp>

namespace fastKMeans
{
    
    class QuadtreeNode
    {
        public: 

        std::shared_ptr<treeNode> point;
        std::vector<std::shared_ptr<treeNode>> pointInSubtree;
        std::vector<std::shared_ptr<QuadtreeNode>> children;
        QuadtreeNode* parent;
        std::vector<std::pair<double,double>> borders;
        double weightOfEdgeToParent;
        bool isMarked {false};

        QuadtreeNode(int dimension, double weightOfEdgeToParent)
        {
            children.resize(std::pow(2,dimension));
            for (int i = 0; i < children.size(); i++)
                children.at(i) = nullptr;
            
            borders.resize(dimension*2);
            parent = nullptr;
            this->point = point;
            this->weightOfEdgeToParent = weightOfEdgeToParent;
        }

        ~QuadtreeNode()
        {
            point.reset();
        }

        void setAllBorders(std::shared_ptr<treeNode> point, double border)
        {
            
            boost::random::mt19937 randEngine;
            boost::random::uniform_real_distribution<> r(0,1);


            double rand {0};
            while (rand == 0)
                rand = r(randEngine)*(border/2);

            if (r(randEngine) > 0.5)
                rand *= -1;

            assert(rand >= -border/2 );
            assert(rand <= border/2 );


            for(int i = 0; i < point->coordinates.size(); i++)
            {
                borders.at(i).first = point->coordinates.at(i)-border-border/2+rand;
                borders.at(i).second = point->coordinates.at(i)+border+border/2+rand;
            }
        }
        bool insertPoint(std::shared_ptr<treeNode> point)
        {
            this->pointInSubtree.push_back(point);
            std::vector<bool> childrenBits(std::pow(2,point->coordinates.size()),false);

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
            

            int childNr {0};
            for (int i = 0; i < dimension; i++)
            {
                double cutOff = (borders.at(i).first + borders.at(i).second) / 2;

                //std::cout << point->coordinates;
                if (point->coordinates.at(i) > cutOff)
                {
                    childNr += pow(2,i);
                    childrenBits.at(i) = true;
                }
            }

            if (children[childNr] == nullptr)
            {
                if (this->point != nullptr)
                {
                    /*
                    double l = this->point->coordinates.at(0);
                    double l2 = this->point->coordinates.at(1);
                    */
                    int i {0};
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
                children[childNr] = std::make_shared<QuadtreeNode>(dimension,this->weightOfEdgeToParent/2);
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
            else
            {
                children[childNr]->insertPoint(point);
            }

            assert(this->point == nullptr);

            return 0;
        }
        bool removePoint(std::shared_ptr<treeNode> point)
        {

            assert(this->point == nullptr);
            auto point_it = std::find(this->pointInSubtree.begin(), this->pointInSubtree.end(), point);
            assert(point_it != this->pointInSubtree.end());
            this->pointInSubtree.erase(point_it);

            std::vector<bool> childrenBits(std::pow(2,point->coordinates.size()),false);

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
            

            int childNr {0};
            for (int i = 0; i < dimension; i++)
            {
                double cutOff = (borders.at(i).first + borders.at(i).second) / 2;

                //std::cout << point->coordinates;
                if (point->coordinates.at(i) > cutOff)
                {
                    childNr += pow(2,i);
                    childrenBits.at(i) = true;
                }
            }
            assert(children[childNr] != nullptr);

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

            children[childNr]->removePoint(point);

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

        double calcMinDist()
        {
            double dist {0};
            QuadtreeNode* currentNode = this;

            while (currentNode->parent != nullptr)
            {
                dist += currentNode->weightOfEdgeToParent;
                currentNode = currentNode->parent;
                if (currentNode->isMarked)
                    break;
            }
            return dist;
            
        }
    };

    class quadTree
    {
    private:
        double maxDist;
        std::shared_ptr<QuadtreeNode> root;
    public:
        quadTree(std::shared_ptr<treeNode> initialPoint, int dimension, double maxDist)
        {
            root = std::make_shared<QuadtreeNode>(dimension,blaze::sqrt(dimension)*2*maxDist);
            root->setAllBorders(initialPoint, maxDist);
        }

        void insertPoint(std::shared_ptr<treeNode> point)
        {
            root->insertPoint(point);
        }

        std::shared_ptr<QuadtreeNode> findLeaf(std::shared_ptr<treeNode> point)
        {
            int dimension = point->coordinates.size();
            std::shared_ptr<QuadtreeNode> currentNode = this->root;
            int childNr {0};

            while(currentNode->point.get() != point.get())
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
                        childNr += pow(2,i);
                }
                assert(currentNode->children[childNr] != nullptr);
                currentNode = currentNode->children[childNr];
            }

            return currentNode;
        }

        QuadtreeNode* markPath(QuadtreeNode* point)
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

        double calcTreeDist(QuadtreeNode* x, treeNode* y, QuadtreeNode* ancestor)
        {
            //!!! Could also find the weight by moving up from x until y is in the subtree

            if (x->point.get() == y)
                return 0;
            int dimension = y->coordinates.size();
            int childNrX {0};
            int childNrY {0};
            double distance {0};

            while(childNrX == childNrY)
            {
                childNrX = 0;
                childNrY = 0;
                for (int i = 0; i < dimension; i++)
                {
                    double cutOff = (ancestor->borders.at(i).first + ancestor->borders.at(i).second) / 2;

                    if (x->point->coordinates.at(i) > cutOff)
                        childNrX += pow(2,i);
                    if (y->coordinates.at(i) > cutOff)
                        childNrY += pow(2,i);
                }
                assert(ancestor->children[childNrX] != nullptr);
                assert(ancestor->children[childNrY] != nullptr);

                if (childNrX == childNrY)
                    ancestor = ancestor->children[childNrX].get();
                else
                    break;
            }

            while (x->parent != ancestor)
            {
                distance += x->weightOfEdgeToParent;
                x = x->parent;
            }
            distance += x->weightOfEdgeToParent;

            return distance;
        }

    };
    
    
}