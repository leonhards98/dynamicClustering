#include "sampleTree.h"


namespace fastKMeans
{
    std::shared_ptr<treeNode> sampleTree::insertPoint(blaze::DynamicVector<double> cooridinates, double pointWeight, double index)
    {
        treeNode* parent = insertionLeafs.front().get();
        insertionLeafs.pop();

        parent->lChild = std::make_shared<treeNode>();
        parent->rChild = std::make_shared<treeNode>();

        parent->lChild->index = parent->index;
        parent->lChild->KMeansWeight = parent->KMeansWeight;
        parent->lChild->parent = parent;
        parent->lChild->coordinates = parent->coordinates;
        parent->lChild->pointWeight = parent->pointWeight;

        parent->rChild->index = index;

        parent->rChild->coordinates.resize(cooridinates.size());
        for (int i = 0; i < cooridinates.size(); i++)
            parent->rChild->coordinates.at(i) = cooridinates.at(i);

        parent->rChild->pointWeight = pointWeight;
        parent->rChild->KMeansWeight = 0;
        parent->rChild->parent = parent;

        parent->coordinates.resize(0);
        parent->index = -1;
        parent->KMeansWeight = 0;

        insertionLeafs.push(parent->lChild);
        insertionLeafs.push(parent->rChild);

        return parent->rChild;
    }


    void sampleTree::removePoint(std::shared_ptr<treeNode> point)
    {
        

    }



    std::vector<std::shared_ptr<treeNode>> sampleTree::initTree(const blaze::DynamicMatrix<double> &data, const blaze::DynamicVector<double> &weights, double M)
    {
        double weight = weights.at(0);
        double index = 0;
        blaze::DynamicVector<double> coordinates;
        coordinates.resize(data.columns());
        for (int j = 0; j < data.columns(); j++)
        {
            coordinates.at(j) = data.at(0,j);
        }
        this->root->coordinates = coordinates;
        this->root->index = index;
        this->root->pointWeight = weight;

        for (int i = 1; i < data.rows(); i++)
        {
            double weight = weights.at(i);
            double index = i;
            blaze::DynamicVector<double> coordinates;
            coordinates.resize(data.columns());
            for (int j = 0; j < data.columns(); j++)
            {
                coordinates.at(j) = data.at(i,j);
            }
            insertPoint(coordinates,weight,index);
        }


        this->resetWeightsWithPointWeights(M);
        //this->root->resetWeightWithoutPointWeight(data.rows(),M);

        return this->getAllLeafs();


    }
    std::shared_ptr<treeNode> sampleTree::searchPoint(int index)
    {
        std::queue<std::shared_ptr<treeNode>> queue;

        queue.push(this->root);

        std::shared_ptr<treeNode> currentNode = nullptr;
        while (!queue.empty())
        {
            currentNode = queue.front();
            queue.pop();
            if (currentNode->index == index)
                return currentNode;
            queue.push(currentNode->lChild);
            queue.push(currentNode->rChild);
        }

        return nullptr;

    }
    std::vector<std::shared_ptr<treeNode>> sampleTree::getAllLeafs()
    {
        std::queue<std::shared_ptr<treeNode>> queue;

        queue.push(this->root);

        std::shared_ptr<treeNode> currentNode = nullptr;
        std::vector<std::shared_ptr<treeNode>> ret;
        int i = 2;
        while (!queue.empty())
        {

            currentNode = queue.front();
            queue.pop();

            //std::cout << currentNode << ": " << currentNode->index << "|";
            //if (blaze::floor(blaze::log2(i)) == blaze::log2(i))
                //std::cout << "\n";
            //i++;

            if (currentNode->lChild == nullptr && currentNode->rChild == nullptr)
            {
                ret.push_back(currentNode);
                continue;
            } 
            queue.push(currentNode->lChild);
            queue.push(currentNode->rChild);
        }

        return ret;

    }


    std::shared_ptr<treeNode> sampleTree::searchLeftmostLeaf()
    {
        std::queue<std::shared_ptr<treeNode>> queue;

        queue.push(this->root);

        std::shared_ptr<treeNode> currentNode = nullptr;
        while (!queue.empty())
        {
            currentNode = queue.front();
            queue.pop();
            if (currentNode->lChild == nullptr && currentNode->rChild == nullptr)
                return currentNode;
            queue.push(currentNode->lChild);
            queue.push(currentNode->rChild);
        }

        return currentNode;

    }

    int sampleTree::getTreeDepth()
    {
        int ret {0};
        std::shared_ptr<treeNode> currentNode = this->root;
        while (currentNode->lChild != nullptr)
        {
            ret++;
            currentNode = currentNode->lChild;
        }

        return ret;
    }

    const std::shared_ptr<treeNode>& sampleTree::getRoot()
    {
        return this->root;
    }

    void sampleTree::updateWeights(treeNode* leaf, double newDist, double oldDist)
    {
        assert(leaf->lChild == nullptr);
        assert(leaf->rChild == nullptr);

        newDist = newDist*leaf->pointWeight;
        treeNode* currentNode = leaf;

        while (currentNode->parent != nullptr)
        {
            currentNode->KMeansWeight = currentNode->KMeansWeight-oldDist+newDist;
            currentNode = currentNode->parent;
        }
    }

    void sampleTree::updateInternalWeights(treeNode* node, double newDist, double oldDist)
    {
        assert(node->lChild != nullptr);
        assert(node->rChild != nullptr);

        if (node->parent == nullptr)
            return;
        treeNode* currentNode = node->parent;

        while (currentNode->parent != nullptr)
        {
            currentNode->KMeansWeight = currentNode->KMeansWeight-oldDist+newDist;
            currentNode = currentNode->parent;
        }
    }
    void sampleTree::resetWeightsWithPointWeights(double M)
    {
        std::vector<std::shared_ptr<treeNode>> leafs = this->getAllLeafs();

        for (int i = 0; i < leafs.size(); i++)
        {
            treeNode* leaf = leafs.at(i).get();
            leaf->KMeansWeight = leaf->pointWeight*M;
            double w {0};

            while (leaf->parent != nullptr)
            {
                w = leaf->KMeansWeight;
                leaf = leaf->parent;
                leaf->KMeansWeight += w;
            }
        }

    }
}