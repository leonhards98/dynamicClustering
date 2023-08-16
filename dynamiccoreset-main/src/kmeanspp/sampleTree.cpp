#include "sampleTree.h"
#include <queue>

namespace fastKMeans
{
    std::shared_ptr<treeNode> sampleTree::insertPoint(blaze::DynamicVector<double> cooridinates, double pointWeight, double index)
    {
        // Pop the next leaf from "insertionLeafs".
        // Two leafs will be added as children to the leaf recieved that way. One will be populated with the data from the
        // original leaf and  the other one will be be used to store the new point. The new leafs will be pushed
        // to the "insertionLeafs" vector
        treeNode *parent = insertionLeafs.front().get();
        insertionLeafs.pop();

        parent->lChild = std::make_shared<treeNode>();
        parent->rChild = std::make_shared<treeNode>();

        // Data of the original leaf (now parent) will be stored in the left child
        parent->lChild->index = parent->index;
        parent->lChild->KMeansWeight = parent->KMeansWeight;
        parent->lChild->parent = parent;
        parent->lChild->coordinates = parent->coordinates;
        parent->lChild->pointWeight = parent->pointWeight;

        // Data of the new point will be stored in the right child
        parent->rChild->index = index;
        parent->rChild->coordinates.resize(cooridinates.size());
        for (int i = 0; i < cooridinates.size(); i++)
            parent->rChild->coordinates.at(i) = cooridinates.at(i);

        parent->rChild->pointWeight = pointWeight;
        parent->rChild->KMeansWeight = 0;
        parent->rChild->parent = parent;

        // Reset the parent
        parent->coordinates.resize(0);
        parent->index = -1;
        parent->KMeansWeight = 0;

        insertionLeafs.push(parent->lChild);
        insertionLeafs.push(parent->rChild);

        return parent->rChild;
    }

    void sampleTree::removePoint(std::shared_ptr<treeNode> point)
    {
        std::shared_ptr<treeNode> rightmostLeaf = this->searchRightmostLeaf();

        treeNode *parent = rightmostLeaf->parent;
        treeNode *grandParent = parent->parent;
        assert(grandParent != nullptr);
        std::shared_ptr<treeNode> sibling = parent->lChild;

        // If the leaf to remove is the sibling of the rightmost leaf, just remove
        // the leaf and swap the rightmost leaf with its parent. Then remove the parent
        // Also need to update the "insertionLeafs" queue accordingly
        if (sibling.get() == point.get())
        {
            double oldWeight = sibling->KMeansWeight + point->KMeansWeight;
            parent->lChild.reset();
            parent->rChild.reset();

            if (parent == grandParent->lChild.get())
            {
                grandParent->lChild.reset();
                grandParent->lChild = rightmostLeaf;
            }
            else
            {
                grandParent->rChild.reset();
                grandParent->rChild = rightmostLeaf;
            }
            rightmostLeaf->parent = grandParent;

            std::queue<std::shared_ptr<treeNode>> insertionLeafsSave;

            insertionLeafs.swap(insertionLeafsSave);

            insertionLeafs.push(insertionLeafsSave.back());
            int size = insertionLeafsSave.size();
            for (int i = 0; i < size - 2; i++)
            {
                insertionLeafs.push(insertionLeafsSave.front());
                insertionLeafsSave.pop();
            }
            this->updateInternalWeights(sibling.get(), sibling->KMeansWeight, oldWeight);

            return;
        }
        // Swap the leaf to be removed with the rightmost leaf. Also need to swap
        // the according pointers in the "insertionLeaf" vector
        if (rightmostLeaf.get() != point.get())
        {
            this->swapLeafs(point, rightmostLeaf);
            int size = insertionLeafs.size();

            this->updateInternalWeights(rightmostLeaf.get(), point->KMeansWeight, rightmostLeaf->KMeansWeight);
            this->updateInternalWeights(point.get(), rightmostLeaf->KMeansWeight, point->KMeansWeight);
            for (int i = 0; i < size; i++)
            {
                auto currentLeaf = insertionLeafs.front();
                insertionLeafs.pop();
                if (currentLeaf == point)
                {
                    insertionLeafs.push(rightmostLeaf);
                    continue;
                }
                if (currentLeaf == rightmostLeaf)
                {
                    insertionLeafs.push(point);
                    continue;
                }
                insertionLeafs.push(currentLeaf);
            }

            assert(insertionLeafs.size() == size);
        }
        double oldWeight = sibling->KMeansWeight + point->KMeansWeight;

        // Remove the rightmost leaf (the leaf to remove. swapped if neccersary)
        // and swap its sibling with its parent. Then remove the parent.
        // Also need to correct the pointers in "insertionLeaf"
        if (parent == grandParent->lChild.get())
        {
            grandParent->lChild.reset();
            grandParent->lChild = sibling;
        }
        else
        {
            grandParent->rChild.reset();
            grandParent->rChild = sibling;
        }
        sibling->parent = grandParent;

        this->updateInternalWeights(sibling.get(), sibling->KMeansWeight, oldWeight);

        std::queue<std::shared_ptr<treeNode>> insertionLeafsSave;

        int size = insertionLeafs.size();
        for (int i = 0; i < size - 1; i++)
        {
            insertionLeafsSave.push(insertionLeafs.front());
            insertionLeafs.pop();
        }

        insertionLeafs.push(insertionLeafsSave.back());
        size = insertionLeafsSave.size();
        for (int i = 0; i < size - 1; i++)
        {
            insertionLeafs.push(insertionLeafsSave.front());
            insertionLeafsSave.pop();
        }

        return;
    }

    std::vector<std::shared_ptr<treeNode>> sampleTree::initTree(const blaze::DynamicMatrix<double> &data, const blaze::DynamicVector<double> &weights, double M)
    {
        // Insert the first point into the root node
        double weight = weights.at(0);
        double index = 0;
        blaze::DynamicVector<double> coordinates;
        coordinates.resize(data.columns());
        for (int j = 0; j < data.columns(); j++)
        {
            coordinates.at(j) = data.at(0, j);
        }
        this->root->coordinates = coordinates;
        this->root->index = index;
        this->root->pointWeight = weight;

        // Insert all remaining points by calling the "insertPoint" function
        for (int i = 1; i < data.rows(); i++)
        {
            double weight = weights.at(i);
            double index = i;
            blaze::DynamicVector<double> coordinates;
            coordinates.resize(data.columns());
            for (int j = 0; j < data.columns(); j++)
            {
                coordinates.at(j) = data.at(i, j);
            }
            insertPoint(coordinates, weight, index);
        }

        this->resetWeightsWithPointWeights(M);

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

    std::shared_ptr<treeNode> sampleTree::searchRightmostLeaf()
    {
        std::queue<std::shared_ptr<treeNode>> queue;

        queue.push(this->root);

        std::shared_ptr<treeNode> currentNode = nullptr;
        while (!queue.empty())
        {
            currentNode = queue.front();
            queue.pop();
            if (currentNode->lChild == nullptr && currentNode->rChild == nullptr)
                continue;
            queue.push(currentNode->lChild);
            queue.push(currentNode->rChild);
        }

        return currentNode;
    }

    bool sampleTree::swapLeafs(std::shared_ptr<treeNode> leaf1, std::shared_ptr<treeNode> leaf2)
    {

        treeNode *parent1 = leaf1->parent;
        treeNode *parent2 = leaf2->parent;

        if (parent1->lChild.get() == leaf1.get())
        {
            parent1->lChild.reset();
            parent1->lChild = leaf2;
        }
        else
        {
            parent1->rChild.reset();
            parent1->rChild = leaf2;
        }

        if (parent2->lChild.get() == leaf2.get())
        {
            parent2->lChild.reset();
            parent2->lChild = leaf1;
        }
        else
        {
            parent2->rChild.reset();
            parent2->rChild = leaf1;
        }

        leaf1->parent = parent2;

        leaf2->parent = parent1;

        return true;
    }

    int sampleTree::getTreeDepth()
    {
        int ret{0};
        std::shared_ptr<treeNode> currentNode = this->root;
        while (currentNode->lChild != nullptr)
        {
            ret++;
            currentNode = currentNode->lChild;
        }

        return ret;
    }

    const std::shared_ptr<treeNode> &sampleTree::getRoot()
    {
        return this->root;
    }

    void sampleTree::updateWeights(treeNode *leaf, double newDist, double oldDist)
    {
        assert(leaf->lChild == nullptr);
        assert(leaf->rChild == nullptr);

        newDist = newDist * leaf->pointWeight;
        treeNode *currentNode = leaf;

        while (currentNode->parent != nullptr)
        {
            currentNode->KMeansWeight = currentNode->KMeansWeight - oldDist + newDist;
            currentNode = currentNode->parent;
        }
    }

    void sampleTree::updateInternalWeights(treeNode *node, double newDist, double oldDist)
    {

        if (node->parent == nullptr)
            return;
        treeNode *currentNode = node->parent;

        while (currentNode->parent != nullptr)
        {
            currentNode->KMeansWeight = currentNode->KMeansWeight - oldDist + newDist;
            currentNode = currentNode->parent;
        }
    }

    void sampleTree::resetWeightsWithPointWeights(double M)
    {
        std::vector<std::shared_ptr<treeNode>> leafs = this->getAllLeafs();

        // Loop over all leafs and then move up the path over all its anchestors and add the edge weight to all
        // other edges on the path.
        for (int i = 0; i < leafs.size(); i++)
        {
            treeNode *leaf = leafs.at(i).get();
            leaf->KMeansWeight = leaf->pointWeight * M;
            double w{0};

            while (leaf->parent != nullptr)
            {
                w = leaf->KMeansWeight;
                leaf = leaf->parent;
                leaf->KMeansWeight += w;
            }
        }
    }
}
