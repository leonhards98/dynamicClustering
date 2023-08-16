#include "fastKMeans.h"
#include <memory>
#include <shared_mutex>

using namespace fastKMeans;

double fastKMeansAlgo::calcMaxDist(const blaze::DynamicMatrix<double> &data)
{
    int dim = data.columns();

    blaze::DynamicVector<double> firstPoint;
    firstPoint.resize(dim);

    for (int i = 0; i < dim; i++)
        firstPoint.at(i) = data.at(0, i);

    double maxDistSquared{0};
    for (int i = 1; i < data.rows(); i++)
    {
        double distSquared{0};
        for (int d = 0; d < dim; d++)
        {
            distSquared += (data.at(i, d) - firstPoint.at(d)) * (data.at(i, d) - firstPoint.at(d));
        }
        if (distSquared > maxDistSquared)
            maxDistSquared = distSquared;
    }

    return sqrt(maxDistSquared) * 2;
}

std::shared_ptr<treeNode> fastKMeansAlgo::multiTreeSample()
{
    std::shared_ptr<treeNode> currentNode = this->sTree->getRoot();

    while (currentNode->lChild != nullptr && currentNode->rChild != nullptr)
    {
        if (static_cast<double>(rand()) / RAND_MAX <= (currentNode->lChild->KMeansWeight / (currentNode->lChild->KMeansWeight + currentNode->rChild->KMeansWeight)))
            currentNode = currentNode->lChild;
        else
            currentNode = currentNode->rChild;
    }

    return currentNode;
}

void fastKMeansAlgo::multiTreeOpen(std::shared_ptr<treeNode> x, quadTree *qTree, std::shared_ptr<clusteringStatic::ClusterAssignmentList> clusters)
{
    std::shared_ptr<QuadtreeNode> leaf = qTree->findLeaf(x);
    QuadtreeNode *topmostNode = qTree->markPath(leaf.get());

    for (const auto &y : topmostNode->pointInSubtree)
    {
        double newDist = qTree->calcTreeDist(leaf.get(), y.get(), topmostNode);
        if (y->KMeansWeight / y->pointWeight > newDist)
        {
            clusters->assign(y->index, leaf->point->index, newDist);
            this->sTree->updateWeights(y.get(), newDist, y->KMeansWeight);
        }
    }
}

std::shared_ptr<clusteringStatic::ClusterAssignmentList> fastKMeansAlgo::fastKMeansPP(int k)
{
    std::shared_ptr<clusteringStatic::ClusterAssignmentList> clusters = std::make_shared<clusteringStatic::ClusterAssignmentList>(this->dim, k);

    for (int i = 0; i < k; i++)
    {
        std::shared_ptr<treeNode> x = this->multiTreeSample();
        this->multiTreeOpen(x, this->qTree1, clusters);
        this->multiTreeOpen(x, this->qTree2, clusters);
        this->multiTreeOpen(x, this->qTree3, clusters);
    }
    return clusters;
}

void fastKMeansAlgo::insertPoint(const blaze::DynamicVector<double> &coordinates, double pointWeight, double index)
{
    std::shared_ptr<treeNode> newLeaf = this->sTree->insertPoint(coordinates, pointWeight, index);
    std::shared_ptr<QuadtreeNode> newQuadLeaf;

    this->qTree1->insertPoint(newLeaf);
    newQuadLeaf = qTree1->findLeaf(newLeaf);
    newLeaf->KMeansWeight = newQuadLeaf->calcMinDist() * newLeaf->pointWeight;
    newLeaf->parent->KMeansWeight = newLeaf->parent->lChild->KMeansWeight + newLeaf->parent->rChild->KMeansWeight;
    this->sTree->updateInternalWeights(newLeaf->parent->parent, newLeaf->parent->KMeansWeight, newLeaf->parent->lChild->KMeansWeight);

    this->qTree2->insertPoint(newLeaf);
    newQuadLeaf = qTree2->findLeaf(newLeaf);
    newLeaf->KMeansWeight = newQuadLeaf->calcMinDist() * newLeaf->pointWeight;
    newLeaf->parent->KMeansWeight = newLeaf->parent->lChild->KMeansWeight + newLeaf->parent->rChild->KMeansWeight;
    this->sTree->updateInternalWeights(newLeaf->parent->parent, newLeaf->parent->KMeansWeight, newLeaf->parent->lChild->KMeansWeight);

    this->qTree3->insertPoint(newLeaf);
    newQuadLeaf = qTree3->findLeaf(newLeaf);
    newLeaf->KMeansWeight = newQuadLeaf->calcMinDist() * newLeaf->pointWeight;
    newLeaf->parent->KMeansWeight = newLeaf->parent->lChild->KMeansWeight + newLeaf->parent->rChild->KMeansWeight;
    this->sTree->updateInternalWeights(newLeaf->parent->parent, newLeaf->parent->KMeansWeight, newLeaf->parent->lChild->KMeansWeight);
}

void fastKMeansAlgo::deletePoint(int index)
{
    std::shared_ptr<treeNode> pointToRemove = sTree->searchPoint(index);
    std::shared_ptr<QuadtreeNode> leafToRemove;

    qTree1->removePoint(pointToRemove);
    qTree2->removePoint(pointToRemove);
    qTree3->removePoint(pointToRemove);

    sTree->removePoint(pointToRemove);

    return;
}
