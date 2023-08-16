#include "../dynamic_coreset/coreset_point.h"
#include <memory>
#include <queue>
namespace fastKMeans
{
    struct treeNode
    {
        blaze::DynamicVector<double> coordinates;
        double pointWeight;
        double KMeansWeight {0};
        int index;
        std::shared_ptr<treeNode> lChild;
        std::shared_ptr<treeNode> rChild;
        treeNode* parent;

        void resetWeightWithoutPointWeight(int nrPoints, double M)
        {
            assert(nrPoints >= 1);
            this->KMeansWeight = nrPoints*M;

            if (nrPoints == 1)
                return;

            double rest = nrPoints - pow(2,blaze::floor(blaze::log2(nrPoints)));
            if (rest == 0)
            {
                lChild->resetWeightWithoutPointWeight(nrPoints/2,M);
                rChild->resetWeightWithoutPointWeight(nrPoints/2,M);
                return;
            }

            double halve = pow(2,blaze::floor(blaze::log2(nrPoints))) / 2;

            if (rest <= halve)
            {
                lChild->resetWeightWithoutPointWeight(halve+rest,M);
                rChild->resetWeightWithoutPointWeight(nrPoints-(halve+rest),M);
            }
            else
            {
                lChild->resetWeightWithoutPointWeight(halve*2,M);
                rChild->resetWeightWithoutPointWeight(nrPoints-(halve*2),M);
            }
            return;
        }

    };

    class sampleTree
    {
    private:
        std::shared_ptr<treeNode> root;
        std::queue<std::shared_ptr<treeNode>> insertionLeafs ;


        
    public:
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

        std::vector<std::shared_ptr<treeNode>> initTree(const blaze::DynamicMatrix<double> &data, const blaze::DynamicVector<double> &weights, double M);
        std::shared_ptr<treeNode> searchPoint(int index);
        std::vector<std::shared_ptr<treeNode>> getAllLeafs();
        std::shared_ptr<treeNode> searchLeftmostLeaf();
        int getTreeDepth();
        const std::shared_ptr<treeNode>& getRoot();
        void updateWeights(treeNode* leaf, double newDist, double oldDist);
        void updateInternalWeights(treeNode* node, double newDist, double oldDist);
        void resetWeightsWithPointWeights(double M);
        std::shared_ptr<treeNode> insertPoint(blaze::DynamicVector<double> cooridinates, double pointWeight, double index);
        void removePoint(std::shared_ptr<treeNode> point);
    };
    
    

}