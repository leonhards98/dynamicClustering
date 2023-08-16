#include "quadTree.h"

#include "../coreset_old/cluster_assignment_list.hpp"

namespace fastKMeans{


    
    class fastKMeansAlgo
    {
    private:
        sampleTree* sTree;
        quadTree* qTree1;
        quadTree* qTree2;
        quadTree* qTree3;
        int dim;

    public:
        fastKMeansAlgo(const blaze::DynamicMatrix<double> &data, const blaze::DynamicVector<double> &weights)
        {
            srand (time(NULL));
            this->dim = data.rows();


            double maxDist = this->calcMaxDist(data);
            double M = 16*data.columns()*maxDist;
            sTree = new sampleTree();
            auto points = sTree->initTree(data, weights, M);
            qTree1 = new quadTree(points.at(rand()%(points.size()-1)),data.columns(),maxDist);
            qTree2 = new quadTree(points.at(rand()%(points.size()-1)),data.columns(),maxDist);
            qTree3 = new quadTree(points.at(rand()%(points.size()-1)),data.columns(),maxDist);
            for (int i = 0; i < points.size(); i++)
            {
                qTree1->insertPoint(points.at(i));
                qTree2->insertPoint(points.at(i));
                qTree3->insertPoint(points.at(i));
            }
        }
        ~fastKMeansAlgo()
        {
            delete qTree1;
            delete qTree2;
            delete qTree3;
            delete sTree;
        }

        double calcMaxDist(const blaze::DynamicMatrix<double> &data);
        std::shared_ptr<treeNode> multiTreeSample();
        void multiTreeOpen(std::shared_ptr<treeNode> x, quadTree* qTree, std::shared_ptr<clusteringStatic::ClusterAssignmentList> clusters);
        void insertPoint(const blaze::DynamicVector<double> &coordinates, double pointWeight, double index);

        std::shared_ptr<clusteringStatic::ClusterAssignmentList> fastKMeansPP(int k);
    };

}