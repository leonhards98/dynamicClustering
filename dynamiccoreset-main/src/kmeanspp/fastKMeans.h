#include "quadTree.h"

#include "../coreset_old/cluster_assignment_list.hpp"

namespace fastKMeans
{

    /**
     * @brief Finds the initial centers for a KMeans Algorithm using FastKMeans++
     * Only faster if a large number of centers need to be determined
     */
    class fastKMeansAlgo
    {
    private:
        /**
         * The sample tree. A balanced binary tree containing every datapoint in the leafs. Each edge weight corresponds to the probability
         * that this path should be chosen when determining the point that should be the next center
         */
        sampleTree *sTree;

        /**
         * The first quad tree. Distances between the points will be estimated by distances within this tree
         */
        quadTree *qTree1;
        /**
         * The second quad tree. Distances between the points will be estimated by distances within this tree
         */
        quadTree *qTree2;
        /**
         * The third quad tree. Distances between the points will be estimated by distances within this tree
         */
        quadTree *qTree3;

        /**
         * the dimenstion of the input points
         */
        int dim;

        /**
         * @brief Calculates an approximation of the maximal distance between any two points of the input data
         * @param data The input points as a matrix where each row corresponds to the coordinates of one point
         * @return
         */
        double calcMaxDist(const blaze::DynamicMatrix<double> &data);

        /**
         * @brief Selects the next center point from the input data points using the sample tree
         * @return the sample tree node that contains the chosen point.
         */
        std::shared_ptr<treeNode> multiTreeSample();

        /**
         * @brief Adds the information to a quad tree, that a point was selected as a new center. This updates the approximations
         * of the distances between the points and its closest centers
         * @param x The sample tree node corresponding to the new center
         * @param qTree the quad tree which should be updated
         * @param clusters Will contain the approximation of the assignment of the points to the closest centers by the given quad tree
         * @return (void)
         */
        void multiTreeOpen(std::shared_ptr<treeNode> x, quadTree *qTree, std::shared_ptr<clusteringStatic::ClusterAssignmentList> clusters);

    public:
        /**
         * @brief Constructor
         * @param data A matrix containing the coordinates of all points. Each row corresponds to one point
         * @param weights A vector containig the weights of all points. Each row contains the weight of the point whose coordinates
         * are stored in the same row of the matrix "data"
         */
        fastKMeansAlgo(const blaze::DynamicMatrix<double> &data, const blaze::DynamicVector<double> &weights)
        {
            srand(time(NULL));
            this->dim = data.rows();

            double maxDist = this->calcMaxDist(data);
            double M = 16 * data.columns() * maxDist;
            sTree = new sampleTree();
            auto points = sTree->initTree(data, weights, M);
            qTree1 = new quadTree(points.at(rand() % (points.size() - 1)), data.columns(), maxDist);
            qTree2 = new quadTree(points.at(rand() % (points.size() - 1)), data.columns(), maxDist);
            qTree3 = new quadTree(points.at(rand() % (points.size() - 1)), data.columns(), maxDist);
            for (int i = 0; i < points.size(); i++)
            {
                qTree1->insertPoint(points.at(i));
                qTree2->insertPoint(points.at(i));
                qTree3->insertPoint(points.at(i));
            }
        }

        /**
         * @brief Destructor
         */
        ~fastKMeansAlgo()
        {
            delete qTree1;
            delete qTree2;
            delete qTree3;
            delete sTree;
        }

        /**
         * @brief Runs the FastKMeans++ algorithm and returs a cluster assignment that contains both the points that were
         * selected as new centers and the assignment of all input points to the clusters
         * @param k Number of clusters to find
         * @return culster assignment containing the centers and assingments
         */
        std::shared_ptr<clusteringStatic::ClusterAssignmentList> fastKMeansPP(int k);

        /**
         * @brief dynamically inserts a new point into the system (sample tree and all quad trees)
         * @param coordinates vector containing the coordinates of the new point
         * @param pointWeight weight of the new point
         * @param index unique index of the new point
         * @return (void)
         */
        void insertPoint(const blaze::DynamicVector<double> &coordinates, double pointWeight, double index);

        /**
         * @brief removes a point from the system, given an index
         * @param index index of the point to remove
         * @return (void)
         */
        void deletePoint(int index);
    };

}
