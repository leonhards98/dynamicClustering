#pragma once

#include <memory>
#include <iostream>
#include <random>
#include <string>

#include <blaze/Math.h>
#include <boost/array.hpp>
#include <boost/range/algorithm_ext/erase.hpp>

#include "cluster_assignment_list.hpp"
#include "clustering_result.hpp"
#include "random.hpp"
#include "stop_watch.hpp"
#include "distances.hpp"
#include "../kmeanspp/fastKMeans.h"

namespace clusteringStatic
{
    /**
     * @brief Implementation of the k-Means clusteringStatic algorithm.
     */
    class KMeans
    {
    public:
        /**
         * @brief Creates a new instance of KMeans.
         * @param numOfClusters The number of clusters to generate.
         * @param initKMeansPlusPlus Initialise centroids using k-Means++.
         * @param maxIterations Maximum number of iterations.
         * @param convergenceDiff The difference in the norms of the centroids when to stop k-Means iteration.
         */
        KMeans(size_t numOfClusters, bool initKMeansPlusPlus = true, size_t maxIterations = 300, double convergenceDiff = 0.0001);

        /**
         * @brief Runs the algorithm.
         * @param data A NxD data matrix containing N data points where each point has D dimensions.
         * @param weights A N data vector containting the weights of the points. Assignment defined by the order
         * @param newKMeansPP If the fast KMeansPP Algorithm sould be used. In only faster for a large number of centers
         */
        std::shared_ptr<ClusteringResult>
        run(const blaze::DynamicMatrix<double> &data, const blaze::DynamicVector<double> &weights, bool newKMeansPP);

        /**
         * @brief Picks `k` points as the initial centers using the k-Means++ initialisation procedure.
         * @param dataMatrix A NxD data matrix containing N data points where each point has D dimensions.
         */
        std::shared_ptr<clusteringStatic::ClusterAssignmentList>
        pickInitialCentersViaKMeansPlusPlus(const blaze::DynamicMatrix<double> &dataMatrix, const blaze::DynamicVector<double> &weights);

        blaze::DynamicMatrix<double>
        copyRows(const blaze::DynamicMatrix<double> &data, const std::vector<size_t> &indicesToCopy);

        /**
         * @brief Run Lloyd's algorithm to perform the clustering of data points.
         * @param dataMatrix A NxD data matrix containing N data points where each point has D dimensions.
         * @param dataMatrix Initial k centroids where k is the number of required clusters.
         */
        std::shared_ptr<ClusteringResult>
        runLloydsAlgorithm(const blaze::DynamicMatrix<double> &dataMatrix, blaze::DynamicMatrix<double> initialCentroids, const blaze::DynamicVector<double> &weights);

    private:
        const size_t NumOfClusters;
        const bool InitKMeansPlusPlus;
        const size_t MaxIterations;
        const double ConvergenceDiff;

    };

}
