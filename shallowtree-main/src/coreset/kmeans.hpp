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
#include "../dynamic_coreset/coreset_point.h"

namespace clustering
{
    /**
     * @brief Implementation of the k-Means clustering algorithm.
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
        KMeans(size_t numOfClusters, utils::IndexFinder i);

        /**
         * @brief Runs the algorithm.
         * @param data A NxD data matrix containing N data points where each point has D dimensions.
         */
        std::shared_ptr<ClusteringResult>
        run(const std::vector<dynClustering::CoresetPoint*> &data1, const std::vector<dynClustering::CoresetPoint*> &data2, const std::vector<double> &weights);

        /**
         * @brief Picks `k` points as the initial centers using the k-Means++ initialisation procedure.
         * @param dataMatrix A NxD data matrix containing N data points where each point has D dimensions.
         */
        std::shared_ptr<clustering::ClusterAssignmentList>
        pickInitialCentersViaKMeansPlusPlus(const std::vector<dynClustering::CoresetPoint*> &data1, const std::vector<dynClustering::CoresetPoint*> &data2, const std::vector<double> &weights);

        blaze::DynamicMatrix<double>
		copyRows(const std::vector<dynClustering::CoresetPoint*> &data1, const std::vector<dynClustering::CoresetPoint*> &data2, const std::vector<size_t> &indicesToCopy);

    private:
        const size_t NumOfClusters;
        utils::IndexFinder i;

        /**
         * @brief Run Lloyd's algorithm to perform the clustering of data points.
         * @param dataMatrix A NxD data matrix containing N data points where each point has D dimensions.
         * @param dataMatrix Initial k centroids where k is the number of required clusters.
         */
        std::shared_ptr<ClusteringResult>
		runLloydsAlgorithm(const std::vector<dynClustering::CoresetPoint*> &data1, const std::vector<dynClustering::CoresetPoint*> &data2, blaze::DynamicMatrix<double> centroids, const std::vector<double> &weights);
    };

}
