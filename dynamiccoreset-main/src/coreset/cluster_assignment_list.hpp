#pragma once

#include <string>
#include <iostream>
#include <set>

#include <blaze/Math.h>
#include <boost/array.hpp>
#include <boost/range/algorithm_ext/erase.hpp>
#include "../dynamic_coreset/coreset_point.h"
#include "distances.hpp"

namespace clustering
{
    /**
     * @brief Represents a collection of points-to-cluster assignments.
     */
    class ClusterAssignmentList
    {
    public:
        /**
         * @brief Creates a new instance of ClusterAssignmentList.
         * @param numOfPoints The total number of points in the dataset.
         * @param numOfClusters The number of clusters that are generated.
         */
        ClusterAssignmentList(size_t numOfPoints, size_t numOfClusters);

        /**
         * @brief Copy constructor.
         * @param other The other cluster assignments to copy from.
         */
        ClusterAssignmentList(const ClusterAssignmentList &other);

        /**
         * @brief Copies cluster assignments from another object.
         * @param other The other cluster assignments to copy from.
         */
        ClusterAssignmentList &
        operator=(const ClusterAssignmentList &other);

        /**
         * @brief Assign all data points to their closest centers.
         */
        void
        assignAll(const std::vector<dynClustering::CoresetPoint *> &data1, const std::vector<dynClustering::CoresetPoint *> &data2, const blaze::DynamicMatrix<double> &centers, utils::IndexFinder &i, const std::vector<double> &weights);

        /**
         * @brief Assign a point to a cluster.
         * @param pointIndex The index of the point to assign the cluster to.
         * @param clusterIndex The index of the cluster to assign the point to.
         * @param distance The distance between the point and the cluster.
         */
        void
        assign(size_t pointIndex, size_t clusterIndex, double distance, double weightedDistance);

        /**
         * @brief Gets the assigned cluster index of a point.
         * @param pointIndex The index of the point for which to return the cluster index.
         */
        size_t
        getCluster(size_t pointIndex) const;

        /**
         * @brief Returns the total number of points in the dataset.
         */
        size_t
        getNumberOfPoints() const;

        /**
         * @brief Returns the number of clusters that are generated.
         */
        size_t
        getNumberOfClusters() const;

        /**
         * @brief Returns the distance of each point to its assigned cluster's centroid.
         */
        blaze::DynamicVector<double> &
        getCentroidDistances();

        /**
         * @brief Returns the number of points in a cluster.
         */
        size_t
        countPointsInCluster(size_t clusterIndex) const;

        size_t
        countWeightedPointsInCluster(size_t clusterIndex, std::vector<double> &weights) const;

        /**
         * @brief Returns the total cost of the cluster assignments.
         *
         * This cost is the sum of pairwise distances between each point and its closest center.
         */
        double
        getTotalCost() const;

        /**
         * @brief Returns the total weighted cost of the cluster assignments.
         *
         * This cost is the sum of pairwise distances between each point and
         * its closest center times the weight of the point.
         */
        double
        getTotalWeightedCost() const;

        /**
         * @brief Returns the cost of a point i.e.: cost(p, A).
         *
         * The cost is the distance between a point and its clostest center.
         */
        double
        getPointCost(size_t pointIndex) const;

        /**
         * @brief Returns the weighted cost of a point i.e.: cost(p, A).
         *
         * The cost is the distance between a point and its clostest center
         * times the weight of the point.
         */
        double
        getWeightedPointCost(size_t pointIndex) const;

        /**
         * @brief Returns the average costs of each cluster.
         */
        std::shared_ptr<blaze::DynamicVector<double>>
        calcAverageClusterCosts() const;

        /**
         * @brief Returns the sum of costs of points in each cluster i.e. cost(C_i) = sum_{p in C_i} cost(p, A).
         */
        std::shared_ptr<blaze::DynamicVector<double>>
        calcClusterCosts() const;

        /**
         * @brief Returns the normalized cost of each point i.e.: cost(p, A) / sum_p { cost(p, A) }.
         */
        blaze::DynamicVector<double>
        getNormalizedCosts() const;

        /**
         * @brief Returns the normalized weighted cost of each point i.e.: weightedCost(p, A) / sum_p { weightedcost(p, A) }.
         */
        blaze::DynamicVector<double>
        getWeightedNormalizedCosts() const;

        /**
         * @brief Compute cluster centers using current assignments.
         */
        void
        calcCenters(const std::vector<dynClustering::CoresetPoint *> &data1, const std::vector<dynClustering::CoresetPoint *> &data2, blaze::DynamicMatrix<double> &newCenters, utils::IndexFinder &i);

        std::shared_ptr<blaze::DynamicVector<double>>
        calcCenter(const std::vector<dynClustering::CoresetPoint *> &data1, const std::vector<dynClustering::CoresetPoint *> &data2, const size_t clusterIndex, utils::IndexFinder &i, std::vector<double> &weights) const;

        std::shared_ptr<std::vector<size_t>>
        getPointsByCluster(size_t clusterIndex) const;

        std::shared_ptr<std::vector<size_t>>
        getClusterIndices() const;

    private:
        /**
         * The total number of points in the dataset.
         */
        size_t numOfPoints;

        /**
         * The number of clusters that are generated.
         */
        size_t numOfClusters;

        /**
         * A vector of size N contain the cluster index for each point.
         */
        blaze::DynamicVector<size_t> clusters;

        /**
         * A vector of size N containing the distance between the
         * assigned cluster of each point in the dataset.
         */
        blaze::DynamicVector<double> distances;

        /**
         * A vector of size N containing the distance between the
         * assigned cluster of each point in the dataset times the weight of the point.
         */
        blaze::DynamicVector<double> weightedDistances;
    };

}
