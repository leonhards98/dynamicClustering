#pragma once

#include <algorithm>
#include <vector>
#include <iostream>

#include "kmeans.hpp"
#include "coreset.hpp"
#include "random.hpp"
#include "../dynamic_coreset/coreset_point.h"

namespace coresets
{
    class SensitivitySampling
    {
    public:
        /**
         * Number of points that the algorithm should aim to include in the coreset: T
         */
        const int* TargetSamplesInCoreset;

        /**
         * Number of clusters to partition the data into: k
         */
        const size_t NumberOfClusters;

        SensitivitySampling(const size_t numberOfClusters, const int* targetSamplesInCoreset);

        std::shared_ptr<Coreset>
        run(const std::vector<dynClustering::CoresetPoint*> &data1, const std::vector<dynClustering::CoresetPoint*> &data2, const std::vector<double> &weights1, const std::vector<double> &weights2);

        utils::Random random;
    private:
        utils::IndexFinder i;

        std::shared_ptr<Coreset>
		generateCoresetPoints(const std::vector<dynClustering::CoresetPoint*> &data1, const std::vector<dynClustering::CoresetPoint*> &data2, const clustering::ClusterAssignmentList &clusterAssignments, utils::IndexFinder &i, std::vector<double> &weights);

        std::shared_ptr<blaze::DynamicVector<double>>
        calcCenterWeights(const clustering::ClusterAssignmentList &clusterAssignments, std::shared_ptr<blaze::DynamicVector<size_t>> sampledIndices, double T, std::vector<double> &weights);
    };
}
