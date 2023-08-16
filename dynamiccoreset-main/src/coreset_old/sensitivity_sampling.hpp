#pragma once

#include <algorithm>
#include <vector>
#include <iostream>

#include "kmeans.hpp"
#include "coreset.hpp"
#include "random.hpp"

namespace coresetsStatic
{
    class SensitivitySampling
    {
    public:
        /**
         * Number of points that the algorithm should aim to include in the coreset: T
         */
        const size_t TargetSamplesInCoreset;

        /**
         * Number of clusters to partition the data into: k
         */
        const size_t NumberOfClusters;

        SensitivitySampling(size_t numberOfClusters, size_t targetSamplesInCoreset);

        std::shared_ptr<Coreset>
        run(const blaze::DynamicMatrix<double> &data, const blaze::DynamicVector<double> &weights);

        utilsStatic::Random random;
    private:

        std::shared_ptr<Coreset>
        generateCoresetPoints(const blaze::DynamicMatrix<double> &data, const clusteringStatic::ClusterAssignmentList &clusterAssignments);

        std::shared_ptr<blaze::DynamicVector<double>>
        calcCenterWeights(const clusteringStatic::ClusterAssignmentList &clusterAssignments, std::shared_ptr<blaze::DynamicVector<size_t>> sampledIndices);
    };
}
