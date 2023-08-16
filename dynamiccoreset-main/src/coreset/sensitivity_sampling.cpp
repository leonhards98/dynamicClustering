#include "sensitivity_sampling.hpp"
#include "distances.hpp"
#include "kmeans.hpp"
#include <vector>
#include <cassert>
#include <numeric>

using namespace coresets;

SensitivitySampling::SensitivitySampling(const size_t numberOfClusters, const int* targetSamplesInCoreset) : TargetSamplesInCoreset(targetSamplesInCoreset),
                                                                                                   NumberOfClusters(numberOfClusters)

{
}

std::shared_ptr<Coreset>
SensitivitySampling::run(const std::vector<dynClustering::CoresetPoint*> &data1, const std::vector<dynClustering::CoresetPoint*> &data2, const std::vector<double> &weights1, const std::vector<double> &weights2)
{
	this->i.init(&data1,&data2);
    std::vector<double> weights = weights1;
    weights.insert(weights.end(),weights2.begin(), weights2.end());

    assert(weights1.size() == data1.size());
    assert(weights2.size() == data2.size());

    clustering::KMeans kMeansAlg(NumberOfClusters, i);

    auto result = kMeansAlg.run(data1, data2, weights);


    auto clusterAssignments = result->getClusterAssignments();

    auto coreset = generateCoresetPoints(data1, data2, clusterAssignments, i, weights);


    return coreset;
}

std::shared_ptr<Coreset>
SensitivitySampling::generateCoresetPoints(const std::vector<dynClustering::CoresetPoint*> &data1, const std::vector<dynClustering::CoresetPoint*> &data2, const clustering::ClusterAssignmentList &clusterAssignments, utils::IndexFinder &i, std::vector<double> &weights)
{

    auto coreset = std::make_shared<Coreset>(*TargetSamplesInCoreset, i);
    

    // Step 2b: compute cost(A). Assume it is the sum of all costs.
    auto sumOfWeightedCosts = clusterAssignments.getTotalWeightedCost();

    // Step 2c: compute the sampling distribution: cost(p, A)/cost(A)
    auto samplingDistribution = clusterAssignments.getWeightedNormalizedCosts();

    auto sampledIndices = random.choice(*TargetSamplesInCoreset, samplingDistribution);

    double T = {0}; //static_cast<double>(*TargetSamplesInCoreset);
    T = static_cast<double>(*TargetSamplesInCoreset);
    

    // Loop through the sampled points and calculate
    // the weight associated with each of these points.
    for (size_t j = 0; j < *TargetSamplesInCoreset; j++)
    {
        size_t sampledPointIndex = (*sampledIndices)[j];

        // We scale the cost of the sampled point by a factor of T i.e. T * cost(p,A)
        double scaledCostPofA = T * clusterAssignments.getPointCost(sampledPointIndex);

        // The weight of the sampled point is now: cost(A) / (T*cost(p,A))
        double weight = sumOfWeightedCosts / scaledCostPofA;

        coreset->addPoint(sampledPointIndex, weight);

        //printf("Sampled point %3ld gets weight %.5f \n", sampledPointIndex, weight);
    }

    auto numberOfClusters = clusterAssignments.getNumberOfClusters();
    auto centerWeights = calcCenterWeights(clusterAssignments, sampledIndices, T, weights);


    for (size_t c = 0; c < numberOfClusters; c++)
    {
        auto weight = (*centerWeights)[c];
        auto center = clusterAssignments.calcCenter(data1, data2, c, i, weights);
        coreset->addCenter(c, center, weight);
    }

    return coreset;
}

std::shared_ptr<blaze::DynamicVector<double>>
SensitivitySampling::calcCenterWeights(
    const clustering::ClusterAssignmentList &clusterAssignments,
    std::shared_ptr<blaze::DynamicVector<size_t>> sampledIndices, double T, std::vector<double> &weights)
{
    auto sumOfWeightedCosts = clusterAssignments.getTotalWeightedCost();
    auto numberOfClusters = clusterAssignments.getNumberOfClusters();

    // Initialise an array to store center weights w_i
    auto centerWeights = std::make_shared<blaze::DynamicVector<double>>(numberOfClusters);
    centerWeights->reset();


    // For each of the T sampled points...
    for (auto &&p : *sampledIndices)
    {
        // Find point p's assigned cluster C_i.
        size_t clusterOfPointP = clusterAssignments.getCluster(p);

        // Find cost(p, A)
        double costPOfA = clusterAssignments.getPointCost(p);

        // Compute cost(A)/(T*cost(p,A))
        double weightContributionOfP = sumOfWeightedCosts / (T * costPOfA);

        // printf("Point %3ld contributes %.5f to cluster %ld  ", p, weightContributionOfP, clusterOfPointP);

        // Sum it up: sum_{p sampled and p in C_i}   cost(A)/(T*cost(p,A))
        (*centerWeights)[clusterOfPointP] += weightContributionOfP;

        // printf("  =>  w_%ld = %.5f\n", clusterOfPointP, (*centerWeights)[clusterOfPointP]);
    }

    // For each of the k' centers, compute the center weights.
    for (size_t c = 0; c < numberOfClusters; c++)
    {
        // Find precomputed center weight: w_i
        double w_i = (*centerWeights)[c];

        // Compute |C_i|
        size_t numberOfWeightedPointsInCluster = clusterAssignments.countWeightedPointsInCluster(c, weights);

        // Compute max(0, |C_i| - w_i)
        double centerWeight = blaze::max(0.0, static_cast<double>(numberOfWeightedPointsInCluster) - w_i);

        // Update the center weight.
        (*centerWeights)[c] = centerWeight;
    }

    return centerWeights;
}
