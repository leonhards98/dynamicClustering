#include "clustering_result.hpp"

using namespace clusteringStatic;

ClusteringResult::ClusteringResult(const ClusterAssignmentList &assignments, blaze::DynamicMatrix<double> &finalCentroids, double totalSquaredCost) :
    clusterAssignments(assignments), centroids(finalCentroids), totalSquaredCost(totalSquaredCost)
{
}

ClusterAssignmentList&
ClusteringResult::getClusterAssignments()
{
    return this->clusterAssignments;
}

blaze::DynamicMatrix<double>&
ClusteringResult::getCentroids()
{
    return this->centroids;
}

double
ClusteringResult::getTotalSquaredCost()
{
    return this->totalSquaredCost;
}
