#include "cluster_assignment_list.hpp"

using namespace clustering;

ClusterAssignmentList::ClusterAssignmentList(size_t n, size_t k) : numOfPoints(n), numOfClusters(k), clusters(n), distances(n), weightedDistances(n)
{
}

ClusterAssignmentList::ClusterAssignmentList(const ClusterAssignmentList &other) : numOfPoints(other.numOfPoints), numOfClusters(other.numOfClusters), clusters(other.clusters), distances(other.distances), weightedDistances(other.weightedDistances)
{
}

void ClusterAssignmentList::assign(size_t pointIndex, size_t clusterIndex, double distance, double weightedDistance)
{
    // TODO: Ensure arguments are not out of range to avoid runtime errors.
    clusters[pointIndex] = clusterIndex;
    distances[pointIndex] = distance;
    weightedDistances[pointIndex] = weightedDistance;
}

void ClusterAssignmentList::assignAll(const std::vector<dynClustering::CoresetPoint*> &data1, const std::vector<dynClustering::CoresetPoint*> &data2, const blaze::DynamicMatrix<double> &centers, utils::IndexFinder &i, const std::vector<double> &weights)
{
    auto n = this->numOfPoints;
    auto k = this->numOfClusters;

    // For each data point, assign the centroid that is closest to it.
    for (size_t p = 0; p < n; p++)
    {
        double bestDistance = std::numeric_limits<double>::max();
        size_t bestCluster = 0;

        // Loop through all the clusters.
        for (size_t c = 0; c < k; c++)
        {
            // Compute the L2 norm between point p and centroid c.
            //const double distance = blaze::norm(blaze::row(dataPoints, p) - blaze::row(centers, c));
            const double distance = blaze::norm(i.getPoint(p)->getCoordinates() - blaze::trans(blaze::row(centers, c)));

            // Decide if current distance is better.
            if (distance < bestDistance)
            {
                bestDistance = distance;
                bestCluster = c;
            }
        }
        

        // Assign cluster to the point p.
        this->assign(p, bestCluster, bestDistance, bestDistance*weights.at(p));
    }
}

size_t
ClusterAssignmentList::getCluster(size_t pointIndex) const
{
    // TODO: Ensure arguments are not out of range to avoid runtime errors.
    return clusters[pointIndex];
}

size_t
ClusterAssignmentList::getNumberOfPoints() const
{
    return this->numOfPoints;
}

size_t
ClusterAssignmentList::getNumberOfClusters() const
{
    return this->numOfClusters;
}

blaze::DynamicVector<double> &
ClusterAssignmentList::getCentroidDistances()
{
    return this->distances;
}

size_t
ClusterAssignmentList::countPointsInCluster(size_t clusterIndex) const
{
    size_t count = 0;
    for (size_t p = 0; p < this->numOfPoints; p++)
    {
        if (clusters[p] == clusterIndex)
        {
            count++;
        }
    }

    return count;
}

size_t
ClusterAssignmentList::countWeightedPointsInCluster(size_t clusterIndex, std::vector<double> &weights) const
{
    size_t count = 0;
    for (size_t p = 0; p < this->numOfPoints; p++)
    {
        if (clusters[p] == clusterIndex)
        {
            count += weights.at(p);
        }
    }

    return count;
}

double
ClusterAssignmentList::getTotalCost() const
{
    return blaze::sum(this->distances);
}

double
ClusterAssignmentList::getTotalWeightedCost() const
{
    return blaze::sum(this->weightedDistances);
}

double
ClusterAssignmentList::getPointCost(size_t pointIndex) const
{
    // TODO: Ensure pointIndex is not out of bounds.

    return this->distances[pointIndex];
}

double
ClusterAssignmentList::getWeightedPointCost(size_t pointIndex) const
{
    // TODO: Ensure pointIndex is not out of bounds.

    return this->weightedDistances[pointIndex];
}

std::shared_ptr<blaze::DynamicVector<double>>
ClusterAssignmentList::calcAverageClusterCosts() const
{
    auto results = std::make_shared<blaze::DynamicVector<double>>(this->numOfClusters);
    results->reset();

    blaze::DynamicVector<double> counts(this->numOfClusters);

    for (size_t p = 0; p < this->numOfPoints; p++)
    {
        auto c = clusters[p];
        (*results)[c] += distances[p];
        counts[c] += 1;
    }

    for (size_t c = 0; c < this->numOfClusters; c++)
    {
        (*results)[c] /= counts[c];
    }

    return results;
}

std::shared_ptr<blaze::DynamicVector<double>>
ClusterAssignmentList::calcClusterCosts() const
{
    auto results = std::make_shared<blaze::DynamicVector<double>>(this->numOfClusters);
    results->reset();

    for (size_t p = 0; p < this->numOfPoints; p++)
    {
        auto c = clusters[p];
        (*results)[c] += distances[p];
    }

    return results;
}

ClusterAssignmentList &
ClusterAssignmentList::operator=(const ClusterAssignmentList &other)
{
    this->numOfPoints = other.numOfPoints;
    this->numOfClusters = other.numOfClusters;
    this->clusters = other.clusters;
    this->distances = other.distances;
    return *this;
}

blaze::DynamicVector<double>
ClusterAssignmentList::getNormalizedCosts() const
{
    return distances / this->getTotalCost();
}

blaze::DynamicVector<double>
ClusterAssignmentList::getWeightedNormalizedCosts() const
{
    return weightedDistances / this->getTotalWeightedCost();
}


size_t
getIndex(const std::vector<size_t> &values, const size_t valueToSearch)
{
    auto it = std::find(values.begin(), values.end(), valueToSearch);
    if (it != values.end())
    {
        auto index = it - values.begin();
        if (index > 0)
        {
            return static_cast<size_t>(index);
        }
    }

    throw std::logic_error("Value not found.");
}

void ClusterAssignmentList::calcCenters(const std::vector<dynClustering::CoresetPoint*> &data1,const std::vector<dynClustering::CoresetPoint*> &data2, blaze::DynamicMatrix<double> &newCenters, utils::IndexFinder &i)
{
    assert(data1.size()+data2.size() == this->numOfPoints);

    const size_t n = data1.size()+data2.size();
    const size_t d = data1.at(0)->getCoordinates().size();
    const size_t k = this->numOfClusters;

    if (newCenters.rows() != k or newCenters.columns() != d)
    {
        newCenters.resize(n, d);
    }

    // Reset variables.
    newCenters = 0;
    std::vector<size_t> clusterMemberCounts(k);
    for (size_t c = 0; c < k; c++)
    {
        clusterMemberCounts[c] = 0;
    }

    auto centerIndices = *this->getClusterIndices();

    for (size_t p = 0; p < n; p++)
    {
        const size_t c = getIndex(centerIndices, this->clusters[p]);

        // Sum all points assigned in cluster c
		for (int r = 0; r < d; r++)
			newCenters(c,r) += i.getPoint(p)->getCoordinates().at(r);

        // Count number of points in each cluster.
        clusterMemberCounts[c] += 1;
    }

    // Divide centers by the number of points in each cluster.
    for (size_t c = 0; c < k; c++)
    {
        const auto count = std::max<size_t>(1, clusterMemberCounts[c]);
        blaze::row(newCenters, c) /= count;
    }
}

std::shared_ptr<blaze::DynamicVector<double>>
ClusterAssignmentList::calcCenter(const std::vector<dynClustering::CoresetPoint*> &data1,const std::vector<dynClustering::CoresetPoint*> &data2, const size_t clusterIndex, utils::IndexFinder &i, std::vector<double> &weights) const
{
    const size_t n = data1.size()+data2.size();
    const size_t d = data1.at(0)->getCoordinates().size();

    auto center = std::make_shared<blaze::DynamicVector<double>>();
    center->resize(d);

    // Reset variables
    center->reset();
    double memberCount = 0.0;

    for (size_t p = 0; p < n; p++)
    {
        const size_t c = this->clusters[p];

        if (c == clusterIndex)
        {
            for (size_t j = 0; j < d; j++)
            {
                center->at(j) += (i.getPoint(p)->getCoordinates().at(j)*weights.at(p));
            }

            memberCount += weights.at(p);
        }
    }

    // Divide centers by the number of points in each cluster.
    if (memberCount == 0)
    {
        for (size_t j = 0; j < d; j++)
        {
            center->at(j) = i.getPoint(clusterIndex)->getCoordinates().at(j);
        }
        return center;
    }
    for (size_t j = 0; j < d; j++)
    {
        center->at(j) /= memberCount;
    }

    return center;
}

std::shared_ptr<std::vector<size_t>>
ClusterAssignmentList::getPointsByCluster(size_t clusterIndex) const
{
    auto pointsToReturn = std::make_shared<std::vector<size_t>>();
    for (size_t p = 0; p < this->numOfPoints; p++)
    {
        if (clusters[p] == clusterIndex)
        {
            pointsToReturn->push_back(p);
        }
    }
    return pointsToReturn;
}

std::shared_ptr<std::vector<size_t>>
ClusterAssignmentList::getClusterIndices() const
{
    auto clusterIndexSet = std::set<size_t>();
    for (size_t p = 0; p < this->numOfPoints; p++)
    {
        clusterIndexSet.insert(clusters[p]);
    }

    // Convert set to a vector.
    auto clusterIndices = std::make_shared<std::vector<size_t>>(clusterIndexSet.begin(), clusterIndexSet.end());
    assert(clusterIndices->size() == this->numOfClusters);
    return clusterIndices;
}
