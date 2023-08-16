#include "kmeans.hpp"
#include <float.h>

using namespace clustering;

KMeans::KMeans(size_t k, utils::IndexFinder i) : NumOfClusters(k),i(i)
{
}

std::shared_ptr<ClusteringResult>
KMeans::run(const std::vector<dynClustering::CoresetPoint*> &data1, const std::vector<dynClustering::CoresetPoint*> &data2, const std::vector<double> &weights)
{
  std::vector<size_t> initialCenters;

  auto clusters = this->pickInitialCentersViaKMeansPlusPlus(data1, data2, weights);
  initialCenters = *clusters->getClusterIndices();

  auto centers = copyRows(data1, data2, initialCenters);
  return this->runLloydsAlgorithm(data1, data2, centers, weights);
}

blaze::DynamicMatrix<double>
KMeans::copyRows(const std::vector<dynClustering::CoresetPoint*> &data1, const std::vector<dynClustering::CoresetPoint*> &data2, const std::vector<size_t> &indicesToCopy)
{
  size_t k = indicesToCopy.size();
  size_t d = data1.at(0)->getCoordinates().size();

  blaze::DynamicMatrix<double> centers(k, d);
  for (size_t c = 0; c < k; c++)
  {
    size_t pointIndex = indicesToCopy[c];
    for (int r = 0; r < d; r++)
      centers.at(c,r) = i.getPoint(pointIndex)->getCoordinates().at(r);
  }
  return centers;
}

void computeSquaredNorms(const blaze::DynamicMatrix<double> &dataPoints, std::vector<double> &squaredNorms)
{
  double val = 0.0;
  for (size_t i = 0; i < dataPoints.rows(); i++)
  {
    double sumOfSquares = 0.0;
    for (size_t j = 0; j < dataPoints.columns(); j++)
    {
      val = dataPoints.at(i, j);
      if (val != 0.0)
      {
        sumOfSquares += val * val;
      }
    }
    squaredNorms[i] = sumOfSquares;
  }
}

std::shared_ptr<clustering::ClusterAssignmentList>
KMeans::pickInitialCentersViaKMeansPlusPlus(const std::vector<dynClustering::CoresetPoint*> &data1, const std::vector<dynClustering::CoresetPoint*> &data2, const std::vector<double> &weights)
{
  utils::Random random;
  size_t n = data1.size()+data2.size();
  size_t k = this->NumOfClusters;
  utils::StopWatch sw(true);


  utils::L2NormCalculator squaredL2Norm(data1,data2,i);

  auto clusters = std::make_shared<clustering::ClusterAssignmentList>(n, k);

  // Declare an array used to maintain the squared distances of
  // every point to the closest center among the set of centers
  // that are being considered for any `c \in {2, 3, ..., k}`.
  // Whenever a new center is picked, we can compare the distance
  // between each point `p1` and the new center to determine whether
  // the array at index `p1` should be updated.
  blaze::DynamicVector<double> smallestDistances(n,DBL_MAX);
  blaze::DynamicVector<double> smallestDistancesWeighted(n);

  size_t centerIndex = 0;
  std::vector<size_t> indices;
  for (size_t c = 0; c < k; c++)
  {
    utils::StopWatch pickCenterSW(true);

    if (c == 0)
    {
      // Pick the first centroid uniformly at random.
      centerIndex = random.choice(weights);
      indices.push_back(centerIndex);
    }
    else
    {
      for (size_t p1 = 0; p1 < n; p1++)
      {
        // Compute dist2(p, C_c) i.e., the squared distance between the point `p1` and 
        // the center `centerIndex` that we picked at the previous iteration.
        double distance = squaredL2Norm.calc(p1, centerIndex);

        // Compute min_dist^2(p, C_c-1)
        // Decide if the current distance is better.
        if (distance < smallestDistances[p1] || p1 == centerIndex)
        {
          // Set the weight of a given point to be the smallest distance
          // to any of the previously selected center points. 
          smallestDistances[p1] = distance;
          clusters->assign(p1, centerIndex, distance, distance* weights.at(p1));
        }
      }

      // Pick the index of a point randomly selected based on the distances.
      // A point with a large distance is more likely to be picked than one with
      // a small distance. We want to select points randomly such that points
      // that are far from any of the selected center points have higher likelihood of
      // being picked as the next candidate center.
      for (int r = 0; r < n; r++)
        smallestDistancesWeighted.at(r) = smallestDistances.at(r)*weights.at(r);

      centerIndex = random.choice(smallestDistancesWeighted);
      while (std::find(indices.begin(), indices.end(),centerIndex) != indices.end())
        centerIndex++;

      indices.push_back(centerIndex);

      assert(centerIndex >= 0);
      assert(centerIndex < n);
    }
  }

  // Final reassignment step.
  for (size_t p1 = 0; p1 < n; p1++)
  {
    double distance = squaredL2Norm.calc(p1, centerIndex);
    if (distance < smallestDistances[p1] || p1 == centerIndex)
    {
      clusters->assign(p1, centerIndex, distance, distance*weights.at(p1));
    }
  }


  return clusters;
}

std::shared_ptr<ClusteringResult>
KMeans::runLloydsAlgorithm(const std::vector<dynClustering::CoresetPoint*> &data1, const std::vector<dynClustering::CoresetPoint*> &data2, blaze::DynamicMatrix<double> centroids, const std::vector<double> &weights)
{
  const size_t n = data1.size()+data2.size();
  const size_t d = data1.at(0)->getCoordinates().size();
  const size_t k = this->NumOfClusters;

  blaze::DynamicVector<size_t> clusterMemberCounts(k);
  ClusterAssignmentList cal(n, k);

    cal.assignAll(data1, data2, centroids, i, weights);


  return std::make_shared<ClusteringResult>(cal, centroids);
}
