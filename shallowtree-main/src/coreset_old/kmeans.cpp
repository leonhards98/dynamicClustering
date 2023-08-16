#include "kmeans.hpp"
#include <float.h>

using namespace clusteringStatic;
using namespace utilsStatic;


KMeans::KMeans(size_t k, bool useKmeansPlusPlus, size_t nIter, double convDiff) : NumOfClusters(k), InitKMeansPlusPlus(useKmeansPlusPlus), MaxIterations(nIter), ConvergenceDiff(convDiff)
{
}

std::shared_ptr<ClusteringResult>
KMeans::run(const blaze::DynamicMatrix<double> &data, const blaze::DynamicVector<double> &weights, bool newKMeansPP)
{
  std::vector<size_t> initialCenters;

  assert(data.rows() == weights.size());


  std::shared_ptr<ClusterAssignmentList> clusters;
  if (newKMeansPP)
  {
    fastKMeans::fastKMeansAlgo kMeansPP(data, weights);
    clusters = kMeansPP.fastKMeansPP(this->NumOfClusters);
  }
  else
    clusters = this->pickInitialCentersViaKMeansPlusPlus(data, weights);

  initialCenters = *clusters->getClusterIndices();

  auto centers = copyRows(data, initialCenters);
  return this->runLloydsAlgorithm(data, centers, weights);
}

blaze::DynamicMatrix<double>
KMeans::copyRows(const blaze::DynamicMatrix<double> &data, const std::vector<size_t> &indicesToCopy)
{
  size_t k = indicesToCopy.size();
  size_t d = data.columns();

  blaze::DynamicMatrix<double> centers(k, d);
  for (size_t c = 0; c < k; c++)
  {
    size_t pointIndex = indicesToCopy[c];
    blaze::row(centers, c) = blaze::row(data, pointIndex);
  }
  return centers;
}
namespace clusteringStatic{
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
}

std::shared_ptr<clusteringStatic::ClusterAssignmentList>
KMeans::pickInitialCentersViaKMeansPlusPlus(const blaze::DynamicMatrix<double> &data, const blaze::DynamicVector<double> &weights)
{
  utilsStatic::Random random;
  size_t n = data.rows();
  size_t k = this->NumOfClusters;
  utilsStatic::StopWatch sw(true);

  utilsStatic::L2NormCalculator squaredL2Norm(data, true);

  auto clusters = std::make_shared<clusteringStatic::ClusterAssignmentList>(n, k);

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
    utilsStatic::StopWatch pickCenterSW(true);

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
          smallestDistancesWeighted[p1] = distance * weights.at(p1);
          clusters->assign(p1, centerIndex, distance, distance * weights.at(p1));
        }
      }

      // Pick the index of a point randomly selected based on the distances.
      // A point with a large distance is more likely to be picked than one with
      // a small distance. We want to select points randomly such that points
      // that are far from any of the selected center points have higher likelihood of
      // being picked as the next candidate center.

      centerIndex = random.choice(smallestDistancesWeighted);
      if (std::find(indices.begin(), indices.end(),centerIndex) != indices.end())
        centerIndex++;

      indices.push_back(centerIndex);
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
KMeans::runLloydsAlgorithm(const blaze::DynamicMatrix<double> &matrix, blaze::DynamicMatrix<double> centroids, const blaze::DynamicVector<double> &weights)
{
  const size_t n = matrix.rows();
  const size_t d = matrix.columns();
  const size_t k = this->NumOfClusters;

  blaze::DynamicVector<double> clusterMemberCounts(k);
  ClusterAssignmentList cal(n, k);
  double totalCostSquared {0};


  if (MaxIterations == 0)
  {
    cal.assignAll(matrix, centroids, weights);
  }

  if (MaxIterations > 0)
  {
    std::vector<double> dataSquaredNorms;
    dataSquaredNorms.resize(n);
    clusteringStatic::computeSquaredNorms(matrix, dataSquaredNorms);

    std::vector<double> centerSquaredNorms;
    centerSquaredNorms.resize(centroids.rows());
    clusteringStatic::computeSquaredNorms(centroids, centerSquaredNorms);

    // Lambda function computes the squared L2 distance between any pair of points.
    // The function will automatically use any precomputed distance if it exists.
    auto calcSquaredL2Norm = [&matrix, &centroids, d, &dataSquaredNorms, &centerSquaredNorms](size_t p, size_t c) -> double
    {
      double dotProd = 0.0;
      for (size_t i = 0; i < d; i++)
      {
        dotProd += matrix.at(p, i) * centroids.at(c, i);
      }

      return dataSquaredNorms[p] + centerSquaredNorms[c] - 2 * dotProd;
    };


    for (size_t i = 0; i < this->MaxIterations; i++)
    {
      utilsStatic::StopWatch iterSW(true);
      // For each data point, assign the centroid that is closest to it.
      for (size_t p = 0; p < n; p++)
      {
        double bestDistance = std::numeric_limits<double>::max();
        size_t bestCluster = 0;

        // Loop through all the clusters.
        for (size_t c = 0; c < k; c++)
        {
          // Compute the L2 norm between point p and centroid c.
          // const double distance = blaze::norm(blaze::row(matrix, p) - blaze::row(centroids, c));
          const double distance = calcSquaredL2Norm(p, c);

          // Decide if current distance is better.
          if (distance < bestDistance)
          {
            bestDistance = distance;
            bestCluster = c;
          }
        }

        // Assign cluster to the point p.
        const double dist = blaze::sqrt(bestDistance);
        cal.assign(p, bestCluster, dist, dist * weights.at(p));
      }

      // Move centroids based on the cluster assignments.

      // First, save a copy of the centroids matrix.
      blaze::DynamicMatrix<double> oldCentrioids(centroids);

      // Set all elements to zero.
      centroids = 0;           // Reset centroids.
      clusterMemberCounts = 0; // Reset cluster member counts.

      for (size_t p = 0; p < n; p++)
      {
        //double x = blaze::row(matrix, p).at(0);
        //double y = blaze::row(matrix, p).at(1);
        //double w = weights.at(p);
        //double a = (blaze::row(matrix, p)*weights.at(p)).at(0);
        //double b = (blaze::row(matrix, p)*weights.at(p)).at(1);
        const size_t c = cal.getCluster(p);
        //double u = blaze::row(centroids, c).at(0);
        //double v = blaze::row(centroids, c).at(1);
        blaze::row(centroids, c) += blaze::row(matrix, p)*weights.at(p);
        //double s = blaze::row(centroids, c).at(0);
        //double t = blaze::row(centroids, c).at(1);
        clusterMemberCounts[c] += weights.at(p);
      }

    //std::cout << centroids << "\n";
      for (size_t c = 0; c < k; c++)
      {
        double f = clusterMemberCounts[c];
        const auto count = std::max<double>(1, clusterMemberCounts[c]);
        blaze::row(centroids, c) /= count;
      }
    //std::cout << centroids << "\n";

      // Recompute the squared distances again.
      clusteringStatic::computeSquaredNorms(centroids, centerSquaredNorms);

      // Compute the Frobenius norm
      auto diffAbsMatrix = blaze::abs(centroids - oldCentrioids);
      auto diffAbsSquaredMatrix = blaze::pow(diffAbsMatrix, 2); // Square each element.
      auto frobeniusNormDiff = blaze::sqrt(blaze::sum(diffAbsSquaredMatrix));

      if (frobeniusNormDiff < this->ConvergenceDiff || i == this->MaxIterations-1)
      {
        totalCostSquared = 0;
        for (size_t p = 0; p < n; p++)
        {
          const size_t c = cal.getCluster(p);
          totalCostSquared += calcSquaredL2Norm(p,c)*weights.at(p);
        }

        break;
      }
    }
  }

  return std::make_shared<ClusteringResult>(cal, centroids, totalCostSquared);
}
