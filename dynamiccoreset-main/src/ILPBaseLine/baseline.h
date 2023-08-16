#include "blaze/Math.h"
#include "/opt/gurobi1001/linux64/include/gurobi_c++.h"

/**
 * @brief Calculates an optimal clustering Solution using an ILP. All centers are required to be original data points.
 * Returns the total squared cost of the solution.
 * @param points Points to calculate the clustering for. Each row of the matrix should contain the coordinates of one point
 * @param weights The weights of the points given by the matrix "pooints"
 * @param nrCenters The number of centers that should be used for the calculation
 * @param runtime Reference to a double where the required runtime will be written into
 * @return The total squared cost of the solution
 */
double getBaseLineCostWithCentersAsPoints (const blaze::DynamicMatrix<double>& points, const blaze::DynamicVector<double>& weights, const int nrCenters, const double timeCutOff, double& runtime);

/**
 * @brief Calculates an optimal clustering Solution using an QP. Centers are not required to be original data points.
 * Returns the total squared cost of the solution.
 * @param points Points to calculate the clustering for. Each row of the matrix should contain the coordinates of one point
 * @param weights The weights of the points given by the matrix "pooints"
 * @param nrCenters The number of centers that should be used for the calculation
 * @param runtime Reference to a double where the required runtime will be written into
 * @return The total squared cost of the solution
 */
double getBaseLineCostWithArbitraryCenters (blaze::DynamicMatrix<double>& points, const blaze::DynamicVector<double>& weights, const int nrCenters, const std::shared_ptr<std::vector<size_t>> centerHints, const double timeCutOff, double& runtime);
