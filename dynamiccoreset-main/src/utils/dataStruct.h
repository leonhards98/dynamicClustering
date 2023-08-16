#pragma once

#include <cstddef>

namespace dynClustering
{
	/**
	 * @brief Struct containing all parameters for the algorithm
	 */
	struct paramStruct
	{
		/**
		 * Number of clusters in the data
		 */
		const size_t numberOfClusters;

		/**
		 * Parameter epsilon for calculation of the coreset cardinality
		 */
		const double epsilon;

		/**
		 * Parameter lambda for the calculation of the coreset cardinality
		 */
		const double lambda;

		/**
		 * Dimension of the Datapoints
		 */
		const int dimensions;

		/**
		 * Multiplied with "numberOfClusters" to get the number of clusters in the inner coresets
		 */
		const double clusterFactor;

		/**
		 * Multiplied with size of the outer coreset to get the size of the inner coresets
		 */
		const double sFactor;

		/**
		 * Multiplied with cutOff of the outer coreset (same as "s") to get the cutOff of the inner coresets
		 */
		const double cutOffFactor;

		/**
		 * The size the interal Coresets should have
		*/
		int s;

		/**
		 * Nr of new coreset points that should be chosen additinally to the centers of the clusters
		 */
		int NrOfNonCenterCoresetPoints;

		/**
		 * Nr of new coreset points that should be chosen additinally to the centers of the clusters in the inner coreset
		 */
		int NrOfNonCenterCoresetPointsNode;

		/**
		 * Size of the final, outer Coreset
		*/
		int sizeOfFinalCoreset;

		/**
		 * The current cutoff. For nodes that contain more points than the cutoff,
		 * a coreset needs to be calculated
		 */
		double cutOff;

		/**
		 * the total number of points at the last reset of the refresh
		 * (the start of the current phase)
		 */
		int np;

		/**
		 * The number of clusters used in the inner coreset
		 */
		const size_t numberOfClustersNode;

		/**
		 * Factor how much bigger than the original size a coreset is allowed to become by just inserting new points.
		 */
		const double insertionFactor;

		/**
		 * If the coreset of a node contains less than this many points, can just insert new points points into the coreset without
		 * recalculating the coreset
		 */
		double insertionCutOff;

		/**
		 *If a point that is not in the final coreset is removed, do nothing. This is done at most as many times as this attribute until a recalvulation takes place
		 */
		const double removeCutOffFactor;

		double removeCutOff;
		int nrPointsInTree;
		paramStruct(const int nrOfClusters,
					const double epsilon, 
					const double lambda,
					const int dimensions,
					const double clusterFactor,
					const double sFactor,
					const double cutOffFactor, 
					const int sizeOfFinalCoreset,
					const double insertionFactor,
					const double removeCutOffFactor,
					double removeCutOff = 1)
					: numberOfClusters(nrOfClusters), epsilon(epsilon), lambda(lambda), dimensions(dimensions), clusterFactor(clusterFactor), sFactor(sFactor), cutOffFactor(cutOffFactor), numberOfClustersNode(numberOfClusters*clusterFactor), insertionFactor(insertionFactor), removeCutOffFactor(removeCutOffFactor), removeCutOff (removeCutOff) {}
	};
}
