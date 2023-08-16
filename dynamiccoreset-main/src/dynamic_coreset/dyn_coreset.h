#pragma once

#include <ostream>
#include <vector>
#include <memory>
#include <boost/unordered/unordered_map.hpp>
#include "coreset_point.h"
#include "coreset_node.h"
#include "coreset_tree.h"
#include "../utils/dataStruct.h"

namespace dynClustering
{
	/**
	 * @brief Maintains a tree where the leafs contain datapoints and the internal nodes
	 * are coresets of the points of their children. Uses this tree to dynamically
	 * maintain a coreset of a given size.
	 */
	class DynCoreset
	{
		/**
		 * Lef node of the tree where current insertions take place
		 */
		std::shared_ptr<CoresetNode> currentNode;

		/**
		 * All points currently contained in the tree, with their indices as keys
		 */
		boost::unordered_map<int, CoresetPoint> data;

		CoresetTree *tree;
		coresets::SensitivitySampling staticCoreset;

		paramStruct &parameters;


		/**
		 * Function that returns the size the coreset must have to gurantee
		 * the given values for lambda and epsilon
		 */
		std::function<int(double epsilon, double lambda, int k, int d)> cardinalityOfCoreset;


		/**
		 * Number of updates since the last remove operation was done where the points were
		 * actually removed from the coreset tree
		 */
		int updatesSinceLastRemoveRecalc {0};

		/**
		 * Contains all the iterators for the hash table to the points that should have been removed already but weren't as part of the remove shortcut
		 */
		//std::vector<boost::unordered_map<int,CoresetPoint>::iterator> PointsToRemoveIter;
		/**
		 * Contains all the indices of the  points that should have been removed already but weren't as part of the remove shortcut
		 */
		std::vector<CoresetPoint*> PointsToRemoveIter;

		int updatesThisPhase;

		/**
		 * empty vector
		 */
		std::vector<CoresetPoint *> data2Proxy{};

		/**
		 * empty vector
		 */
		std::vector<double> weights2Proxy{};

		/**
		 * Contains the centers calculated by the outer coreset. These poins need to be deleted when they are not needed
		 */
		std::vector<CoresetPoint *> centers{};

		/**
		 * Contains all points contained in the outer coreset, including centers.
		 */
		std::vector<CoresetPoint *> coresetPoints{};

		/**
		 * Contains the weights of all points contained in the outer coreset
		 */
		std::vector<double> coresetWeights{};

		/**
		 * @brief Inserts a point into the tree and updates all data structures accoringly. Calculates an updated final coreset
		 * @param index index of the point to be inserted
		 * @param weight initial weight of the point to be inserted
		 * @param coordinates coordinates of the point to be inserted
		 * @return
		 */
		bool insertPoint(int index, double weight, std::vector<double> &coordinates);

		/**
		 * @brief Inserts a point into the tree and updates all data structures accoringly. Calculates an updated final coreset
		 * @param point Point to be inserted as a matrix with one row that
		 * contains the index, if the point in iserted (1) or deleted (0), the weight, the coordinates
		 * @return
		 */
		bool insertPoint(blaze::DynamicMatrix<double> &point);

		/**
		 * @brief Remove a point from the tree and updates all data structures accoringly. Calculates an updated final coreset
		 * @param index Index of the point to remove
		 * @return
		 */
		bool removePoint(int index);

		/**
		 * @brief Reset the refresh procedure. Change n0 to the current number of points, move the refresh pointer to the starting position
		 * @param nrPointsIncreased set true if the last operation increased the number of points in the datastructure
		 * @return (void)
		 */
		void goToNextPhase(bool nrPointsIncreased);

		/**
		 * @brief run the static coreset algorithm on the points on the root of the tree to get the final result (outer coreset)
		 * @return (void)
		 */
		void runOuterCoreset();

	public:
		/**
		 * @brief Writes all points and centers containen in the final, outer coreset to the given stream
		 * @param out Stream to write the coreset to
		 * @return (void)
		 */
		void WriteCoresetToStream(std::ostream &out);

		/**
		 * @brief Constructor
		 * @param parameters Struct containing all parameters dor the algorithm
		 * @param cardinalityOfCoreset Function used to get the size of the final Coreset
		 */
		DynCoreset(paramStruct &parameters, std::function<int(double epsilon, double lambda, int k, int d)> cardinalityOfCoreset) : parameters(parameters), staticCoreset({parameters.numberOfClusters,&parameters.NrOfNonCenterCoresetPoints})
		{
			this->updatesThisPhase = 0;
			this->parameters.sizeOfFinalCoreset = cardinalityOfCoreset(parameters.epsilon / 3, parameters.lambda / 2, parameters.numberOfClusters, parameters.dimensions);
			this->cardinalityOfCoreset = cardinalityOfCoreset;
		}
		/**
		 * @brief Destructor. Deletes all heap allocated memory
		 */
		~DynCoreset()
		{
			delete tree;
			for (int i = centers.size() - 1; i >= 0; i--)
				delete centers.at(i);
		}

		/**
		 * @brief Write the structure of the tree to the given stream
		 * @param out stream to write to
		 * @return (void)
		 */
		void writeTreeToStream(std::ostream &out)
		{
			out << "Final Coreset: " << this->coresetPoints.size() << "\n";
			tree->writeToStream(out);
		}

		/**
		 * @brief Constructs the tree as a root and two children and fills it with the given points.
		 * The number of points needs to be larger than s. These are evenly split on the two leafs.
		 * @param firstSPoints coordinates of the points to initialize the tree with as a Matrix
		 * Each row is one point. First column is the index, second 1 if the point gets inserted or 0 if deleted, third row the weight followed by the coordinates
		 * @return
		 */
		bool initTree(std::shared_ptr<blaze::DynamicMatrix<double>> firstSPoints);

		/**
		 * @brief Apply an Operation: Inserts or Remove a point from the tree and updates all data structures accoringly. Calculates an updated final coreset
		 * @param point Point of the operation as a matrix with one row that
		 * contains the index, the weight, if the point in iserted (1) or deleted (0), the coordinates. If the point is removed,
		 * the cooridinates will not be used
		 * @return
		 */
		bool applyOperation(blaze::DynamicMatrix<double> &point);

		const std::vector<CoresetPoint*>& getCoresetPoints() const
		{
			return this->coresetPoints;
		}
		const std::vector<double>& getCoresetWeights() const
		{
			return this->coresetWeights;
		}
		const int getCoresetSize() const
		{
			return coresetPoints.size();
		}
	};
}
