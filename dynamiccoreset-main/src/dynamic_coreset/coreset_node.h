#pragma once
#include <algorithm>
#include <vector>
#include <memory>
#include <unordered_set>
#include "coreset_point.h"
#include "../coreset/sensitivity_sampling.hpp"
#include "../utils/dataStruct.h"

namespace dynClustering
{
	class CoresetTree;
	/**
	 * @brief Represets a Node in the coreset Tree that contains a list of Points and their weights. This list of points is, except for
	 * leafs, always the union of the points contains in its children. If this union is grater than a cutoff, a coreset is calculated.
	 * The centers of this coreset are also stored in this class
	 */
	class CoresetNode
	{
		/**
		 * A vector containing pointers to the data points that are part of the node. The weights stored in these points are the
		 * original weights and are not the weights the points have in the context of this node
		 */
		std::vector<CoresetPoint *> data;

		/**
		 * The weights that the points have in the context of this node. The actual weight is the weight in this vector
		 * divided by the WeightDenominator. Which weight is connected to which node is determined by
		 * the order in the vector
		 */
		std::vector<double> weights;

		/**
		 * A pointer to the parent node
		 */
		CoresetNode *parent;

		/**
		 * A shared pointer to the left child node
		 */
		std::shared_ptr<CoresetNode> lChild;
		/**
		 * A shared pointer to the right child node
		 */
		std::shared_ptr<CoresetNode> rChild;

		/**
		 * The static coreset class used to calculate the coresets of the unions of the points in the children if
		 * this union is to large
		 */
		coresets::SensitivitySampling staticCoreset;

		/**
		 * A vector containing pointer to the centers that were created by the coreset algorithm. If no coreset is run, this
		 * will be empty. This class owns this points and needs to delete them if it is dectructed or new centers are needed
		 */
		std::vector<CoresetPoint *> centers;

		/**
		 * Struct containing all parameters for the algorithm
		 */
		const paramStruct &parameters;

		/**
		 * Contant used for weights rounding
		 */
		int cForRounding;

		/**
		 * The denomitor of the weights in this node. This means the acutal weights are the ones stored in the weights vector
		 * divided by the weightDenominator
		 */
		int weightDenominator;



		friend CoresetTree;

	public:
		/**
		 * @brief Contructor
		 * @param parameters Struct containing all parameters for the algorithm
		 */
		CoresetNode(const paramStruct &parameters) : parameters (parameters), parent(nullptr), lChild(nullptr), rChild(nullptr), staticCoreset({parameters.numberOfClustersNode, &parameters.NrOfNonCenterCoresetPointsNode})
		{
			this->cForRounding = 0;
			this->resetWeightDenominator();
		}

		/**
		 * @brief Desctructor
		 */
		~CoresetNode()
		{
			// Deletes all centers and resets the shared pointer to the cildren. This destructs the subtree
			for (int i = centers.size() - 1; i >= 0; i--)
				delete centers.at(i);
			lChild.reset();
			rChild.reset();
		}

		/**
		 * @brief Inserts a point into the node. Uses the initial weight of the point. If this is a leaf, set the relevant leaf
		 * of the point to this node
		 * @param point Point to insert
		 * @return
		 */
		bool insertPoint(CoresetPoint *point);

		/**
		 * @brief Collects all points from the children of this node and stores them in the points and weights vector. If there are
		 * to many, calculates a creset first. Then calls the  "runStaticCoreset" function of its parent until the root is hit.
		 * Instantly calls the function of its parent, if this node is a leaf
		 * @param pointToInsert Pointer to the point that is currently inserted. Nullptr if no point gets inserted
		 * @return (void)
		 */
		bool runStaticCoreset(CoresetPoint *pointToInsert, std::unordered_set<CoresetNode*>* newNodesForRemove = nullptr);

		int getLevel()
		{
			int counter {0};
			auto current = this;

			while (current->parent != nullptr)
			{
				current = current->parent;
				counter++;
			}

			return counter;

		}
		bool isLeaf()
		{
			return (lChild == nullptr && rChild == nullptr);
		}
		bool isInternal()
		{
			return !isLeaf();
		}
		bool isRoot()
		{
			return (parent == nullptr);
		}
		const std::shared_ptr<CoresetNode> getLChild()
		{
			return lChild;
		}
		const std::shared_ptr<CoresetNode> getRChild()
		{
			return rChild;
		}
		CoresetNode *getParent()
		{
			return parent;
		}
		void setLChild(std::shared_ptr<CoresetNode> newChild)
		{
			this->lChild = newChild;
		}
		void setRChild(std::shared_ptr<CoresetNode> newChild)
		{
			this->rChild = newChild;
		}
		void setParent(CoresetNode *newParent)
		{
			this->parent = newParent;
		}
		const int getNrPoints()
		{
			return this->data.size();
		}
		const int getWeightsDenominator()
		{
			return this->weightDenominator;
		}

		/**
		 * @brief Sets the weights of the points (weights vector) to their actual weights by dividing them by the weightDenominator.
		 * The sets the weightDenominator to 1
		 * @return (void)
		 */
		void resetWeightDenominator()
		{
			std::for_each(this->weights.begin(), this->weights.end(), [this](double &w)
						  { w /= this->getWeightsDenominator(); });
			this->weightDenominator = 1; // blaze::pow(*np,cForRounding+1)*blaze::ceil(1/eps);
		}
		const std::vector<CoresetPoint *> &getData() { return data; }
		const std::vector<double> &getWeights() { return weights; }
		const std::vector<CoresetPoint *> &getCenters() { return centers; }

		/**
		 * @brief Write the coordinates and weights of all points contained in this node to the given stream
		 * @param out The stream to write to
		 * @return (void)
		 */
		void writePointsToStream(std::ostream &out);

		/**
		 * @brief Write the indices of all points contained in this node to the given stream
		 * @param out The stream to write to
		 * @return (void)
		 */
		void writeIndicesToStream(std::ostream &out);

		/**
		 * @brief Changes the left child of the current node to "newChild" and the parent of the newChild to "newParent"
		 * @param newChild Node that will become the left child of this node
		 * @param newParent Node that will become the parent of "newChild" (most often this current node)
		 * @return (void)
		 */
		void changeLChild(std::shared_ptr<CoresetNode> newChild, CoresetNode *newParent);
		/**
		 * @brief Changes the right child of the current node to "newChild" and the parent of the newChild to "newParent"
		 * @param newChild Node that will become the right child of this node
		 * @param newParent Node that will become the parent of "newChild" (most often this current node)
		 * @return (void)
		 */
		void changeRChild(std::shared_ptr<CoresetNode> newChild, CoresetNode *newParent);

		/**
		 * @brief Removes a point from this node using its index in the points vector. Also removes its weight. The weight will be lost
		 * @param i Index of the point in the points vector to remove. Not the real, unique index of the point.
		 * @return (void)
		 */
		void erasePoint(int i);

		/**
		 * @brief Appends a vector of weights to this node. Does not consider weight denominators
		 * This should only be used to insert Points into leaf nodes.
		 * @param newWeights Vector containing the new weights
		 * @return (void)
		 */
		void appendWeights(const std::vector<double> &newWeights);

		/**
		 * @brief Appends a vector of points to this node. Also changes the relevant leaf of the points to this node.
		 * This should only be used to insert Points into leaf nodes.
		 * @param newPoints
		 * @return (void)
		 */
		void appendPoints(const std::vector<CoresetPoint *> &newPoints);

		/**
		 * @brief Removes all data from the vectors "points", "weights", "centers". Does not delete the points in the vector "centers"
		 * @return (void)
		 */
		void clearData();

		std::shared_ptr<CoresetNode> getSibling()
		{
			assert(this->parent != nullptr);

			CoresetNode *parent = this->parent;

			if (this == parent->lChild.get())
				return parent->rChild;
			else
				return parent->lChild;
		}
	};
}
