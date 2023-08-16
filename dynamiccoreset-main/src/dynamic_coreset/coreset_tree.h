#pragma once

#include "boost/unordered/unordered_map_fwd.hpp"
#include "coreset_point.h"
#include "coreset_node.h"
#include <memory>
#include <ostream>
#include <functional>
#include <vector>
#include "../utils/dataStruct.h"
#include <boost/unordered/unordered_map.hpp>


namespace dynClustering
{
	/**
	 * @brief A tree where the leafs contain datapoints and the internal nodes
	 * are coresets of the points of their children. This tree can be dynamically
	 * maintained
	 */
	class CoresetTree
	{
		std::shared_ptr<CoresetNode> root;

		/**
		 * Leaf where new points are currently inserted
		 */
		CoresetNode *currentLeaf;

		/**
		 * Struct containing all parameters for the algorithm
		 */
		paramStruct &parameters;


		/**
		 * Points to the leaf that needs to be refreshed in the next step
		 */
		CoresetNode *refreshPointer;

		/**
		 * Points to a leaf that is to small in the current phase due to points being removed.
		 * Will be combined with other leafs at refreshed
		 */
		CoresetNode *toSmallLeaf;

		/**
		 * @brief Inserts a leaf into the tree. Keeps the tree a balances binary tree
		 * Leaf is inserted at the rightmost position in the lowest level
		 * @return
		 */
		std::shared_ptr<CoresetNode> insertLeaf();

		/**
		 * @brief removed a given leaf. Swaps the leaf with the rightmost leaf on the lowest level
		 * and then removed that leaf. All data in the leaf ist lost.
		 * @param LeafToRemove
		 * @return (void)
		 */
		void removeLeaf(std::shared_ptr<CoresetNode> LeafToRemove);

		/**
		 * @brief Inserts a leaf into the tree. The given parent leaf becomes the
		 * sibling of the new leaf and a new internal node is inserted as a parent of
		 * both the new leaf and the given parent
		 * @param parent Leaf where the new leaf is inserted. This leaf will be the sibling
		 * of the new leaf with a new parent
		 * @return
		 */
		std::shared_ptr<CoresetNode> insertChild(std::shared_ptr<CoresetNode> &parent);

		/**
		 * @brief Recalculates all coresets from the leaf pointed to by the refresh pointer
		 * to the root. Does not run the outer coreset. Then moves the refresh pointer to the right
		 * @param onlyMove if true, only moves the refresh pointer without calculating the coresets
		 * @return (void)
		 */
		void doRefresh(bool onlyMove = false);

		/**
		 * @brief Move all points and weights from the first to the second leaf. Also considers weight denomitors
		 * @param from
		 * @param to
		 * @return
		 */
		bool movePointsBetweenLeafs(CoresetNode *from, CoresetNode *to);

		/**
		 * @brief Swaps the position of the two given leafs in the tree
		 * @param leaf1
		 * @param leaf2
		 * @return
		 */
		bool swapLeafs(std::shared_ptr<CoresetNode> leaf1, std::shared_ptr<CoresetNode> leaf2);

		/**
		 * @brief Returns the rightmost leaf on the lowest level
		 * @return Rightmost leaf
		 */
		std::shared_ptr<CoresetNode> getRightmostLeaf();

		/**
		 * @brief Returns the leftmost leaf on the lowest incomplete level.
		 * If all levels are complete, returns the leftmost leaf on the lowest level
		 * @return Leftmost leaf
		 */
		std::shared_ptr<CoresetNode> getLeftmostLeaf();

		/**
		 * @brief Returns the leaf next to the given leaf on the right side
		 * @param leaf
		 * @return Leaf to the right
		 */
		CoresetNode *getLeafToRight(CoresetNode *leaf);

		/**
		 * @brief Checks if any leaf starting from the leaf pointed to by the refresh pointer is smaller that allowed,
		 * which means it has less points than s/2. Marks to small leafs by pointing the pointer
		 * "toSmallLef" to it and combines them. If the pointer toSmallLeaf is not a nullpointer after this function returns,
		 * the number of leafs that are to small was odd.
		 * @return (void)
		 */
		void checkIfLeafsToSmall();

		/**
		 * @brief Checks if any leaf starting from the leaf pointed to by the refresh pointer is larger than allowed,
		 * which means it has more than s points. To large leafs are split into two leafs.
		 * @return (void)
		 */
		void checkIfLeafsToLarge();

	public:
		/**
		 * @brief Constructor
		 * @param parameters Struct containing all parameters for the algorithm
		 */
		CoresetTree(paramStruct &parameters) : parameters(parameters)
		{

			std::shared_ptr<CoresetNode> lChild = std::make_shared<CoresetNode>(parameters);
			std::shared_ptr<CoresetNode> rChild = std::make_shared<CoresetNode>(parameters);
			this->root = std::make_shared<CoresetNode>(parameters);
			root->changeLChild(lChild, root.get());
			root->changeRChild(rChild, root.get());

			this->refreshPointer = lChild.get();
			this->toSmallLeaf = nullptr;
		}
		bool setCurrentLeaf(CoresetNode *leaf)
		{
			if (leaf == nullptr)
			{
				this->currentLeaf = leaf;
				return 0;
			}
			if (!leaf->isLeaf())
				return 1;
			this->currentLeaf = leaf;
			return 0;
		}

		/**
		 * @brief Insert a new Point the leaf pointed to by the "currentLeaf" pointer.
		 * If this leaf has to many points, insert a new one
		 * @param newPoint Point to insert
		 * @return True if the insert shortcut was used
		 */
		bool insertPoint(CoresetPoint *newPoint);

		/**
		 * @brief Remove a point. If the number of Points in the leaf becomes to small,
		 * add the remaining points to the leaf pointed to by "currentLeaf", then delete the leaf
		 * @param point Point to delete
		 * @return
		 */
		bool removePoint(CoresetPoint *point);

		/**
		 * @brief Removes points. If the number of Points in the leaf becomes to small,
		 * add the remaining points to the leaf pointed to by "currentLeaf", then delete the leaf
		 * @param point Point to delete
		 * @return
		 */
		//bool removePoints(std::vector<boost::unordered_map<int,CoresetPoint>::iterator> points);
		bool removePoints(std::vector<CoresetPoint*> points);

		/**
		 * @brief Reset the refresh procedure. Change n0 to the current number of points, move the refresh pointer to the
		 *  starting position. Then checks if any leafs are to small or to large.
		 * @param nrPointsIncreased set true if the last operation increased the number of points in the datastructure
		 * @return (void)
		 */
		void resetRefresh(bool nrPointsIncreased);

		/**
		 * @brief Write the structure of the tree to the given stream
		 * @param out Strem to write to
		 * @return (void)
		 */
		void writeToStream(std::ostream &out);

		/**
		 * @brief Write the coreset of the root node to the given stream
		 * @param out the stream to write to
		 * @return
		 */
		void writeCoresetToStream(std::ostream &out) { this->root->writePointsToStream(out); }

		std::shared_ptr<CoresetNode> getRoot() { return this->root; }

		int getTreeHeight()
		{
			int counter {0};
			auto current = this->root;

			while (current->isLeaf() == false)
			{
				current = current->getLChild();
				counter++;
			}
			return counter;

		}
	};

}
