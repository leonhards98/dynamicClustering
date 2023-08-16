#include "coreset_tree.h"
#include <memory>
#include <queue>
#include <unordered_map>
#include <vector>
#include <unordered_set>

namespace dynClustering
{

	std::shared_ptr<CoresetNode> CoresetTree::insertLeaf()
	{
		std::shared_ptr<CoresetNode> leftmostLeaf = this->getLeftmostLeaf();

		return this->insertChild(leftmostLeaf);
	}

	std::shared_ptr<CoresetNode> CoresetTree::getRightmostLeaf()
	{
		std::queue<std::shared_ptr<CoresetNode>> queue;

		queue.push(this->root);

		std::shared_ptr<CoresetNode> rightmostLeaf = nullptr;
		while (!queue.empty())
		{
			rightmostLeaf = queue.front();
			queue.pop();
			if (rightmostLeaf->isLeaf())
				continue;
			queue.push(rightmostLeaf->getLChild());
			queue.push(rightmostLeaf->getRChild());
		}

		return rightmostLeaf;
	}

	std::shared_ptr<CoresetNode> CoresetTree::getLeftmostLeaf()
	{
		std::queue<std::shared_ptr<CoresetNode>> queue;

		queue.push(this->root);

		std::shared_ptr<CoresetNode> leftmostLeaf = nullptr;
		while (!queue.empty())
		{
			leftmostLeaf = queue.front();
			queue.pop();
			if (leftmostLeaf->isLeaf())
				break;
			queue.push(leftmostLeaf->getLChild());
			queue.push(leftmostLeaf->getRChild());
		}

		return leftmostLeaf;
	}

	void CoresetTree::removeLeaf(std::shared_ptr<CoresetNode> LeafToRemove)
	{
		std::shared_ptr<CoresetNode> rightmostLeaf = this->getRightmostLeaf();

		CoresetNode *parent = rightmostLeaf->getParent();
		CoresetNode *grandParent = parent->getParent();
		assert(grandParent != nullptr);
		std::shared_ptr<CoresetNode> sibling = parent->getLChild();

		// If the leaf to remove is the sibling of the rightmost leaf, just remove
		// the leaf and swap the rightmost leaf with its parent. Then remove the parent
		if (sibling.get() == LeafToRemove.get())
		{
			parent->lChild.reset();
			parent->rChild.reset();

			if (parent == grandParent->lChild.get())
			{
				grandParent->lChild.reset();
				grandParent->lChild = rightmostLeaf;
			}
			else
			{
				grandParent->rChild.reset();
				grandParent->rChild = rightmostLeaf;
			}
			rightmostLeaf->parent = grandParent;

			rightmostLeaf->runStaticCoreset(nullptr);

			return;
		}
		if (rightmostLeaf.get() != LeafToRemove.get())
			this->swapLeafs(LeafToRemove, rightmostLeaf);

		// Remove the rightmost leaf (the leaf to remove. swapped if neccersary)
		// and swap its sibling with its parent. Then remove the parent.
		if (parent == grandParent->lChild.get())
		{
			grandParent->lChild.reset();
			grandParent->lChild = sibling;
		}
		else
		{
			grandParent->rChild.reset();
			grandParent->rChild = sibling;
		}
		sibling->parent = grandParent;

		sibling->runStaticCoreset(nullptr);
		if (rightmostLeaf != LeafToRemove)
			rightmostLeaf->runStaticCoreset(nullptr);
		if (LeafToRemove.get() == this->refreshPointer)
			refreshPointer = sibling.get();

		return;
	}

	std::shared_ptr<CoresetNode> CoresetTree::insertChild(std::shared_ptr<CoresetNode> &parent)
	{
		// Create a new node and make it a child of the parent node. Then swap the
		// parent and the new child. Insert a new leaf as a sibling to the original parent
		// and a child to the new parent.
		// Return this new leaf
		CoresetNode *grandParent = parent->getParent();
		parent->resetWeightDenominator();
		std::shared_ptr<CoresetNode> newInternalNode = std::make_shared<CoresetNode>(parameters);
		std::shared_ptr<CoresetNode> newLeafNode = std::make_shared<CoresetNode>(parameters);
		if (grandParent->getLChild()->isLeaf())
		{
			grandParent->changeLChild(newInternalNode, grandParent);
			newInternalNode->changeLChild(parent, newInternalNode.get());
			newInternalNode->changeRChild(newLeafNode, newInternalNode.get());
		}
		else
		{
			grandParent->changeRChild(newInternalNode, grandParent);
			newInternalNode->changeLChild(parent, newInternalNode.get());
			newInternalNode->changeRChild(newLeafNode, newInternalNode.get());
		}
		return newLeafNode;
	}

	bool CoresetTree::insertPoint(CoresetPoint *newPoint)
	{
		this->doRefresh();
		bool usedInsertShortCut {false};

		// If there is no leafs that is smaller than required, create a new leaf and
		// insert the new point there
		if (this->currentLeaf != nullptr)
		{
			currentLeaf->insertPoint(newPoint);

			usedInsertShortCut = currentLeaf->runStaticCoreset(newPoint);

			if (currentLeaf->getNrPoints() >= parameters.s / 2)
				currentLeaf = nullptr;
			return usedInsertShortCut;
		}
		else
		{
			currentLeaf = this->insertLeaf().get();
			currentLeaf->insertPoint(newPoint);

			usedInsertShortCut = currentLeaf->runStaticCoreset(newPoint);
			return usedInsertShortCut;
		}
	}
	bool CoresetTree::removePoint(CoresetPoint *point)
	{
		this->doRefresh();

		CoresetNode *leaf = point->getRelevantLeaf();
		CoresetNode *parent = point->getRelevantLeaf()->parent;
		std::shared_ptr<CoresetNode> leaf_shared;
		if (leaf == parent->lChild.get())
			leaf_shared = parent->lChild;
		else
			leaf_shared = parent->rChild;

		for (int i = 0; i < leaf->getData().size(); i++)
		{
			if (leaf->getData().at(i)->getIndex() == point->getIndex())
			{
				leaf->erasePoint(i);
				break;
			}
		}

		if (leaf->getNrPoints() == 0)
		{
			if (leaf == currentLeaf)
				currentLeaf = nullptr;
			this->removeLeaf(leaf_shared);
			return true;
		}

		if (leaf->getNrPoints() >= parameters.s / 4 || leaf == currentLeaf)
		{
			leaf->runStaticCoreset(nullptr);
			return true;
		}

		if (this->currentLeaf == nullptr)
		{
			leaf->runStaticCoreset(nullptr);
			this->currentLeaf = leaf;
			return true;
		}

		if (leaf->parent->isRoot())
		{
			leaf->runStaticCoreset(nullptr);
			return true;
		}

		// If the function is here, the leaf where the point originally was now has
		// not enough points. All points are moved to the leaf where insertons currently
		// take place and the empty leaf is deleted
		this->movePointsBetweenLeafs(leaf, currentLeaf);
		this->removeLeaf(leaf_shared);
		currentLeaf->runStaticCoreset(nullptr);

		return true;
	}
	bool CoresetTree::removePoints(std::vector<CoresetPoint*> points)
	{
		if (this->parameters.removeCutOffFactor != -1)
			this->parameters.removeCutOff = this->parameters.nrPointsInTree * this->parameters.removeCutOffFactor;

		std::unordered_set<CoresetNode*> leafsToConsider;
		std::vector<CoresetNode*> leafsToConsiderAfter;
		std::unordered_set<CoresetNode*> internalToConsider;
		this->doRefresh();

		for (int i = 0; i < points.size(); i++)
		{
			CoresetNode *leaf = points[i]->getRelevantLeaf();
			CoresetNode *parent = points[i]->getRelevantLeaf()->parent;
			leafsToConsider.insert(leaf);

			std::shared_ptr<CoresetNode> leaf_shared;
			if (leaf == parent->lChild.get())
				leaf_shared = parent->lChild;
			else
				leaf_shared = parent->rChild;

			for (int j = 0; j < leaf->getData().size(); j++)
			{
				if (leaf->getData().at(j)->getIndex() == points[i]->getIndex())
				{
					leaf->erasePoint(j);
					break;
				}
			}
		}

		for (CoresetNode* leaf : leafsToConsider)
		{
			CoresetNode *parent = leaf->parent;
			std::shared_ptr<CoresetNode> leaf_shared;
			if (leaf == parent->lChild.get())
				leaf_shared = parent->lChild;
			else
				leaf_shared = parent->rChild;

			if (leaf->getNrPoints() == 0)
			{
				if (leaf == currentLeaf)
					currentLeaf = nullptr;
				this->removeLeaf(leaf_shared);
				continue;
			}

			if (leaf->getNrPoints() >= parameters.s / 2 || leaf == currentLeaf)
			{
				leafsToConsiderAfter.push_back(leaf);
				continue;
			}

			if (this->currentLeaf == nullptr)
			{
				leafsToConsiderAfter.push_back(leaf);
				this->currentLeaf = leaf;
				continue;
			}

			if (leaf->parent->isRoot())
			{
				leafsToConsiderAfter.push_back(leaf);
				continue;
			}

			// If the function is here, the leaf where the point originally was now has
			// not enough points. All points are moved to the leaf where insertons currently
			// take place and the empty leaf is deleted
			this->movePointsBetweenLeafs(leaf, currentLeaf);
			this->removeLeaf(leaf_shared);
			leafsToConsiderAfter.push_back(currentLeaf);
			if (currentLeaf->getNrPoints() >= parameters.s / 2)
				currentLeaf = nullptr;

		}
		if (leafsToConsiderAfter.empty())
			return 0;

		int treeHeight = this->getTreeHeight();
		for (CoresetNode* leaf : leafsToConsiderAfter)
		{
			if (leaf->getLevel() != treeHeight)
			{
				internalToConsider.insert(leaf);
				continue;
			}
			leaf->runStaticCoreset(nullptr,&internalToConsider);

		}
		if (internalToConsider.find(nullptr) != internalToConsider.end())
		{
			assert(internalToConsider.size() == 1);
			return true;
		}
		while (true)
		{
			leafsToConsider.clear();
			for (CoresetNode* leaf : internalToConsider)
				leafsToConsider.insert(leaf);
			internalToConsider.clear();

			for (CoresetNode* leaf : leafsToConsider)
			{
				leaf->runStaticCoreset(nullptr,&internalToConsider);
			}

			if (internalToConsider.find(nullptr) != internalToConsider.end())
			{
				assert(internalToConsider.size() == 1);
				break;
			}

		}
		return true;

	}

	bool CoresetTree::movePointsBetweenLeafs(CoresetNode *from, CoresetNode *to)
	{
		if (from->isInternal() || to->isInternal())
			return false;

		if (from->weightDenominator != to->weightDenominator)
		{

			std::for_each(to->weights.begin(), to->weights.end(), [&from](double &w)
						  { w *= from->weightDenominator; });

			std::for_each(from->weights.begin(), from->weights.end(), [&to](double &w)
						  { w *= to->weightDenominator; });

			to->weightDenominator *= from->weightDenominator;
		}

		to->appendWeights(from->getWeights());
		to->appendPoints(from->getData());
		from->weights.clear();
		from->data.clear();

		return true;
	}

	bool CoresetTree::swapLeafs(std::shared_ptr<CoresetNode> leaf1, std::shared_ptr<CoresetNode> leaf2)
	{
		if (leaf1->isInternal() || leaf2->isInternal())
			return false;

		CoresetNode *parent1 = leaf1->parent;
		CoresetNode *parent2 = leaf2->parent;

		if (parent1->lChild.get() == leaf1.get())
		{
			parent1->lChild.reset();
			parent1->lChild = leaf2;
		}
		else
		{
			parent1->rChild.reset();
			parent1->rChild = leaf2;
		}

		if (parent2->lChild.get() == leaf2.get())
		{
			parent2->lChild.reset();
			parent2->lChild = leaf1;
		}
		else
		{
			parent2->rChild.reset();
			parent2->rChild = leaf1;
		}

		leaf1->parent = parent2;

		leaf2->parent = parent1;

		return true;
	}

	void CoresetTree::writeToStream(std::ostream &out)
	{
		std::queue<std::shared_ptr<CoresetNode>> queue;

		queue.push(this->root);

		int counter{1};
		int counter2{1};
		while (!queue.empty())
		{
			auto i = queue.front();
			if (counter2 == 1)
				out << "Level: ";
			// out << "Node: " << i->getLChild().get() << "-" << i->getRChild().get() << " | ";
			out << "Node " << i.get() << ", " << (this->refreshPointer == i.get() ? "r" : "") << ": " << i->getNrPoints();
			// i->writeIndicesToStream(out);
			if (i->isLeaf())
				out << ", l | ";
			else
				out << ", i | ";

			queue.pop();
			if (!i->isLeaf())
			{
				queue.push(i->getLChild());
				queue.push(i->getRChild());
			}
			if (counter2 == counter)
			{
				out << "\n";
				counter2 = 0;
				counter *= 2;
			}
			counter2++;
		}
		out << "\n";
	}

	void CoresetTree::doRefresh(bool onlyMove)
	{
		// If the refresh pointer is a nullpointer, all leafs have already been
		// refreshed this phase. Then return
		if (refreshPointer == nullptr)
			return;

		if (!onlyMove)
			(*this->refreshPointer).runStaticCoreset(nullptr);

		refreshPointer = this->getLeafToRight(refreshPointer);
	}

	void CoresetTree::checkIfLeafsToSmall()
	{
		while (refreshPointer != nullptr)
		{
			if (this->refreshPointer->getNrPoints() < parameters.s / 2 && this->refreshPointer != this->currentLeaf)
			{
				if (this->toSmallLeaf == nullptr)
				{
					toSmallLeaf = refreshPointer;
					refreshPointer = this->getLeafToRight(refreshPointer);
				}
				// There are now to leafs that are to small (the one pointed to by the refresh pointer and the one pointed to by the "toSmallLeaf" pointer).
				// Those two leafs are now combined by swapping those two leafs with the positions of the rightmost leaf and its sibling,
				// copying the points from one of them to the other and then removing the empty leaf.
				else
				{
					if (this->currentLeaf == refreshPointer)
						currentLeaf = nullptr;

					std::shared_ptr<CoresetNode> RightmostLeaf1 = this->getRightmostLeaf();
					std::shared_ptr<CoresetNode> RightmostLeaf2 = RightmostLeaf1->getSibling();

					// Get the shared pointers to the two small leafs
					std::shared_ptr<CoresetNode> refreshLeaf_shared;
					std::shared_ptr<CoresetNode> toSmallLeaf_shared;
					if (refreshPointer == refreshPointer->getParent()->lChild.get())
						refreshLeaf_shared = refreshPointer->getParent()->lChild;
					else
						refreshLeaf_shared = refreshPointer->getParent()->rChild;

					if (toSmallLeaf == toSmallLeaf->getParent()->lChild.get())
						toSmallLeaf_shared = toSmallLeaf->getParent()->lChild;
					else
						toSmallLeaf_shared = toSmallLeaf->getParent()->rChild;

					CoresetNode *rightParent = RightmostLeaf1->parent;

					// Swap the two leafs that are to small with the two rightmost leafs on the lowest level
					if (toSmallLeaf != RightmostLeaf1.get() && toSmallLeaf != RightmostLeaf2.get() && refreshPointer != RightmostLeaf1.get() && refreshPointer != RightmostLeaf2.get())
					{
						this->swapLeafs(toSmallLeaf_shared, RightmostLeaf2);
						this->swapLeafs(refreshLeaf_shared, RightmostLeaf1);

						refreshPointer = RightmostLeaf2.get();
					}
					else if (toSmallLeaf != RightmostLeaf2.get() && toSmallLeaf != RightmostLeaf1.get() && refreshPointer == RightmostLeaf2.get())
					{
						this->swapLeafs(toSmallLeaf_shared, RightmostLeaf1);

						refreshPointer = RightmostLeaf1.get();
					}
					else if (toSmallLeaf != RightmostLeaf2.get() && toSmallLeaf != RightmostLeaf1.get() && refreshPointer == RightmostLeaf1.get())
					{
						this->swapLeafs(toSmallLeaf_shared, RightmostLeaf2);

						refreshPointer = RightmostLeaf2.get();
					}
					else if (toSmallLeaf == RightmostLeaf2.get() && refreshPointer != RightmostLeaf1.get() && refreshPointer != RightmostLeaf2.get())
					{
						this->swapLeafs(refreshLeaf_shared, RightmostLeaf1);

						refreshPointer = toSmallLeaf;
					}
					else if (toSmallLeaf == RightmostLeaf1.get() && refreshPointer != RightmostLeaf1.get() && refreshPointer != RightmostLeaf2.get())
					{
						this->swapLeafs(refreshLeaf_shared, RightmostLeaf2);

						refreshPointer = toSmallLeaf;
					}
					else
						refreshPointer = toSmallLeaf;

					for (auto &d : refreshLeaf_shared->data)
						d->setRelevantLeaf(toSmallLeaf);

					// The centers need to be empty because leafs do not calculate coresets
					// Then move all Points from the refreshLeaf to the toSmallLeaf
					assert(toSmallLeaf_shared->centers.size() == 0);
					assert(refreshLeaf_shared->centers.size() == 0);

					toSmallLeaf_shared->data.insert(toSmallLeaf_shared->data.end(), refreshLeaf_shared->data.begin(), refreshLeaf_shared->data.end());
					refreshLeaf_shared->data.clear();

					// Change the weights if the weight denominators were different then move the weights between the leafs
					if (toSmallLeaf_shared->weightDenominator != refreshLeaf_shared->weightDenominator)
					{

						std::for_each(toSmallLeaf_shared->weights.begin(), toSmallLeaf_shared->weights.end(), [&refreshLeaf_shared](double &w)
									  { w *= refreshLeaf_shared->weightDenominator; });

						std::for_each(refreshLeaf_shared->weights.begin(), refreshLeaf_shared->weights.end(), [&toSmallLeaf_shared](double &w)
									  { w *= toSmallLeaf_shared->weightDenominator; });

						toSmallLeaf_shared->weightDenominator *= refreshLeaf_shared->weightDenominator;
					}
					toSmallLeaf_shared->weights.insert(toSmallLeaf_shared->weights.end(), refreshLeaf_shared->weights.begin(), refreshLeaf_shared->weights.end());
					refreshLeaf_shared->weights.clear();

					// remove the now empyt refresh leaf, swap the toSmallLeaf with its parent and delete the old parent
					CoresetNode *rightGrandParent = rightParent->parent;

					if (rightGrandParent->lChild.get() == rightParent)
					{
						rightGrandParent->lChild.reset();
						rightGrandParent->lChild = toSmallLeaf_shared;
					}
					else
					{
						rightGrandParent->rChild.reset();
						rightGrandParent->rChild = toSmallLeaf_shared;
					}

					toSmallLeaf_shared->parent = rightGrandParent;

					toSmallLeaf_shared->runStaticCoreset(nullptr);

					toSmallLeaf = nullptr;
				}
			}
			else
				refreshPointer = this->getLeafToRight(refreshPointer);
		}
	}

	void CoresetTree::checkIfLeafsToLarge()
	{
		while (refreshPointer != nullptr)
		{
			// If the number of Points in a leaf are to large, swap the leaf with the leftmost leaf on the second lowest level
			// or with the leftmost leaf on the lowest level if there is only one level that contains leafs.
			// Then create a new internal node that gets inserted between the leaf that is to large ans its parent
			// and create a sibling for the leaf. Then move halve of the poits from the to large leaf to its sibling.
			if (this->refreshPointer->getNrPoints() > parameters.s)
			{
				if (this->currentLeaf == refreshPointer)
					currentLeaf = nullptr;

				std::shared_ptr<CoresetNode> LeftmostLeaf = this->getLeftmostLeaf();

				std::shared_ptr<CoresetNode> refreshLeaf_shared;
				if (refreshPointer == refreshPointer->getParent()->lChild.get())
					refreshLeaf_shared = refreshPointer->getParent()->lChild;
				else
					refreshLeaf_shared = refreshPointer->getParent()->rChild;

				if (refreshLeaf_shared != LeftmostLeaf)
				{
					this->swapLeafs(refreshLeaf_shared, LeftmostLeaf);

					refreshPointer = LeftmostLeaf.get();
				}
				else
					refreshPointer = this->getLeafToRight(refreshPointer);

				std::shared_ptr<CoresetNode> newInternal = std::make_shared<CoresetNode>(parameters);
				std::shared_ptr<CoresetNode> newRChild = std::make_shared<CoresetNode>(parameters);

				assert(refreshLeaf_shared->centers.size() == 0);

				// Move the points and weights to the new leaf
				int nrPointsInRChild = refreshLeaf_shared->data.size() / 2;
				newRChild->weightDenominator = refreshLeaf_shared->weightDenominator;

				newRChild->data.insert(newRChild->data.end(), refreshLeaf_shared->data.begin() + nrPointsInRChild, refreshLeaf_shared->data.end());
				refreshLeaf_shared->data.erase(refreshLeaf_shared->data.begin() + nrPointsInRChild, refreshLeaf_shared->data.end());

				// Set the correst relevant leaf for the points that were moved between leafs
				for (auto &d : newRChild->data)
					d->setRelevantLeaf(newRChild.get());

				newRChild->weights.insert(newRChild->weights.end(), refreshLeaf_shared->weights.begin() + nrPointsInRChild, refreshLeaf_shared->weights.end());
				refreshLeaf_shared->weights.erase(refreshLeaf_shared->weights.begin() + nrPointsInRChild, refreshLeaf_shared->weights.end());

				CoresetNode *refreshParent = refreshLeaf_shared->parent;

				if (refreshParent->lChild == refreshLeaf_shared)
				{
					refreshParent->lChild.reset();
					refreshParent->lChild = newInternal;
				}
				else
				{
					refreshParent->rChild.reset();
					refreshParent->rChild = newInternal;
				}

				newInternal->parent = refreshParent;
				newInternal->lChild = refreshLeaf_shared;
				newInternal->rChild = newRChild;

				refreshLeaf_shared->parent = newInternal.get();
				newRChild->parent = newInternal.get();

				newRChild->runStaticCoreset(nullptr);
				refreshLeaf_shared->runStaticCoreset(nullptr);
			}
			else
				refreshPointer = this->getLeafToRight(refreshPointer);
		}
	}

	void CoresetTree::resetRefresh(bool nrPointsIncreased)
	{

		// Find the leftmost leaf on the lowest level
		CoresetNode *parent = this->root.get();
		CoresetNode *child = parent->lChild.get();
		while (child->lChild != nullptr)
		{
			parent = child;
			child = parent->lChild.get();
		}
		refreshPointer = child;
		this->toSmallLeaf = nullptr;
		if (nrPointsIncreased)
			this->checkIfLeafsToSmall();
		else
			this->checkIfLeafsToLarge();
		
		refreshPointer = child;
	}

	CoresetNode *CoresetTree::getLeafToRight(CoresetNode *leaf)
	{
		CoresetNode *ret;
		CoresetNode *refreshParent = leaf->getParent();
		CoresetNode *refreshChild = leaf;

		// Moves up in the tree, until the last node is a left child of its parent.
		while (refreshParent->lChild.get() != refreshChild)
		{
			if (refreshParent->parent == nullptr)
			{
				return nullptr;
			}
			refreshChild = refreshParent;
			refreshParent = refreshChild->getParent();
		}

		// Then move to the right child of that parent and continue moving down alsways to the left
		// child until a leaf is hit
		refreshChild = refreshParent->rChild.get();
		while (refreshChild->lChild.get() != nullptr)
		{
			refreshParent = refreshChild;
			refreshChild = refreshParent->lChild.get();
		}

		return refreshChild;
	}
}
