#include "coreset_tree.h"
#include "blaze/math/expressions/DVecDVecMapExpr.h"
#include <cmath>
#include <memory>
#include <queue>

namespace dynClustering
{

	std::shared_ptr<CoresetNode> CoresetTree::insertLeaf()
	{
		std::shared_ptr<CoresetNode> newLeafNode = std::make_shared<CoresetNode>(numberOfClusters*clusterFactor,s, numberOfNotCenterCoresetPoints, cutOff, np, eps);

		std::shared_ptr<CoresetNode> newInternalNode1;
		std::shared_ptr<CoresetNode> newInternalNode2;
		std::shared_ptr<CoresetNode> newInternalNode3;
		std::shared_ptr<CoresetNode> newInternalNode4;
		std::shared_ptr<CoresetNode> newInternalNode5;

		if (leafs.size() % (this->cardOfInternalNodes*this->cardOfInternalNodes*this->cardOfInternalNodes*this->cardOfInternalNodes*this->cardOfInternalNodes) == 0)
		{
			newInternalNode5 = std::make_shared<CoresetNode>(numberOfClusters*clusterFactor,s, numberOfNotCenterCoresetPoints, cutOff, np, eps);
			intermediate_level5.push_back(newInternalNode5);
		}
		else
			newInternalNode5 = this->intermediate_level5.at(this->intermediate_level5.size()-1);

		if (leafs.size() % (this->cardOfInternalNodes*this->cardOfInternalNodes*this->cardOfInternalNodes*this->cardOfInternalNodes) == 0)
		{
			newInternalNode4 = std::make_shared<CoresetNode>(numberOfClusters*clusterFactor,s, numberOfNotCenterCoresetPoints, cutOff, np, eps);
			intermediate_level4.push_back(newInternalNode4);
			newInternalNode4->setRoot(newInternalNode5.get());
			newInternalNode5->children.push_back(newInternalNode4.get());
		}
		else
			newInternalNode4 = this->intermediate_level4.at(this->intermediate_level4.size()-1);

		if (leafs.size() % (this->cardOfInternalNodes*this->cardOfInternalNodes*this->cardOfInternalNodes) == 0)
		{
			newInternalNode3 = std::make_shared<CoresetNode>(numberOfClusters*clusterFactor,s, numberOfNotCenterCoresetPoints, cutOff, np, eps);
			intermediate_level3.push_back(newInternalNode3);
			newInternalNode3->setRoot(newInternalNode4.get());
			newInternalNode4->children.push_back(newInternalNode3.get());
		}
		else
			newInternalNode3 = this->intermediate_level3.at(this->intermediate_level3.size()-1);

		if (leafs.size() % (this->cardOfInternalNodes*this->cardOfInternalNodes) == 0)
		{
			newInternalNode2 = std::make_shared<CoresetNode>(numberOfClusters*clusterFactor,s, numberOfNotCenterCoresetPoints, cutOff, np, eps);
			intermediate_level2.push_back(newInternalNode2);
			newInternalNode2->setRoot(newInternalNode3.get());
			newInternalNode3->children.push_back(newInternalNode2.get());
		}
		else
			newInternalNode2 = this->intermediate_level2.at(this->intermediate_level2.size()-1);

		if (leafs.size() % this->cardOfInternalNodes == 0)
		{
			newInternalNode1 = std::make_shared<CoresetNode>(numberOfClusters*clusterFactor,s, numberOfNotCenterCoresetPoints, cutOff, np, eps);
			intermediate_level1.push_back(newInternalNode1);
			newInternalNode1->setRoot(newInternalNode2.get());
			newInternalNode2->children.push_back(newInternalNode1.get());
		}
		else
			newInternalNode1 = this->intermediate_level1.at(this->intermediate_level1.size()-1);


			
		newLeafNode->setRoot(newInternalNode1.get());
		newInternalNode1->children.push_back(newLeafNode.get());

		this->leafs.push_back(newLeafNode);
		return newLeafNode;
	}


	void CoresetTree::removeLeaf(std::shared_ptr<CoresetNode> LeafToRemove)
	{
		assert(LeafToRemove->data.size() == 0);
		for (int i = 0; i < leafs.size(); i++)
		{
			if (LeafToRemove == leafs.at(i))
			{
				leafs.erase(leafs.begin()+i);
				break;
			}
		}
		return;


	}


	bool CoresetTree::insertPoint(CoresetPoint *newPoint, bool isInit)
	{


		if (currentLeaf->getNrPoints() < this->maxPointsInLeafs)
		{
			currentLeaf->insertPoint(newPoint);
			this->runStaticCoreset(currentLeaf, isInit);

			this->n++;

			return 0;

		}
		else
		{

			bool needTreeUpdate {true};
			for(auto leaf = leafs.rbegin(); leaf != leafs.rend(); leaf++)
			{
				if ((*leaf)->getNrPoints() < maxPointsInLeafs)
				{
					if ((*leaf)->getNrPoints() == maxPointsInLeafs - 1)
					{
						(*leaf)->insertPoint(newPoint);

						this->runStaticCoreset((*leaf).get(), isInit);
						this->n++;
						return 0;
					}
					currentLeaf = (*leaf).get();
					needTreeUpdate = false;
					break;
				}
			}
			if (needTreeUpdate)
			{
				currentLeaf = this->insertLeaf().get();
				this->maxPointsInLeafs = blaze::max(1,std::floor(this->leafs.size()*this->maxPointsFactor));
			}

			currentLeaf->insertPoint(newPoint);

			this->runStaticCoreset(currentLeaf, isInit);
			this->n++;
			return 0;
		}

		
	}
	bool CoresetTree::removePoint(CoresetPoint* point)
	{

		auto leaf = point->getRelevantLeaf();

		for (int i = 0; i < leaf->getData().size(); i++)
		{
			if (leaf->getData().at(i)->getIndex() == point->getIndex())
			{
				leaf->erasePoint(i);
				break;
			}
		}
		if (this->n < (this->leafs.size()-1)*blaze::max(1, std::floor((this->leafs.size()-1) * this->maxPointsFactor)))
		{
			this->maxPointsInLeafs = blaze::max(1, std::floor((this->leafs.size()-1) * this->maxPointsFactor));
			this->reduceTreeSize();
			this->currentLeaf = this->leafs.at(this->leafs.size()-1).get();
			for (const std::shared_ptr<CoresetNode> l : this->leafs)
				this->runStaticCoreset(l.get());

			this->n--;
			return true;
		}

		this->runStaticCoreset(leaf);
		if (this->currentLeaf == nullptr)
			currentLeaf = leaf;

		this->n--;
		return true;

	}

	void CoresetTree::writeToStream(std::ostream &out)
	{
		out << "n: " << this->n << "\n";
		out << "Tree:\nRoot: " << this->root->data.size() << "\nLeafs: ";
		for (const auto& leaf : leafs)
			out << leaf->getNrPoints() << " | ";
		out << std::endl;

	}

	void CoresetTree::reduceTreeSize()
	{


		int openLeaf {0};
		for (auto& point : leafs.at(leafs.size()-1)->data)
		{
			for (; openLeaf < leafs.size()-1; openLeaf++)
			{
				if (leafs.at(openLeaf)->getNrPoints() < this->maxPointsInLeafs)
				{
					leafs.at(openLeaf)->insertPoint(point);
					break;
				}
			}

		}
		leafs.at(leafs.size()-1)->clearData();

		if (leafs.at(leafs.size()-1)->root->children.size() == 1)
		{
			leafs.at(leafs.size()-1)->root->children.clear();
			leafs.at(leafs.size()-1)->root->data.clear();
			leafs.at(leafs.size()-1)->root->weights.clear();
			leafs.at(leafs.size()-1)->root->resultData.clear();
			leafs.at(leafs.size()-1)->root->resultWeights.clear();
		}
		else
		{
			auto child = std::find(leafs.at(leafs.size()-1)->root->children.begin(), leafs.at(leafs.size()-1)->root->children.end(),leafs.at(leafs.size()-1).get());
			leafs.at(leafs.size()-1)->root->children.erase(child);
			this->runStaticCoreset(leafs.at(leafs.size()-1)->root->children.at(0));
		}

		leafs.erase(leafs.begin() + leafs.size()-1);



		
		for (int i = 0; i < leafs.size(); i++)
		{
			if (leafs.at(i)->getNrPoints() > this->maxPointsInLeafs)
			{
				for (int j = this->maxPointsInLeafs; j < leafs.at(i)->getNrPoints(); j++)
				{
					for (; openLeaf < leafs.size(); openLeaf++)
					{
						if (leafs.at(openLeaf)->getNrPoints() < this->maxPointsInLeafs)
						{
							leafs.at(openLeaf)->insertPoint(leafs.at(i)->data.at(j));
							break;
						}
					}
				}
				while (leafs.at(i)->getNrPoints() != this->maxPointsInLeafs)
					leafs.at(i)->erasePoint(leafs.at(i)->getNrPoints()-1);
				assert(leafs.at(i)->getNrPoints() == this->maxPointsInLeafs);
				assert(leafs.at(i)->data.size() == leafs.at(i)->weights.size());
			}
		}

	for (const auto &leaf : leafs)
			assert(leaf->getNrPoints() <= maxPointsInLeafs);

	}

	void CoresetTree::runStaticCoreset(CoresetNode* startLeaf, bool isInit)
	{
		if (isInit || this->numberOfClusters >= startLeaf->data.size())
		{
			startLeaf->resultData = startLeaf->data;
			startLeaf->resultWeights = startLeaf->weights;

			startLeaf->root->clearVectors();
			startLeaf->root->root->clearVectors();
			startLeaf->root->root->root->clearVectors();
			startLeaf->root->root->root->root->clearVectors();
			startLeaf->root->root->root->root->root->clearVectors();

			for (auto &child : startLeaf->root->children)
			{
				startLeaf->root->data.insert(startLeaf->root->data.end(),child->resultData.begin(),child->resultData.end());
				startLeaf->root->weights.insert(startLeaf->root->weights.end(),child->resultWeights.begin(),child->resultWeights.end());
			}
			startLeaf->root->runStaticCoreset();

			for (auto &child : startLeaf->root->root->children)
			{
				startLeaf->root->root->data.insert(startLeaf->root->root->data.end(),child->resultData.begin(),child->resultData.end());
				startLeaf->root->root->weights.insert(startLeaf->root->root->weights.end(),child->resultWeights.begin(),child->resultWeights.end());
			}
			startLeaf->root->root->runStaticCoreset();

			for (auto &child : startLeaf->root->root->root->children)
			{
				startLeaf->root->root->root->data.insert(startLeaf->root->root->root->data.end(),child->resultData.begin(),child->resultData.end());
				startLeaf->root->root->root->weights.insert(startLeaf->root->root->root->weights.end(),child->resultWeights.begin(),child->resultWeights.end());
			}
			startLeaf->root->root->root->runStaticCoreset();

			for (auto &child : startLeaf->root->root->root->root->children)
			{
				startLeaf->root->root->root->root->data.insert(startLeaf->root->root->root->root->data.end(),child->resultData.begin(),child->resultData.end());
				startLeaf->root->root->root->root->weights.insert(startLeaf->root->root->root->root->weights.end(),child->resultWeights.begin(),child->resultWeights.end());
			}
			startLeaf->root->root->root->root->runStaticCoreset();

			for (auto &child : startLeaf->root->root->root->root->root->children)
			{
				startLeaf->root->root->root->root->root->data.insert(startLeaf->root->root->root->root->root->data.end(),child->resultData.begin(),child->resultData.end());
				startLeaf->root->root->root->root->root->weights.insert(startLeaf->root->root->root->root->root->weights.end(),child->resultWeights.begin(),child->resultWeights.end());
			}
			startLeaf->root->root->root->root->root->runStaticCoreset();

			root->clearData();
			//!!!Change here also
			for (const auto& internalNode : this->intermediate_level5)
			{
				this->root->data.insert(this->root->data.end(),internalNode->resultData.begin(),internalNode->resultData.end());
				this->root->weights.insert(this->root->weights.end(),internalNode->resultWeights.begin(),internalNode->resultWeights.end());
			}
			if (isInit || root->data.size() <= this->sizeOfFinalCoreset)
			{
				this->root->resultData = this->root->data;
				this->root->resultWeights = this->root->weights;
				return;
			}
			root->runStaticCoreset();
			return;
			
		}
		startLeaf->runStaticCoreset();

		startLeaf->root->clearVectors();
		startLeaf->root->root->clearVectors();
		startLeaf->root->root->root->clearVectors();
		startLeaf->root->root->root->root->clearVectors();
		startLeaf->root->root->root->root->root->clearVectors();

		for (auto &child : startLeaf->root->children)
		{
			startLeaf->root->data.insert(startLeaf->root->data.end(),child->resultData.begin(),child->resultData.end());
			startLeaf->root->weights.insert(startLeaf->root->weights.end(),child->resultWeights.begin(),child->resultWeights.end());
		}
		startLeaf->root->runStaticCoreset();

		for (auto &child : startLeaf->root->root->children)
		{
			startLeaf->root->root->data.insert(startLeaf->root->root->data.end(),child->resultData.begin(),child->resultData.end());
			startLeaf->root->root->weights.insert(startLeaf->root->root->weights.end(),child->resultWeights.begin(),child->resultWeights.end());
		}
		startLeaf->root->root->runStaticCoreset();

		for (auto &child : startLeaf->root->root->root->children)
		{
			startLeaf->root->root->root->data.insert(startLeaf->root->root->root->data.end(),child->resultData.begin(),child->resultData.end());
			startLeaf->root->root->root->weights.insert(startLeaf->root->root->root->weights.end(),child->resultWeights.begin(),child->resultWeights.end());
		}
		startLeaf->root->root->root->runStaticCoreset();

		for (auto &child : startLeaf->root->root->root->root->children)
		{
			startLeaf->root->root->root->root->data.insert(startLeaf->root->root->root->root->data.end(),child->resultData.begin(),child->resultData.end());
			startLeaf->root->root->root->root->weights.insert(startLeaf->root->root->root->root->weights.end(),child->resultWeights.begin(),child->resultWeights.end());
		}
		startLeaf->root->root->root->root->runStaticCoreset();

		for (auto &child : startLeaf->root->root->root->root->root->children)
		{
			startLeaf->root->root->root->root->root->data.insert(startLeaf->root->root->root->root->root->data.end(),child->resultData.begin(),child->resultData.end());
			startLeaf->root->root->root->root->root->weights.insert(startLeaf->root->root->root->root->root->weights.end(),child->resultWeights.begin(),child->resultWeights.end());
		}
		startLeaf->root->root->root->root->root->runStaticCoreset();
		
		this->root->clearData();

		//!!Also change here
		for (const auto& internalNode : this->intermediate_level5)
		{
			this->root->data.insert(this->root->data.end(),internalNode->resultData.begin(),internalNode->resultData.end());
			this->root->weights.insert(this->root->weights.end(),internalNode->resultWeights.begin(),internalNode->resultWeights.end());
		}

		root->runStaticCoreset();
	}


}
