#pragma once

#include "coreset_point.h"
#include "coreset_node.h"
#include <memory>
#include <ostream>
#include <functional>
#include <vector>

namespace dynClustering
{
	class CoresetTree
	{
		std::vector<std::shared_ptr<CoresetNode>> leafs;
		std::vector<std::shared_ptr<CoresetNode>> intermediate_level1;
		std::vector<std::shared_ptr<CoresetNode>> intermediate_level2;
		std::vector<std::shared_ptr<CoresetNode>> intermediate_level3;
		std::vector<std::shared_ptr<CoresetNode>> intermediate_level4;
		std::vector<std::shared_ptr<CoresetNode>> intermediate_level5;
		int cardOfInternalNodes = 3;
		std::shared_ptr<CoresetNode> root;
		CoresetNode* currentLeaf;
		int* s;
		int* cutOff;
		int* np;
		int* numberOfNotCenterCoresetPoints;
		double eps;
		int numberOfClusters;
		int sizeOfFinalCoreset;
		int dimensions;
		int clusterFactor;
		int maxPointsInLeafs;
		double maxPointsFactor;
		int n;

		std::shared_ptr<CoresetNode> insertLeaf();
		void removeLeaf(std::shared_ptr<CoresetNode> LeafToRemove);
		std::shared_ptr<CoresetNode> insertChild(std::shared_ptr<CoresetNode> &parent);


	public:
		CoresetTree(int numberOfClusters, int dimensions, int sizeOfFinalCoreset, int *s, int *numberOfNotCenterCoresetPoints, int *cutOff, int* np, double eps, int clusterFactor, double maxPointsFactor)
		{
			this->dimensions = dimensions;
			this->numberOfClusters = numberOfClusters;
			this->sizeOfFinalCoreset = sizeOfFinalCoreset;
			this->s = s;
			this->numberOfNotCenterCoresetPoints = numberOfNotCenterCoresetPoints;
			this->np = np;
			this->n = 0;
			this->eps = eps;
			this->clusterFactor = clusterFactor;
			this->cutOff = cutOff;
			this->maxPointsFactor = maxPointsFactor;
			auto x = leafs.begin();

			root = std::make_shared<CoresetNode>(numberOfClusters*clusterFactor, s, numberOfNotCenterCoresetPoints, cutOff, np, eps);
			root->setRoot(nullptr);
			std::shared_ptr<CoresetNode> Child = std::make_shared<CoresetNode>(numberOfClusters*clusterFactor, s, numberOfNotCenterCoresetPoints, cutOff, np, eps);
			std::shared_ptr<CoresetNode> internalNode1 = std::make_shared<CoresetNode>(numberOfClusters*clusterFactor, s, numberOfNotCenterCoresetPoints, cutOff, np, eps);
			std::shared_ptr<CoresetNode> internalNode2 = std::make_shared<CoresetNode>(numberOfClusters*clusterFactor, s, numberOfNotCenterCoresetPoints, cutOff, np, eps);
			std::shared_ptr<CoresetNode> internalNode3 = std::make_shared<CoresetNode>(numberOfClusters*clusterFactor, s, numberOfNotCenterCoresetPoints, cutOff, np, eps);
			std::shared_ptr<CoresetNode> internalNode4 = std::make_shared<CoresetNode>(numberOfClusters*clusterFactor, s, numberOfNotCenterCoresetPoints, cutOff, np, eps);
			std::shared_ptr<CoresetNode> internalNode5 = std::make_shared<CoresetNode>(numberOfClusters*clusterFactor, s, numberOfNotCenterCoresetPoints, cutOff, np, eps);

			Child->setRoot(internalNode1.get());
			leafs.push_back(Child);

			internalNode1->children.push_back(Child.get());
			internalNode1->setRoot(internalNode2.get());
			intermediate_level1.push_back(internalNode1);

			internalNode2->children.push_back(internalNode1.get());
			internalNode2->setRoot(internalNode3.get());
			intermediate_level2.push_back(internalNode2);

			internalNode3->children.push_back(internalNode2.get());
			internalNode3->setRoot(internalNode4.get());
			intermediate_level3.push_back(internalNode3);

			internalNode4->children.push_back(internalNode3.get());
			internalNode4->setRoot(internalNode5.get());
			intermediate_level4.push_back(internalNode4);

			internalNode5->children.push_back(internalNode4.get());
			intermediate_level5.push_back(internalNode5);

			this->currentLeaf = Child.get();
			this->maxPointsInLeafs = 1;

		}
		void setS(int *s)
		{
			this->s = s;
		}

		int getNrLeafs()
		{
			return this->leafs.size();
		}

		//bool insertPoint(unsigned int index, double weight, std::vector<double> &coordinates);
		bool insertPoint(CoresetPoint* newPoint, bool isInit = false);
		bool removePoint(CoresetPoint* point);
		void reduceTreeSize();
		void writeToStream(std::ostream &out);
		//void writeCoresetToStream(std::ostream &out);
		void runStaticCoreset(CoresetNode* startLeaf, bool isInit = false);
		const std::shared_ptr<CoresetNode>& getRoot() {return root;}
	};

}
