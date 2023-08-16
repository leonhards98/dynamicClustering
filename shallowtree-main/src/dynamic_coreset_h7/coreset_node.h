#pragma once
#include <algorithm>
#include <vector>
#include <memory>
#include "coreset_point.h"
#include "../coreset/sensitivity_sampling.hpp"

namespace dynClustering
{
	class CoresetTree;
	class CoresetNode
	{
		std::vector<CoresetPoint*> data;
		std::vector<CoresetPoint*> resultData;
		std::vector<double> weights;
		std::vector<double> resultWeights;
		std::vector<CoresetNode*> children;

		coresets::SensitivitySampling staticCoreset;

		std::vector<CoresetPoint*> centers;
		CoresetNode* root;

		int *sizeOfCoreset;
		int* cutOff;
		int *np;
		int* numberOfNotCenterCoresetPoints;

		double eps;
		friend CoresetTree;


		public:

		CoresetNode(size_t numberOfclusters, int *sizeOfCoreset, int* numberOfNotCenterCoresetPoints, int *cutOff, int* np, double eps) : sizeOfCoreset(sizeOfCoreset), numberOfNotCenterCoresetPoints(numberOfNotCenterCoresetPoints), cutOff(cutOff), np(np), eps(eps), staticCoreset({numberOfclusters,numberOfNotCenterCoresetPoints}) 
		{
		}

		~CoresetNode()
		{
			for(int i = centers.size()-1; i >= 0; i--)
				delete centers.at(i);
		}

		bool insertPoint(CoresetPoint* point);
		void runStaticCoreset();

		void clearVectors()
		{
			this->resultData.clear();
			this->data.clear();
			this->resultWeights.clear();
			this->weights.clear();
		}


		int getNrPoints()
		{
			return this->data.size();
		}
		void setSizeOfCoreset(int *s)
		{
			this->sizeOfCoreset = s;
		}
		const std::vector<CoresetPoint*>& getData() {return data;}
		const std::vector<CoresetPoint*>& getResultData() {return resultData;}
		const std::vector<double>& getWeights() {return weights;}
		const std::vector<double>& getResultWeights() {return resultWeights;}
		const std::vector<CoresetPoint*>& getCenters() {return centers;}
		void setRoot(CoresetNode* root) {this->root = root;}
		void writePointsToStream(std::ostream &out);
		void writeIndicesToStream(std::ostream &out);

		void changeLChild(std::shared_ptr<CoresetNode> newChild, CoresetNode* newParent);
		void changeRChild(std::shared_ptr<CoresetNode> newChild, CoresetNode* newParent);
		void erasePoint(int i);
		void appendWeights(const std::vector<double> &newWeights);
		void appendPoints(const std::vector<CoresetPoint*> &newPoints);
		void clearData();


	};
}
