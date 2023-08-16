#pragma once

#include <ostream>
#include <vector>
#include <memory>
#include <boost/unordered/unordered_map.hpp>
#include "coreset_point.h"
#include "coreset_node.h"
#include "coreset_tree.h"

namespace dynClustering
{
	class DynCoreset
	{
		public:
		std::shared_ptr<CoresetNode> currentNode;
		boost::unordered_map<unsigned int,CoresetPoint> data;
		CoresetTree* tree;
		coresets::SensitivitySampling staticCoreset;
		int s;
		int numberOfNotCenterCoresetPoints;
		int numberOfNotCenterCoresetPointsNode;
		size_t numberOfClusters;
		int sizeOfFinalCoreset;
		int dimensions;
		std::function<int(double epsilon, double lambda, int k, int d)> cardinalityOfCoreset;
		double epsilon;
		double lambda;
		int np;
		int updatesThisPhase;
		int clusterFactor;
		double sFactor;
		double cutOffFactor;
		int cutOff;
		double maxPointsFactor;

		std::vector<CoresetPoint*> data2Proxy {};
		std::vector<double> weights2Proxy {};
		std::vector<CoresetPoint*> centers {};
		std::vector<CoresetPoint*> coresetPoints {};
		std::vector<double> coresetWeights {};
		void WriteCoresetToStream(std::ostream &out);
		
		

		
		DynCoreset(size_t numberOfClusters,std::function<int(double epsilon, double lambda, int k, int d)> cardinalityOfCoreset, double epsilon, double lambda, int dimensions, int clusterFactor, double sFactor, double cutOffFactor, double maxPointsFactor) : staticCoreset({numberOfClusters,&numberOfNotCenterCoresetPoints})
		{
			this->updatesThisPhase = 0;
			this->dimensions = dimensions;
			this->epsilon = epsilon;
			this->lambda = lambda;
			this->numberOfClusters = numberOfClusters;
			this->sizeOfFinalCoreset = cardinalityOfCoreset(epsilon/3,lambda/2,numberOfClusters,dimensions);
			this->cardinalityOfCoreset = cardinalityOfCoreset;
			this->clusterFactor = clusterFactor;
			this->sFactor = sFactor;
			this->cutOffFactor = cutOffFactor;
			this->maxPointsFactor = maxPointsFactor;

		}
		~DynCoreset()
		{
			delete tree;

			for(int i = centers.size()-1; i >= 0; i--)
				delete centers.at(i);
		}

		void writeTreeToStream(std::ostream& out)
		{
			out << "Final Coreset: " << this->coresetPoints.size() << "\n";
			tree->writeToStream(out);
		}


		bool initTree(std::shared_ptr<blaze::DynamicMatrix<double>> firstSPoints);

		bool insertPoint(int index, double weight, std::vector<double>& coordinates);
		bool insertPoint(blaze::DynamicMatrix<double>& point, bool isInit = false);
		bool removePoint(int index);
		bool applyOperation(blaze::DynamicMatrix<double>& point);
		void goToNextPhase(bool nrPointsIncreased);
		void runOuterCoreset();

		

	};
}
