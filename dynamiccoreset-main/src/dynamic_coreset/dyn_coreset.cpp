#include "dyn_coreset.h"
#include "boost/unordered/unordered_map_fwd.hpp"
#include "coreset_point.h"
#include <algorithm>
#include <memory>
#include <vector>
#include <iostream>

namespace dynClustering
{
	bool DynCoreset::initTree(std::shared_ptr<blaze::DynamicMatrix<double>> firstSPoints)
	{
		
		this->staticCoreset.random.initialize(time(nullptr));
		int n0 = firstSPoints->rows();
		this->parameters.np = 4*n0;

		double s_outer = this->cardinalityOfCoreset(parameters.epsilon/(6*blaze::ceil(blaze::log(parameters.np))),parameters.lambda/(2*parameters.np), parameters.numberOfClusters,parameters.dimensions);
		this->parameters.s = this->cardinalityOfCoreset(parameters.epsilon/(6*blaze::ceil(blaze::log(parameters.np))),parameters.lambda/(2*parameters.np), parameters.numberOfClusters,parameters.dimensions) * parameters.sFactor;
		this->parameters.NrOfNonCenterCoresetPoints = s_outer - this->parameters.numberOfClusters;
		this->parameters.NrOfNonCenterCoresetPointsNode = this->parameters.s - this->parameters.numberOfClustersNode;
		this->parameters.cutOff = this->parameters.s*this->parameters.cutOffFactor;
		this->parameters.insertionCutOff = this->parameters.s * this->parameters.insertionFactor;
		this->parameters.nrPointsInTree = n0;

		if (this->parameters.removeCutOffFactor != -1)
			this->parameters.removeCutOff = this->parameters.nrPointsInTree * this->parameters.removeCutOffFactor;

		assert(n0 > parameters.s);
		assert(n0 < 4*parameters.s);



		tree = new CoresetTree(parameters);

		//Insert all points one at a time
		for(size_t i = 0; i < firstSPoints->rows(); i++)
		{
			blaze::DynamicMatrix<double> coordinates = blaze::submatrix(*firstSPoints,i,3,1,firstSPoints->columns()-3);
			int index = firstSPoints->at(i,0);
			auto it = data.insert({static_cast<int>(firstSPoints->at(i,0)),{index,firstSPoints->at(i,2),coordinates}});
			if (i < firstSPoints->rows()/2)
				this->tree->getRoot()->getLChild()->insertPoint(&it.first->second);
			else
				this->tree->getRoot()->getRChild()->insertPoint(&it.first->second);
		}
		


		tree->setCurrentLeaf(nullptr);
		return true;
	}
	bool DynCoreset::insertPoint(int index, double weight, std::vector<double>& coordinates)
	{
		//Check if the index of the point is unique
		//if (this->data.count(index) > 0)
			//return false;

		auto it = this->data.insert({index,{index,weight,coordinates}});
		bool insertShortCut = this->tree->insertPoint(&it.first->second);

		this->updatesSinceLastRemoveRecalc++;

		this->runOuterCoreset();

		//Loop over all points that were marked as removed but were not yet removed.
		//If one of them is part of the coreset, delete it.
		for (const CoresetPoint* pointToRemove : PointsToRemoveIter)
		{
			auto pointInCoreset = std::find(this->coresetPoints.begin(), this->coresetPoints.end(), pointToRemove);

			if (pointInCoreset != this->coresetPoints.end())
			{
				int PointIndexInCoreset = pointInCoreset - this->coresetPoints.begin();
				this->coresetPoints.erase(pointInCoreset);
				this->coresetWeights.erase(this->coresetWeights.begin()+PointIndexInCoreset);
			}
		}

		this->updatesThisPhase++;
		if (updatesThisPhase >= this->parameters.np/8)
			this->goToNextPhase(true);
		return true;

	}
	bool DynCoreset::applyOperation(blaze::DynamicMatrix<double>& point)
	{
		if (point.at(0,1) == 1)
			insertPoint(point);
		else if (point.at(0,1) == 0)
			removePoint(point.at(0,0));
		else
			return false;

		if (this->updatesSinceLastRemoveRecalc >= this->parameters.removeCutOff)
		{
			if (!this->PointsToRemoveIter.empty())
			{
				this->tree->removePoints(this->PointsToRemoveIter);
				for (auto p : PointsToRemoveIter)
					data.erase(p->getIndex());
				this->PointsToRemoveIter.clear();
				this->runOuterCoreset();
			}
			this->updatesSinceLastRemoveRecalc = 0;
		}
		return true;
		
	}

	bool DynCoreset::insertPoint(blaze::DynamicMatrix<double>& point)
	{
		this->parameters.nrPointsInTree++;
		int index = point.at(0,0);

		//Check if the index of the point is unique
		if (this->data.count(index) > 0)
		{
			std::cout << "Double Point " << point << "\n";
			return false;
		}
			
		double weight = point.at(0,2);
		blaze::DynamicMatrix<double> coordinates = blaze::submatrix(point,0,3,1,point.columns()-3);
		auto it = this->data.insert({index,{index,weight,coordinates}});
		bool insertShortCut = this->tree->insertPoint(&it.first->second);



		//Loop over all points that were marked as removed but were not yet removed.
		//If one of them is part of the coreset, delete it.
		for (const CoresetPoint* pointToRemove : PointsToRemoveIter)
		{
			auto pointInCoreset = std::find(this->coresetPoints.begin(), this->coresetPoints.end(), pointToRemove);

			if (pointInCoreset != this->coresetPoints.end())
			{
				int PointIndexInCoreset = pointInCoreset - this->coresetPoints.begin();
				this->coresetPoints.erase(pointInCoreset);
				this->coresetWeights.erase(this->coresetWeights.begin()+PointIndexInCoreset);
			}
		}


		this->updatesThisPhase++;
		if (std::abs(updatesThisPhase) >= this->parameters.np/8)
			this->goToNextPhase(true);

		this->runOuterCoreset();
		return true;

	}

	bool DynCoreset::removePoint(int index)
	{

		
		this->parameters.nrPointsInTree--;
		this->updatesThisPhase--;
		boost::unordered_map<int,CoresetPoint>::iterator it = this->data.find(index);
		if (it == this->data.end())
			return false;
		auto inFinalCoreset = std::find(this->coresetPoints.begin(), this->coresetPoints.end(),&it->second);

		if (inFinalCoreset == this->coresetPoints.end() && std::abs(updatesThisPhase) < this->parameters.np/8)
		{
			this->updatesSinceLastRemoveRecalc++;
			this->PointsToRemoveIter.push_back(&it->second);
			return true;
		}
		this->PointsToRemoveIter.push_back(&it->second);

		this->tree->removePoints(this->PointsToRemoveIter);

		for (auto p : PointsToRemoveIter)
			data.erase(p->getIndex());
		this->PointsToRemoveIter.clear();
		this->updatesSinceLastRemoveRecalc = 0;

		if (std::abs(updatesThisPhase) >= this->parameters.np/8)
			this->goToNextPhase(false);

		this->runOuterCoreset();


		return true;

	}

	void DynCoreset::goToNextPhase(bool nrPointsIncreased)
	{
		this->updatesThisPhase = 0;
		if (nrPointsIncreased)
			this->parameters.np = 4*this->data.size();
		else
			this->parameters.np = this->parameters.np - (this->parameters.np/2);
		double s_outer = this->cardinalityOfCoreset(parameters.epsilon/(6*blaze::ceil(blaze::log(parameters.np))),parameters.lambda/(2*parameters.np), parameters.numberOfClusters,parameters.dimensions);
		this->parameters.s = this->cardinalityOfCoreset(parameters.epsilon/(6*blaze::ceil(blaze::log(parameters.np))),parameters.lambda/(2*parameters.np), parameters.numberOfClusters,parameters.dimensions) * parameters.sFactor;
		this->parameters.NrOfNonCenterCoresetPoints = s_outer - this->parameters.numberOfClusters;
		this->parameters.NrOfNonCenterCoresetPointsNode = this->parameters.s - this->parameters.numberOfClustersNode;
		this->parameters.cutOff = this->parameters.s*this->parameters.cutOffFactor;
		this->parameters.insertionCutOff = this->parameters.s * this->parameters.insertionFactor;
		this->tree->resetRefresh(nrPointsIncreased);

		std::cout << "\nNow s = " << parameters.s << "\n";
	}


	void DynCoreset::runOuterCoreset()
	{
		//Directly write the coreset of the tree root to the result vectors if it is smaller than the required size of the outer coreset
		if (tree->getRoot()->getData().size() <= this->parameters.sizeOfFinalCoreset || tree->getRoot()->getData().size() <= this->parameters.insertionCutOff)
		//if (tree->getRoot()->getData().size() <= this->parameters.sizeOfFinalCoreset)
		{
			this->coresetPoints = tree->getRoot()->getData();
			this->coresetWeights = tree->getRoot()->getWeights();
			std::for_each(coresetWeights.begin(), coresetWeights.end(), [this](double& w){w /= tree->getRoot()->getWeightsDenominator();});
			return;
		}
		std::shared_ptr<coresets::Coreset> results = this->staticCoreset.run(tree->getRoot()->getData(),this->data2Proxy,tree->getRoot()->getWeights(), this->weights2Proxy);

		//Delete the old centers, because the coreset algorithm calculated new ones. Then write points, centers and weights to the result vectors
		results->writePointsToVector(&coresetPoints,&centers);
		results->writeWeightsToVector(&coresetWeights);
		std::for_each(coresetWeights.begin(), coresetWeights.end(), [this](double& w){w /= tree->getRoot()->getWeightsDenominator();});

	}

	void DynCoreset::WriteCoresetToStream(std::ostream &out)
	{
		out << "Coreset:\n";
		for (int i = 0; i < this->coresetPoints.size(); i++)
		{
			out << this->coresetWeights.at(i)/this->tree->getRoot()->getWeightsDenominator() << "\n" << this->coresetPoints.at(i)->getCoordinates() << "\n\n";
		}
	}
}
