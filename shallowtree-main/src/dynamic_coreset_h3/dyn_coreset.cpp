#include "dyn_coreset.h"
#include "coreset_point.h"
#include <memory>
#include <vector>
#include <iostream>

namespace dynClustering
{
	bool DynCoreset::initTree(std::shared_ptr<blaze::DynamicMatrix<double>> firstSPoints)
	{
		int n0 = 2;
		this->np = 4*n0;

		double s_outer = this->cardinalityOfCoreset(epsilon/(6*blaze::ceil(blaze::log(np))),lambda/(2*np), numberOfClusters,dimensions);
		this->s = s_outer * sFactor;
		this->numberOfNotCenterCoresetPoints = s_outer - this->numberOfClusters;
		this->numberOfNotCenterCoresetPointsNode = s - this->numberOfClusters;
		this->cutOff = this->s*this->cutOffFactor;




		tree = new CoresetTree(numberOfClusters, dimensions, sizeOfFinalCoreset, &s, &numberOfNotCenterCoresetPointsNode, &cutOff, &np, this->epsilon, clusterFactor, maxPointsFactor);

		for(auto i = 0; i < firstSPoints->rows(); i++)
		{
			blaze::DynamicMatrix<double> point = blaze::submatrix(*firstSPoints,i,0,1,firstSPoints->columns());

			this->insertPoint(point, true);
		}
		


		return true;
	}
	bool DynCoreset::insertPoint(int index, double weight, std::vector<double>& coordinates)
	{
		//!!!!For debug
		if (this->data.count(index) > 0)
			return false;

		auto it = this->data.insert({index,{index,weight,coordinates}});
		this->tree->insertPoint(&it.first->second);

		this->runOuterCoreset();

		this->updatesThisPhase++;
		if (updatesThisPhase >= this->np/8)
			this->goToNextPhase(true);
		return true;

	}
	bool DynCoreset::applyOperation(blaze::DynamicMatrix<double>& point)
	{
		if (point.at(0,1) == 1)
			return insertPoint(point);
		else if (point.at(0,1) == 0)
			return removePoint(point.at(0,0));
		else
			return false;
		
	}

	bool DynCoreset::insertPoint(blaze::DynamicMatrix<double>& point, bool isInit)
	{
		int index = point.at(0,0);
		double weight = point.at(0,2);
		blaze::DynamicMatrix<double> coordinates = blaze::submatrix(point,0,3,1,point.columns()-3);

		auto it = this->data.insert({index,{index,weight,coordinates}});
		this->tree->insertPoint(&it.first->second, isInit);
		

		this->runOuterCoreset();

		this->updatesThisPhase++;
		if (std::abs(updatesThisPhase) >= this->np/8)
			this->goToNextPhase(true);
		return true;

	}

	bool DynCoreset::removePoint(int index)
	{
		auto it = this->data.find(index);
		if (it == this->data.end())
			return false;

		this->tree->removePoint(&it->second);

		this->runOuterCoreset();

		this->updatesThisPhase--;
		if (std::abs(updatesThisPhase) >= this->np/8)
			this->goToNextPhase(false);

		return true;

	}

	void DynCoreset::goToNextPhase(bool nrPointsIncreased)
	{
		this->updatesThisPhase = 0;
		if (nrPointsIncreased)
			this->np = 4*this->data.size();
		else
			this->np = this->np - (this->np/2);
		this->s = blaze::max(static_cast<int>(this->numberOfClusters), this->cardinalityOfCoreset(epsilon/(6*blaze::ceil(blaze::log(np))),lambda/(2*np), numberOfClusters,dimensions) * sFactor);

		std::cout << "\nNow s = " << s << "\n";
	}


	void DynCoreset::runOuterCoreset()
	{
		for(int i = centers.size()-1; i >= 0; i--)
			delete centers.at(i);
		centers.clear();

		if (tree->getRoot()->getResultData().size() <= this->sizeOfFinalCoreset || sizeOfFinalCoreset == cutOff)
		{
			coresetPoints = tree->getRoot()->getResultData();
			coresetWeights = tree->getRoot()->getResultWeights();
			return;
		}
		auto results = this->staticCoreset.run(tree->getRoot()->getResultData(),this->data2Proxy,tree->getRoot()->getResultWeights(), this->weights2Proxy);
		results->writePointsToVector(&coresetPoints,&centers);
		results->writeWeightsToVector(&coresetWeights);

	}

	void DynCoreset::WriteCoresetToStream(std::ostream &out)
	{
		out << "Coreset:\n";
		for (int i = 0; i < this->coresetPoints.size(); i++)
		{
			out << this->coresetWeights.at(i) << "\n" << this->coresetPoints.at(i)->getCoordinates() << "\n\n";
		}
	}
}
