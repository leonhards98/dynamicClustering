#include "coreset_node.h"
#include <memory>
#include <algorithm>
#include <cassert>

namespace dynClustering
{
	bool CoresetNode::insertPoint(CoresetPoint* point)
	{
		this->data.push_back(point);
		this->weights.push_back(point->getInitialWeight());
		point->setRelevantLeaf(this);
		return true;

	}



	void CoresetNode::runStaticCoreset()
	{

		for(int i = centers.size()-1; i >= 0; i--)
			delete centers.at(i);
		centers.clear();

		if (this->data.size() <= *this->cutOff)
		{
			this->resultWeights = this->weights;
			this->resultData = this->data;

			assert(resultData.size() == resultWeights.size());

		}
		else
		{
			std::vector<CoresetPoint*> proxyData {};
			std::vector<double> proxyWeight {};
			auto result = this->staticCoreset.run(this->data, proxyData, this->weights, proxyWeight);
			result->writePointsToVector(&resultData,&this->centers);
			result->writeWeightsToVector(&resultWeights);

			assert(resultData.size() == resultWeights.size());
		}

	}


	void CoresetNode::writePointsToStream(std::ostream &out)
	{
		for (int i = 0; i < this->data.size(); i++)
		{
			std::cout << "weight: " << this->weights.at(i) << "\n";
			std::cout << this->data.at(i)->getCoordinates() << "\n";

		}
	}

	void CoresetNode::writeIndicesToStream(std::ostream &out)
	{
		out << "; ";
		for (int i = 0; i < this->data.size(); i++)
		{
			out << this->data.at(i)->getIndex() << ", ";

		}
	}

	void CoresetNode::erasePoint(int i)
	{
		this->data.erase(this->data.begin()+i);
		this->weights.erase(this->weights.begin()+i);

		assert(data.size() == weights.size());
	}

	void CoresetNode::appendWeights(const std::vector<double> &newWeights)
	{
		this->weights.insert(this->weights.end(),newWeights.begin(),newWeights.end());
	}

	void CoresetNode::appendPoints(const std::vector<CoresetPoint*> &newPoints)
	{
		for (int i = 0; i < newPoints.size(); i++)
		{
			newPoints.at(i)->setRelevantLeaf(this);
		}
		this->data.insert(this->data.end(),newPoints.begin(),newPoints.end());
	}

	void CoresetNode::clearData()
	{
		this->weights.clear();
		this->data.clear();

		for(int i = centers.size()-1; i >= 0; i--)
			delete centers.at(i);
		this->centers.clear();
	}

}
