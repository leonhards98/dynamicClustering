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
		if (this->lChild == nullptr && this->rChild == nullptr)
			point->setRelevantLeaf(this);
		return true;

	}

	void CoresetNode::changeLChild(std::shared_ptr<CoresetNode> newChild, CoresetNode* newParent)
	{
		lChild.reset();
		lChild = newChild;
		newChild->parent = newParent;
	}
	void CoresetNode::changeRChild(std::shared_ptr<CoresetNode> newChild, CoresetNode* newParent)
	{
		rChild.reset();
		rChild = newChild;
		newChild->parent = newParent;
	}


	bool CoresetNode::runStaticCoreset(CoresetPoint *pointToInsert, std::unordered_set<CoresetNode*>* newNodesForRemove)
	{


		if (this->isInternal())
		{
			int lWeight {this->getLChild()->weightDenominator};
			int rWeight {this->getRChild()->weightDenominator};
			
			//If the union of the points contained in the two children is small enough, copy all points and weights to
			//this node. If the weight denominators of the two children are not identical, set the weight denominator of this
			//node to the product of the two and multiply the weights of each child by the weight denominator of the other child
			if (this->getLChild()->getData().size() + this->getRChild()->getData().size() <= this->parameters.cutOff)
			{
				this->weights = this->getLChild()->weights;
				std::vector<double>::iterator border = this->weights.insert(this->weights.end(),this->getRChild()->weights.begin(), this->getRChild()->weights.end());


				if ((lWeight == 1) && (rWeight == 1))
					this->weightDenominator = 1;

				else if (lWeight == rWeight)
					this->weightDenominator = this->getRChild()->weightDenominator;
				else
				{
					this->weightDenominator = lWeight*rWeight;
					std::for_each(weights.begin(), border,[rWeight] (double &w) { w *= rWeight; });
					std::for_each(border, weights.end(),[lWeight] (double &w) { w *= lWeight; });

				}
				this->data = this->getLChild()->data;
				this->data.insert(this->data.end(),this->getRChild()->data.begin(), this->getRChild()->data.end());

				assert(data.size() == weights.size());

			}
			//If the number of points in the coreset is less than the insertionCutOff, just insert the point into
			//the coreset without any calculations using a weight of 1
			else if (this->getNrPoints() < this->parameters.insertionCutOff && pointToInsert != nullptr)
			{
				this->data.push_back(pointToInsert);
				this->weights.push_back(this->weightDenominator * pointToInsert->getInitialWeight());

			}
			//If the union of the points contained in the two children is to large, calcualte a coreset and save it to this node.
			//Procede with the weight denominators the same way as if no coreset is calculated
			else
			{
				//delete all centers, because new ones will be calculated. This function will then be called on all anchestors of
				//this node, ensuring that these centers are not needed elsewere
				for(int i = centers.size()-1; i >= 0; i--)
					delete centers.at(i);
				centers.clear();

				if ((lWeight == 1) && (rWeight == 1))
					this->weightDenominator = 1000;// blaze::pow(*np,cForRounding+1)*blaze::ceil(1/eps);

				else if (lWeight == rWeight)
					this->weightDenominator = this->getRChild()->weightDenominator;
				else
				{
					this->weightDenominator = lWeight*rWeight;
					this->getLChild()->weightDenominator = this->weightDenominator;
					this->getRChild()->weightDenominator = this->weightDenominator;

					std::for_each(this->getLChild()->weights.begin(), this->getLChild()->weights.end(), [rWeight] (double &w) { w *= rWeight; });
					std::for_each(this->getRChild()->weights.begin(), this->getRChild()->weights.end(), [lWeight] (double &w) { w *= lWeight; });

				}
				std::shared_ptr<coresets::Coreset> result = this->staticCoreset.run(this->getLChild()->getData(), this->getRChild()->getData(), this->getLChild()->getWeights(), this->getRChild()->getWeights());
				result->writePointsToVector(&data,&this->centers);
				result->writeWeightsToVector(&weights);

				int weightDenominator = this->weightDenominator;

				if ((lWeight == 1) && (rWeight == 1))
					std::for_each(this->weights.begin(), this->weights.end(), [weightDenominator] (double &w) { w = blaze::floor(w * weightDenominator); });

				else
					std::for_each(this->weights.begin(), this->weights.end(), [weightDenominator] (double &w) { w = blaze::floor(w); });

				assert(data.size() == weights.size());

				//If one coreset needs to be recalculated, all coresets in the anchestor nodes
				//also need to be recalculated -> pass a nullptr
				if (newNodesForRemove != nullptr)
					newNodesForRemove->insert(parent);
				else if (this->parent != nullptr)
					parent->runStaticCoreset(nullptr);
				return false;
			}
		}
		if (newNodesForRemove != nullptr)
			newNodesForRemove->insert(parent);
		else if (this->parent != nullptr)
			return parent->runStaticCoreset(pointToInsert);
		return true;


		
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
		this->centers.clear();
	}

}

//delete all centers, because new ones will be calculated. This function will then be called on all anchestors of
				//this node, ensuring that these centers are not needed elsewere
