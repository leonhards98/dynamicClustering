#pragma once

#include <blaze/Math.h>
#include <vector>
#include <iostream>


namespace dynClustering
{

	class CoresetNode;
	class CoresetPoint
	{
		double weight;
		int index;
		double sumOfSquared;
		CoresetNode* relevantLeaf;
		blaze::DynamicVector<double> coordinates;

		public:

		CoresetPoint (int index, double weight,blaze::DynamicMatrix<double>& coordinates) 
		{
			sumOfSquared = 0;
			this->weight = weight;
			this->index = index;
			this->coordinates.resize(coordinates.columns());
			for (int i = 0; i < coordinates.columns(); i++)
			{
				this->coordinates.at(i) = coordinates.at(0,i);
				sumOfSquared += coordinates.at(0,i)*coordinates.at(0,i);
			}
			this->relevantLeaf = nullptr;
		}
		CoresetPoint (int index, double weight,blaze::DynamicVector<double>& coordinates) 
		{
			sumOfSquared = 0;
			this->weight = weight;
			this->index = index;
			this->coordinates = coordinates;
			for (int i = 0; i < coordinates.size(); i++)
			{
				sumOfSquared += coordinates.at(i)*coordinates.at(i);
			}
			this->relevantLeaf = nullptr;
		}
		CoresetPoint (int index, double weight,std::vector<double>& coordinates) 
		{
			sumOfSquared = 0; 
			this->weight = weight;
			this->index = index;
			this->coordinates.resize(coordinates.size());
			for (int i = 0; i < coordinates.size(); i++)
			{
				this->coordinates.at(i) = coordinates.at(i);
				sumOfSquared += coordinates.at(i)*coordinates.at(i);
			}
			this->relevantLeaf = nullptr;
		}
		void setRelevantLeaf (CoresetNode* relevantLeaf)
		{
			this->relevantLeaf = relevantLeaf;
		}
		const int getIndex()
		{
			return this->index;
		}
		CoresetNode* getRelevantLeaf()
		{
			return relevantLeaf;
		}
		const blaze::DynamicVector<double>& getCoordinates () {return coordinates;}
		double getInitialWeight () {return weight;}
		double getSumOfSquares () {return sumOfSquared;}
		
		

	};
}
