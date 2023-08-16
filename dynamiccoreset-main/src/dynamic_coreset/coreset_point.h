#pragma once

#include <blaze/Math.h>
#include <vector>
#include <iostream>

namespace dynClustering
{

	class CoresetNode;
	/**
	 * @brief Represets a single Point with an index, weight and coordinates. The weight in this class is always the original weight
	 * and will not change or be the weight of the point in the final coreset
	 */
	class CoresetPoint
	{
		/**
		 * Weight of the point
		 */
		double weight;

		/**
		 * Index of the point. This should be unique for every point
		 */
		int index;

		/**
		 * The sum of squares of each dimension of the coordinates. Used to calculate distances faster
		 */
		double sumOfSquared;

		/**
		 * The leaf of the coreset tree in which this point is stored. A point is always stored in exactly one leaf, but can be
		 * part of an arbitratry number of internal nodes.
		 */
		CoresetNode *relevantLeaf;

		/**
		 * The coordinates of the point
		 */
		blaze::DynamicVector<double> coordinates;

	public:
		/**
		 * @brief Constructor. Also calculates the sum the squares of the coordinates
		 * @param index The unique and identifying index of the point.
		 * @param weight The original weight of the point
		 * @param coordinates The coordinates of the point
		 * @return
		 */
		CoresetPoint(int index, double weight, blaze::DynamicMatrix<double> &coordinates)
		{
			sumOfSquared = 0;
			this->weight = weight;
			this->index = index;
			this->coordinates.resize(coordinates.columns());
			for (int i = 0; i < coordinates.columns(); i++)
			{
				this->coordinates.at(i) = coordinates.at(0, i);
				sumOfSquared += coordinates.at(0, i) * coordinates.at(0, i);
			}
			this->relevantLeaf = nullptr;
		}

		/**
		 * @brief Constructor. Also calculates the sum the squares of the coordinates
		 * @param index The unique and identifying index of the point.
		 * @param weight The original weight of the point
		 * @param coordinates The coordinates of the point
		 * @return
		 */
		CoresetPoint(int index, double weight, blaze::DynamicVector<double> &coordinates)
		{
			sumOfSquared = 0;
			this->weight = weight;
			this->index = index;
			this->coordinates = coordinates;
			for (int i = 0; i < coordinates.size(); i++)
			{
				sumOfSquared += coordinates.at(i) * coordinates.at(i);
			}
			this->relevantLeaf = nullptr;
		}

		/**
		 * @brief Constructor. Also calculates the sum the squares of the coordinates
		 * @param index The unique and identifying index of the point.
		 * @param weight The original weight of the point
		 * @param coordinates The coordinates of the point
		 * @return
		 */
		CoresetPoint(int index, double weight, std::vector<double> &coordinates)
		{
			sumOfSquared = 0;
			this->weight = weight;
			this->index = index;
			this->coordinates.resize(coordinates.size());
			for (int i = 0; i < coordinates.size(); i++)
			{
				this->coordinates.at(i) = coordinates.at(i);
				sumOfSquared += coordinates.at(i) * coordinates.at(i);
			}
			this->relevantLeaf = nullptr;
		}

		void setRelevantLeaf(CoresetNode *relevantLeaf)
		{
			this->relevantLeaf = relevantLeaf;
		}
		void setIndex(int index)
		{
			this->index = index;
		}
		const int getIndex()
		{
			return this->index;
		}
		CoresetNode *getRelevantLeaf()
		{
			return relevantLeaf;
		}
		const blaze::DynamicVector<double> &getCoordinates()
		{
			return coordinates;
		}
		double getInitialWeight()
		{
			return weight;
		}
		double getSumOfSquares()
		{ 
			return sumOfSquared;
		}
	};
}
