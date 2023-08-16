#pragma once

#include <algorithm>
#include <iostream>
#include <memory>
#include <random>
#include <vector>
#include <float.h>

#include <blaze/Math.h>
#include <boost/array.hpp>
#include "../dynamic_coreset/coreset_point.h"

namespace utils
{
	class IndexFinder
	{
        const std::vector<dynClustering::CoresetPoint*>* X1;
        const std::vector<dynClustering::CoresetPoint*>* X2;
        int size;

		public:

        IndexFinder(const std::vector<dynClustering::CoresetPoint*>* X1, const std::vector<dynClustering::CoresetPoint*>* X2) : X1(X1), X2(X2), size(X1->size()) {}
        IndexFinder() {}

        void init(const std::vector<dynClustering::CoresetPoint*>* X1, const std::vector<dynClustering::CoresetPoint*>* X2)
        {
            this->X1 = X1;
            this->X2 = X2;
            this->size = X1->size();
        }

		dynClustering::CoresetPoint* getPoint(int index)
		{
			if (index < size)
				return X1->at(index);

			return X2->at(index-size);
		}



	};

    class L2NormCalculator
    {
    private:
        const std::vector<dynClustering::CoresetPoint*> X1;
        const std::vector<dynClustering::CoresetPoint*> X2;
		IndexFinder i;
        
    public:
        L2NormCalculator(const std::vector<dynClustering::CoresetPoint*> &X1, const std::vector<dynClustering::CoresetPoint*> &X2, IndexFinder &i) : X1(X1), X2(X2), i(i)
        {
        }


        double calc(size_t xIndex, size_t yIndex)
        {
            if (xIndex == yIndex)
                return 0.0;
            if (i.getPoint(xIndex)->getCoordinates() == i.getPoint(yIndex)->getCoordinates())
                return 0.0;

            
			double dotProd = blaze::dot(i.getPoint(xIndex)->getCoordinates(),i.getPoint(yIndex)->getCoordinates());

            double squaredNorm = i.getPoint(xIndex)->getSumOfSquares() + i.getPoint(yIndex)->getSumOfSquares() - 2 * dotProd;

            assert (squaredNorm < DBL_MAX);

            return  squaredNorm;
        }
    };
}
