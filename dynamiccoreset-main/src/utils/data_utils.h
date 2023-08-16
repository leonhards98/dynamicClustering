#pragma once

#include "blaze/math/dense/DynamicVector.h"
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <fstream>
#include <blaze/Math.h>
#include <vector>
#include <iostream>
#include <float.h>
#include <unordered_map>
#include "../dynamic_coreset/dyn_coreset.h"
#include "../coreset_old/coreset.hpp"

/**
 * @brief Parse the Points in the given file. Each line in the file should correspond to one point and its elements should be
 * seperated by a tab, whitespace or comma. The first element should be the unique index of the point (unsigned int), followed by
 * if the point should be inserted (1) or deleted (0). The next element should be the weight of the point, followed by the coordinates
 * @param filePath Path to the file containing the points
 * @param nr_of_points The number of points (rows in the file) that should be read
 * @param dimenstion The dimension of the points (number of elements that are coordinates)
 * @return A shared pointer to a (nr_of_points x dimenstion+3) matrix containing all elements of the file in the same order
 */
std::shared_ptr<blaze::DynamicMatrix<double>> parse(const std::string &filePath, int nr_of_points, int dimenstion)
{
	std::ifstream fin(filePath);
	if(!fin) throw std::runtime_error("File not found : " + filePath);
	std::string linebuffer;

	int line_nr = 0;
	auto data = std::make_shared<blaze::DynamicMatrix<double>>(nr_of_points, dimenstion+3);

	while(fin && getline(fin,linebuffer)){
		if(nr_of_points == line_nr) break;
		std::vector<std::string> tokens;
		boost::split(tokens, linebuffer, boost::is_any_of("\t ,;"));
		
		for (int i = 0; i < dimenstion+3; i++)
			data->at(line_nr,i) = stod(tokens[i]);
		line_nr++;
	}
	fin.close();
	return data;

}

/**
 * @brief Parses a file containg the ground truth (true cluster assignments) of the input points. Each row should represent the
 * assignment for one point where the order of the points is the same as in the file containg the points themselves. Each row
 * should contain one number that corresponds to the index of the cluster that the respective point is part of
 * @param filePath Path to the file containing the assignment
 * @param nr_of_points The number of points (rows in the file) that should be read
 * @param nr_colums The number of columns to read
 * @return A vector containing the indices of the clusters in the same order as they occred in the file
 */
std::shared_ptr<blaze::DynamicMatrix<double>> parseGroundTruth(const std::string &filePath, int nr_of_points, int nr_columns)
{
	auto data = std::make_shared<blaze::DynamicMatrix<double>>(nr_of_points,nr_columns);
	if (filePath == "")
		return data;

	std::ifstream fin(filePath);
	if(!fin) throw std::runtime_error("File not found : " + filePath);
	std::string linebuffer;

	int line_nr = 0;

	while(fin && getline(fin,linebuffer)){
		if(nr_of_points == line_nr) break;
		std::vector<std::string> tokens;
		boost::split(tokens, linebuffer, boost::is_any_of("\t ,;"));


		for (int c = 0; c < nr_columns; c++)
			data->at(line_nr,c) = std::stod(tokens[c]);

		line_nr++;
	}
	fin.close();
	return data;

}



/**
 * @brief Calculates the Rand-Index of two cluster assignments. Calculated by sampling 10 000 random pairs of points.
 * @param centers1 Coordintates of the centers of the first cluster assingment. A matrix where each row should contain the
 * coordinates of one center
 * @param centers2 Coordintates of the centers of the second cluster assingment. A matrix where each row should contain the
 * coordinates of one center
 * @param points Coordinates of the points that were clustered. Each row should contain the coordinates of one point. 
 * @return The Rand-Index of the two cluster assignments that can be created by assigning each point to the center that is closest to it
 */
double calcRandIndex(const blaze::DynamicMatrix<double> &centers1, const blaze::DynamicMatrix<double> &centers2, const blaze::DynamicMatrix<double> &points, const blaze::DynamicVector<double> &weights)
{

	int total {10000};
	double sumWeights {0};
	double agree {0};
	for (int i = 0; i < total; i++)
	{
		//Two random points
		int p1Index = rand()%points.rows();
		int p2Index = rand()%points.rows();
		auto p1 = blaze::row(points, p1Index);
		auto p2 = blaze::row(points, p2Index);

		//"minDistp1c1" is the distance between the point p1 and the points in centers1 that is closest to it
		double minDistp1c1 {DBL_MAX};
		double minDistp1c2 {DBL_MAX};
		double minDistp2c1 {DBL_MAX};
		double minDistp2c2 {DBL_MAX};

  		//"p1Closestc1" is the index of the point in centers1 that is closest to the point p1
		int p1Closestc1;
		int p2Closestc1;
		int p1Closestc2;
		int p2Closestc2;

		for (int j = 0; j < centers1.rows(); j++)
		{
			auto c = blaze::row(centers1,j);
			auto distp1c1 = blaze::norm(c-p1);
			auto distp2c1 = blaze::norm(c-p2);

			if (distp1c1 < minDistp1c1)
			{
				minDistp1c1 = distp1c1;
				p1Closestc1 = j;
			}
			if (distp2c1 < minDistp2c1)
			{
				minDistp2c1 = distp2c1;
				p2Closestc1 = j;
			}
		}

		for (int j = 0; j < centers2.rows(); j++)
		{
			auto c = blaze::row(centers2,j);
			auto distp1c2 = blaze::norm(c-p1);
			auto distp2c2 = blaze::norm(c-p2);

			if (distp1c2 < minDistp1c2)
			{
				minDistp1c2 = distp1c2;
				p1Closestc2 = j;
			}
			if (distp2c2 < minDistp2c2)
			{
				minDistp2c2 = distp2c2;
				p2Closestc2 = j;
			}
		}
		//Increase the counter if both assignments assign both points to the same or both points to different clusters
		if (p1Closestc1 == p2Closestc1 && p1Closestc2 == p2Closestc2)
			agree += weights.at(p1Index) + weights.at(p2Index);
		else if (p1Closestc1 != p2Closestc1 && p1Closestc2 != p2Closestc2)
			agree += weights.at(p1Index) + weights.at(p2Index);

		sumWeights += weights.at(p1Index) + weights.at(p2Index);


	}
	//Return the fraction of sampled points pairs for which both assignments assinged the same
	//or for which both assignments assinged different clusters
	return agree/sumWeights;

}

/**
 * @brief Calculates the Rand-Index of two cluster assignments. Calculated by sampling 10 000 random pairs of points.
 * @param centers1 Coordintates of the centers of the first cluster assingment. A matrix where each row should contain the
 * coordinates of one center
 * @param assignment2 A vector containing for each point the cluster it belongs to (integer)
 * @param points Coordinates of the points that were clustered. Each row should contain the coordinates of one point. 
 * @return The Rand-Index of the two cluster assignments that can be created by assigning each point to the center that is closest to it
 */
double calcRandIndexWithAssignment(const blaze::DynamicMatrix<double> &centers1, const std::vector<int> &assignment2, const blaze::DynamicMatrix<double> &points, const blaze::DynamicVector<double> &weights)
{
	if (assignment2.size() == 0)
		return 0;

	int total {10000};
	double sumWeights {0};
	double agree {0};
	for (int i = 0; i < total; i++)
	{
		int p1Index = rand()%points.rows();
		int p2Index = rand()%points.rows();
		auto p1 = blaze::row(points, p1Index);
		auto p2 = blaze::row(points, p2Index);

		//"minDistp1c1" is the distance between the point p1 and the points in centers1 that is closest to it
		double minDistp1c1 {DBL_MAX};
		double minDistp2c1 {DBL_MAX};

  		//"p1Closestc1" is the index of the point in centers1 that is closest to the point p1
		int p1Closestc1;
		int p2Closestc1;

		for (int j = 0; j < centers1.rows(); j++)
		{
			auto c = blaze::row(centers1,j);
			auto distp1c1 = blaze::norm(c-p1);
			auto distp2c1 = blaze::norm(c-p2);

			if (distp1c1 < minDistp1c1)
			{
				minDistp1c1 = distp1c1;
				p1Closestc1 = j;
			}
			if (distp2c1 < minDistp2c1)
			{
				minDistp2c1 = distp2c1;
				p2Closestc1 = j;
			}
		}

		//Increase the counter if both assignments assign both points to the same or both points to different clusters
		if (p1Closestc1 == p2Closestc1 && assignment2.at(p1Index) == assignment2.at(p2Index))
			agree += weights.at(p1Index) + weights.at(p2Index);
		else if (p1Closestc1 != p2Closestc1 && assignment2.at(p1Index) != assignment2.at(p2Index))
			agree += weights.at(p1Index) + weights.at(p2Index);

		sumWeights += weights.at(p1Index) + weights.at(p2Index);

	}
	//Return the fraction of sampled points pairs for which both assignments assinged the same
	//or for which both assignments assinged different clusters
	return agree/sumWeights;

}

/**
 * @brief Calculates the sum of the weighted squared cost of a given cluster assignment
 * @param points Coordinates of the points that were clustered. Each row should contain the coordinates of one point. 
 * @param centers Coordintates of the centers of the cluster assingment. A matrix where each row should contain the
 * coordinates of one center
 * @param weights Weights of the points. A vector where each row corresponds to the weight of the point that is stored at the
 * same position (row) in the "points" matrix
 * @return The total weighted sqaured cost of the given clustering
 */
double calcTotalSquaredCost(const blaze::DynamicMatrix<double> &points, const blaze::DynamicMatrix<double> & centers, const blaze::DynamicVector<double> &weights)
{
	blaze::DynamicVector<double> distancesSquared(points.rows(), DBL_MAX);

	for (int c = 0; c < centers.rows(); c++)
	{
		for (int p = 0; p < points.rows(); p++)
		{
			const double dist = blaze::sqrNorm(blaze::row(centers,c)-blaze::row(points,p));
			if (dist < distancesSquared.at(p))
				distancesSquared.at(p) = dist;

		}
	}
	return blaze::dot(distancesSquared,weights);
}


/**
 * @brief Calculates the sum of the weighted squared cost of a given cluster assignment
 * @param points Coordinates of the points that were clustered. Each row should contain the coordinates of one point. 
 * @param assignment A vector containing for each point the cluster it belongs to (integer)
 * @param weights Weights of the points. A vector where each row corresponds to the weight of the point that is stored at the
 * same position (row) in the "points" matrix
 * @param k Number of clusters
 * @return The total weighted sqaured cost of the given clustering
 */
double calcTotalSquaredCostWithAssignment(const blaze::DynamicMatrix<double> &points, const std::vector<int> &assignment, const blaze::DynamicVector<double> &weights, int k)
{
	blaze::DynamicMatrix<double> newCenters(k,points.columns(),0);
	std::vector<int> clusterMemberCounts(k,0);

	for (int r = 0; r < points.rows(); r++)
	{
		// The indices from the input files start from 1, but we need indices from 0 to k here
		const int cIndex = assignment.at(r) - 1;

		for (int c = 0; c < points.columns(); c++)
			newCenters.at(cIndex,c) += (points.at(r,c) * weights.at(r));

		clusterMemberCounts.at(cIndex) += weights.at(r);
	}

	for (int c = 0; c < k; c++)
	{
		const auto count = std::max<size_t>(1, clusterMemberCounts[c]);
		blaze::row(newCenters, c) /= count;
	}

	//std::cout << newCenters << "\n";

	return calcTotalSquaredCost(points,newCenters,weights);


}

/**
 * @brief Writes the coreset points calculated by the dynamic coreset algorithm into a points matrix and a weights vector
 * @param coresetPoints Matrix into which the points will be written. One row per point
 * @param coresetWeights Vector containing the weights of the points. Same order as in the point matrix
 * @param coreset coreset calcualted by the dynamic algorithm
 * @param dimensions Dimension of the points in the dataset
 * @return (void)
 */
void writeDynCoresetToMatrix(blaze::DynamicMatrix<double>& coresetPoints, blaze::DynamicVector<double> &coresetWeights, const dynClustering::DynCoreset &coreset, int dimensions)
{
		int size = coreset.getCoresetSize();

		coresetPoints.resize(size, dimensions);
		coresetWeights.resize(size);

		// Writes the output of the dynamic algorithm to "coresetPoints" and "coresetWeights"
		for (int r = 0; r < size; r++)
		{
			coresetWeights.at(r) = coreset.getCoresetWeights().at(r);
			for (int j = 0; j < dimensions; j++)
			{
				coresetPoints.at(r, j) = coreset.getCoresetPoints().at(r)->getCoordinates().at(j);
			}
		}
}

/**
 * @brief Writes the coreset points calculated by the dynamic coreset algorithm into a points matrix and a weights vector
 * @param coresetPoints Matrix into which the points will be written. One row per point
 * @param coresetWeights Vector containing the weights of the points. Same order as in the point matrix
 * @param Points All points used to calcaulte the static coreset
 * @param CoresetResult Result from the static coreset algorithm
 * @param dimensions Dimension of the used datapoints
 * @return (void)
 */
void writeStatCoresetToMatrix(blaze::DynamicMatrix<double>& coresetPoints, blaze::DynamicVector<double> &coresetWeights, const blaze::DynamicMatrix<double> &staticPoints, const std::shared_ptr<coresetsStatic::Coreset> staticCoresetResult, int dimensions)
{
			int size = staticCoresetResult->size();
			coresetPoints.resize(size, dimensions);
			coresetWeights.resize(size);
			for (int i = 0; i < size; i++)
			{
				coresetWeights.at(i) = staticCoresetResult->at(i)->Weight;
				if (staticCoresetResult->at(i)->IsCenter != true)
				{
					for (int j = 0; j < dimensions; j++)
						coresetPoints.at(i, j) = staticPoints.at(staticCoresetResult->at(i)->Index, j);
				}
				else
				{
					const blaze::DynamicVector<double> center = staticCoresetResult->getCenterCoordinates(staticCoresetResult->at(i)->Index);
					for (int j = 0; j < dimensions; j++)
						coresetPoints.at(i, j) = center.at(j);
				}
			}

}

/**
 * @brief Appies the next operation (point insertion/removal) to the vector containing all current points.
 * @param Points All points currently used by the algorithms which will be modified.
 * @param Weights Weights of the points in the same order which will be modified.
 * @param point new point to insert or delete.
 * @param dimensions dimension of the used points.
 * @return 
 */
double updateTotalDataset(blaze::DynamicMatrix<double> &staticPoints, blaze::DynamicVector<double> &staticWeights, blaze::DynamicVector<double> &staticIndices, const blaze::DynamicMatrix<double> &point, const int dimensions)
{

		int insertionCounter {0};
		if (point.at(0, 1) == 1)
		{
			staticPoints.resize(staticPoints.rows() + 1, dimensions - 3);
			for (int j = 0; j < dimensions - 3; j++)
				staticPoints.at(staticPoints.rows() - 1, j) = point.at(0, j + 3);
			staticWeights.resize(staticWeights.size() + 1);
			staticIndices.resize(staticWeights.size() + 1);
			staticWeights.at(staticPoints.rows() - 1) = point.at(0, 2);
			staticIndices.at(staticPoints.rows() - 1) = point.at(0, 0);
			insertionCounter++;
		}
		// Removes the current point from "staticPoints", "staticWeights" and "staticGroundTruth"
		// Does this by clearing the matrices and the inserting all points except the one that should be removed
		else if (point.at(0, 1) == 0)
		{
			// intermediate storage of the staticPoints
			blaze::DynamicMatrix<double> zw = staticPoints;
			// intermediate storage of the staticWeights
			blaze::DynamicVector<double> zwWeight = staticWeights;
			// intermediate storage of the staticIndices
			blaze::DynamicVector<double> zwIndices = staticIndices;
			staticPoints.clear();
			staticWeights.clear();
			staticIndices.clear();
			staticPoints.resize(zw.rows() - 1, dimensions - 3);
			staticWeights.resize(zwWeight.size() - 1);
			staticIndices.resize(zwIndices.size() - 1);
			int counter{0};
			int j;
			for (j = 0; j < zw.rows(); j++)
			{
				if (zwIndices[j] == point.at(0,0))
				{
					/*
					bool correct = true;
					for (int i = 0; i < zw.columns(); i++)
						if (blaze::row(zw, j)[i] != blaze::row(blaze::submatrix(point, 0, 3, 1, point.columns() - 3), 0)[i])
							correct = false;
					std::cout << std::setprecision (20) << blaze::row(zw, j) << blaze::row(blaze::submatrix(point, 0, 3, 1, point.columns() - 3), 0) << std::endl;
					if (correct == true)
					*/
					continue;
				}

				staticWeights.at(counter) = zwWeight.at(j);
				staticIndices.at(counter) = zwIndices.at(j);
				for (int q = 0; q < zw.columns(); q++)
					staticPoints.at(counter, q) = zw.at(j, q);

				counter++;
			}
			assert(counter == (j - 1));
		}
		return insertionCounter;

}




















