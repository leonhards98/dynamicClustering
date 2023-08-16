#include <algorithm>
#include <blaze/math/dense/DynamicMatrix.h>
#include <blaze/math/views/Forward.h>
#include <iostream>
#include <memory>
#include <ostream>
#include <vector>
#include <string>
#include <chrono>
#include <unordered_map>

#include "dynamic_coreset/dyn_coreset.h"
#include "coreset/sensitivity_sampling.hpp"
#include "utils/data_utils.h"
#include "utils/CLI11.hpp"
#include "utils/dataStruct.h"
#include "blaze/Math.h"

#include "coreset_old/kmeans.hpp"
#include "coreset_old/sensitivity_sampling.hpp"

//Uncomment this to also run the ILP as a baseline
//#include "ILPBaseLine/baseline.h"

using namespace dynClustering;
/**
 * @brief Uses the theoretical Formula of the static coreset agorithm to determine
 * the output size of the coreset required to the the desired accuracy of lambda and epsilon.
 * If no epsilon and lambda is given, uses the fixed size given by the global variable "sizeOfFinalCoreset"
 * @param epsilon
 * @param lambda
 * @param k Number of centers
 * @param d Dimesion
 * @return required coreset size
 */
int cardinalityOfCoreset(double epsilon, double lambda, int k, int d);

//double checkKMeansError(DynCoreset &coreset, blaze::DynamicMatrix<double> &points, blaze::DynamicVector<double> &weights, int k, int d, int nrPointsInFile);

// Golbal variables for the size of the final coreset in case lambda and epsilon are not used
int sizeOfFinalCoreset;
bool useFinalSize;

int main(int argc, char *argv[])
{
	srand(time(nullptr));
	
	// Variables for commandline arguments
	std::string epsilon_str;
	std::string lambda_str;
	std::string sizeOfFinalCoresest_str;
	std::string infile;
	std::string outfilePath;
	std::string dimension_str;
	std::string nrOfPointsInFile_str;
	std::string nrClusters_str;
	std::string n0_str;
	std::string clusterFactor_str;
	std::string sFactor_str;
	std::string cutOffFactor_str;
	std::string skip_str;
	std::string insertionFactor_str;
	std::string removeCutOffFactor_str;
	std::string removeCutOffFixed_str;
	std::string PTimeFactor_str;


	// Parsing commandline arguments
	CLI::App inputParser("CommandLineArguments");
	inputParser.add_option("-f", infile, "File containing the input data")->required();
	inputParser.add_option("-o", outfilePath, "File to write the time measurements into")->required();
	inputParser.add_option("-l", lambda_str, "Epsilon for the epsilon-coreset. If not set, the final size of the coreset needs to be provided (-c)");
	inputParser.add_option("-e", epsilon_str, "Lambda. If not set, the final size of the coreset needs to be provided (-c)");
	inputParser.add_option("-c", sizeOfFinalCoresest_str, "Size of the final Coreset. If not set, need to provide value for lambda (-l) and epsilon (-e)");
	inputParser.add_option("-k", nrClusters_str, "Number of points to read from file. If empty, read all points")->required();
	inputParser.add_option("-d", dimension_str, "Dimesion of the points the in input file. If empty, read all dimensions")->required();
	inputParser.add_option("-i", n0_str, "Number of Points used to initalize the tree. If empty, it is set to 2*s");
	inputParser.add_option("-n", nrOfPointsInFile_str, "Number of points to read from file.")->required();
	inputParser.add_option("-y", clusterFactor_str, "Factor used to determine the number of clusters for the inner coresets. Default 1");
	inputParser.add_option("-s", sFactor_str, "Factor used to multiply the cutoff point s, that is also the size, for inner coreset algorithms. Default 1");
	inputParser.add_option("-t", cutOffFactor_str, "Factor used to multipy the size of the inner coresets (including sFactor)to get the cutoff Point for calculating a static coreset. Default 1");
	inputParser.add_option("-x", insertionFactor_str, "Factor used to multiply the size of inner coresets. Only recalculate if inner coresets would get bigger than that");
	inputParser.add_option("-v", removeCutOffFactor_str, "Remove cut-off Delta given in factor of all points in the dataset");
	inputParser.add_option("--vfixed", removeCutOffFixed_str, "Remove cut-off Delta given as an absolute number");
	inputParser.add_option("-p", PTimeFactor_str, "Factor used to multiply with the runtime of the dynamic algorithm to get the time allowed for the ILP and MIQP");
	inputParser.add_option("--skip", skip_str, "Skip for output");

	CLI11_PARSE(inputParser, argc, argv);

	if ((epsilon_str == "" || lambda_str == "") && sizeOfFinalCoresest_str == "")
	{
		std::cout << "Need to provide either lambda and epsilon or provide the number of points in the final coreset";
		return 0;
	}
	if (removeCutOffFactor_str != "" && removeCutOffFixed_str != "")
		std::cout << "Need to give only a relative or a abolute cutoff but not both";

	//Store the provided input

	::useFinalSize = sizeOfFinalCoresest_str == "" ? false : true;
	double epsilon{epsilon_str == "" ? 0 : std::stod(epsilon_str)};
	double lambda{lambda_str == "" ? 0 : std::stod(lambda_str)};
	int skip{skip_str == "" ? 20 : std::stoi(skip_str)};
	int nrClusters{std::stoi(nrClusters_str)};
	int dimensions{std::stoi(dimension_str)};
	int n0{0};
	int nrOfPointsInFile{std::stoi(nrOfPointsInFile_str)};
	::sizeOfFinalCoreset = sizeOfFinalCoresest_str == "" ? 0 : std::stoi(sizeOfFinalCoresest_str);
	int nrPointsFinalCoreset{cardinalityOfCoreset(epsilon / 3, lambda / 2, nrClusters, dimensions)};
	double clusterFactor{clusterFactor_str == "" ? 1 : std::stod(clusterFactor_str)};
	double sFactor{sFactor_str == "" ? 1 : std::stod(sFactor_str)};
	double cutOffFactor{cutOffFactor_str == "" ? 1 : std::stod(cutOffFactor_str)};
	double insertionFactor{insertionFactor_str == "" ? 1 : std::stod(insertionFactor_str)};
	double removeCutOffFactor{removeCutOffFactor_str == "" ? -1 : std::stod(removeCutOffFactor_str)/100};
	double removeCutOffFixed{removeCutOffFixed_str == "" ? 0 : std::stod(removeCutOffFixed_str)};
	double PTimeFactor{PTimeFactor_str == "" ? 1 : std::stod(PTimeFactor_str)};

	// If no start size is given, calcualte it s.t. it is as small as possible while still requiring coreset calculations in the tree. 
	if (n0_str == "")
	{
		if (::useFinalSize)
			n0 = (::sizeOfFinalCoreset + 1)*sFactor;
		else
		{
			int maxN0{50000};
			for (n0 = 10; n0 <= maxN0; n0 += 10)
			{
				int np = 4 * n0;
				int s = cardinalityOfCoreset(epsilon / (6 * blaze::ceil(blaze::log(np))), lambda / (2 * np), nrClusters, dimensions) * sFactor;
				if (n0 > s && n0 < 4 * s)
					break;
			}
			assert(n0 < maxN0);
		}
	}
	else
		n0 = std::stoi(n0_str);

	//Struct containing all parameters for the dynamic algorithm
	paramStruct parameters(nrClusters,epsilon,lambda,dimensions,clusterFactor,sFactor,cutOffFactor,sizeOfFinalCoreset, -1, -1);
	paramStruct parametersOpt(nrClusters,epsilon,lambda,dimensions,clusterFactor,sFactor,cutOffFactor,sizeOfFinalCoreset, insertionFactor, removeCutOffFactor, removeCutOffFixed);
	
	// Writes header to the output file
	std::ofstream outFile(outfilePath);
	outFile << "# Input: " << infile << ", Nr of Points: " << nrOfPointsInFile << "\n";
	outFile << "# k: " << nrClusters << ", d: " << dimensions << ", Nr of final Coreset Points: " << nrPointsFinalCoreset << "\n";
	outFile << "# epsilon: " << epsilon << ", lambda: " << lambda << ", n0: " << n0 << "\n";
	outFile << "# sFactor: " << sFactor << ", clusterFactor: " << clusterFactor << ", insertionFactor: " << insertionFactor << ", cutOffFactor: " << cutOffFactor << "\n";
	outFile << "# Operation DC_t DK_t D_c DC_q DK_q OC_t OK_t O_c OC_q OK_q SC_t SK_t S_c SC_q SK_q R_t T_c RC_q RK_q T_t T_c ILP_t ILP_c MIQP_t MIQP_c" << std::endl;

	std::shared_ptr<blaze::DynamicMatrix<double>> data = parse(infile, nrOfPointsInFile, dimensions);

	//Data points used to initialze the dynamic coreset algorithm
	std::shared_ptr<blaze::DynamicMatrix<double>> initData = std::make_shared<blaze::DynamicMatrix<double>>(blaze::submatrix(*data, 0, 0, n0, data->columns()));


	DynCoreset coreset(parameters, cardinalityOfCoreset);
	DynCoreset coresetOpt(parametersOpt, cardinalityOfCoreset);
	clusteringStatic::KMeans staticKMeans(nrClusters, true);

	// Contain the output of the dynamic or static coreset algorithm
	blaze::DynamicMatrix<double> coresetPoints;
	blaze::DynamicVector<double> coresetWeights;
	blaze::DynamicMatrix<double> coresetPointsOpt;
	blaze::DynamicVector<double> coresetWeightsOpt;
	blaze::DynamicMatrix<double> coresetPointsStat;
	blaze::DynamicVector<double> coresetWeightsStat;
	// Nr of points in the coreset
	int size;

	coreset.initTree(initData);
	coresetOpt.initTree(initData);

	// Sum of the runtimes for all dynamic coreset calcualtions since statistics were last calculated
	std::chrono::microseconds dynamicCoresetTime{0};
	// Sum of the runtimes for all dynamic coreset calcualtions since statistics were last calculated
	std::chrono::microseconds dynamicCoresetTimeOpt{0};
	// Runtimes for last static coreset calcualtion
	std::chrono::microseconds staticCoresetTime{0};
	// Sum of the runtimes for all calcualtions of KMeans usind dynamic coreset since statistics were last calculated
	std::chrono::microseconds dynamicKMeansTime{0};
	// Sum of the runtimes for all calcualtions of KMeans usind dynamic coreset since statistics were last calculated
	std::chrono::microseconds dynamicKMeansTimeOpt{0};
	// Runtimes for last KMeans calculation using static coreset
	std::chrono::microseconds staticKMeansTime{0};
	// Average runtime of selecting a random coreset and applying KMeans
	std::chrono::microseconds randomTime{0};
	// Runtime for the ILP in seconds
	double ILPTime {0};
	// Runtime for the MIQP in seconds
	double MIQPTime {0};


	// Contain the raw points/weights that are currently used in the dynamic and static algorithm
	blaze::DynamicMatrix<double> staticPoints;
	blaze::DynamicVector<double> staticWeights;
	blaze::DynamicVector<double> staticIndices;
	staticPoints = blaze::submatrix(*data, 0, 3, n0, data->columns() - 3);
	staticWeights.resize(n0);
	staticIndices.resize(n0);
	for (int j = 0; j < n0; j++)
	{
		staticWeights.at(j) = data->at(j, 2);
		staticIndices.at(j) = data->at(j, 0);
	}

	//Nr of points that have been inserted
	int insertionCounter {n0};

	for (int i = n0; i < nrOfPointsInFile; i++)
	{
		if (i % 1000 == 0)
			std::cout << "Now i = " << i << std::endl;

		// next Point to be inserted or deleted
		blaze::DynamicMatrix<double> point = blaze::submatrix(*data, i, 0, 1, data->columns());
		coresetsStatic::SensitivitySampling staticCoreset(nrClusters, cardinalityOfCoreset(epsilon, lambda, nrClusters, dimensions)-nrClusters);

		// Run the dynamic algorithm
		auto startDynamicCoreset = std::chrono::high_resolution_clock::now();
		coreset.applyOperation(point);
		auto stopDynamicCoreset = std::chrono::high_resolution_clock::now();

		// Run the dynamic algorithm
		auto startDynamicCoresetOpt = std::chrono::high_resolution_clock::now();
		coresetOpt.applyOperation(point);
		auto stopDynamicCoresetOpt = std::chrono::high_resolution_clock::now();


		std::shared_ptr<clusteringStatic::ClusteringResult> staticResult;
		std::shared_ptr<clusteringStatic::ClusteringResult> dynamicResult;
		std::shared_ptr<clusteringStatic::ClusteringResult> dynamicResultOpt;

		writeDynCoresetToMatrix(coresetPoints,coresetWeights,coreset, dimensions);


		// Runs the KMeans on the coreset resulting for the dynamic coreset algorithm
		auto startDynamicKMeans = std::chrono::high_resolution_clock::now();
		dynamicResult = staticKMeans.run(coresetPoints, coresetWeights, false);
		auto stopDynamicKMeans = std::chrono::high_resolution_clock::now();


		writeDynCoresetToMatrix(coresetPointsOpt,coresetWeightsOpt,coresetOpt, dimensions);


		// Runs the KMeans on the coreset resulting for the dynamic coreset algorithm
		auto startDynamicKMeansOpt = std::chrono::high_resolution_clock::now();
		dynamicResultOpt = staticKMeans.run(coresetPointsOpt, coresetWeightsOpt, false);
		auto stopDynamicKMeansOpt = std::chrono::high_resolution_clock::now();

		// Inserts the current point from "staticPoints", "staticWeights"
		updateTotalDataset(staticPoints, staticWeights, staticIndices, point, data->columns());

		//Use the static coreset
		//if (i % skip == 0 && i >= nrOfPointsInFile - 300 * skip)
		if (i % skip == 0 && i >= nrOfPointsInFile - 625 * skip)
		//if (i % skip == 0)
		{

			// Runs the static coreset algorithm on all points
			auto startStaticCoreset = std::chrono::high_resolution_clock::now();
			std::shared_ptr<coresetsStatic::Coreset> staticCoresetResult = staticCoreset.run(staticPoints, staticWeights);
			auto stopStaticCoreset = std::chrono::high_resolution_clock::now();


			// Writes the output of the static algorithm to "coresetPoints" and "coresetWeights"
			writeStatCoresetToMatrix(coresetPointsStat,coresetWeightsStat,staticPoints,staticCoresetResult,dimensions);

			// Runs the KMeans on the coreset resulting for the static coreset algorithm
			auto startStaticKMeans = std::chrono::high_resolution_clock::now();
			staticResult = staticKMeans.run(coresetPointsStat, coresetWeightsStat, false);
			auto stopStaticKMeans = std::chrono::high_resolution_clock::now();

			staticCoresetTime += duration_cast<microseconds>(stopStaticCoreset - startStaticCoreset);
			staticKMeansTime += duration_cast<microseconds>(stopStaticKMeans - startStaticKMeans);

		}
		dynamicCoresetTime += duration_cast<microseconds>(stopDynamicCoreset - startDynamicCoreset);
		dynamicKMeansTime += duration_cast<microseconds>(stopDynamicKMeans - startDynamicKMeans);
		
		dynamicCoresetTimeOpt += duration_cast<microseconds>(stopDynamicCoresetOpt - startDynamicCoresetOpt);
		dynamicKMeansTimeOpt += duration_cast<microseconds>(stopDynamicKMeansOpt - startDynamicKMeansOpt);


		// Every skip steps, calculate statistics and write them to the output file
		//if (i % skip == 0 && i >= nrOfPointsInFile - 300 * skip)
		if (i % skip == 0 && i >= nrOfPointsInFile - 625 * skip)
		//if (i % skip == 0)
		{

			// Total cost when using KMeans++ on all Points
			double CostWithoutCoresetWithoutLoyd{0};
			// Total cost when using KMeans on all Points
			double CostWithoutCoreset{0};
			// Total cost when using KMeans on a random Coreset
			double RonT {0};
			double TonR {0};
			double RonR {0};
			// Total cost when using an ILP on all points. Each center is a point from the input
			double ILPCost{0};
			// Total cost when using an MIQP on all points. Calcuates new centers
			double MIQPCost{0};
			// Average runtime of running KMeans on the total Dataset without using a coreset
			std::chrono::microseconds fullTime{0};

			// For the last statistics calculation, also calculate the cost without using a coreset,
			// The cost when using a random coreset, the cost when using an ILP
			std::shared_ptr<clusteringStatic::ClusteringResult> fullResult;
			auto startFull = std::chrono::high_resolution_clock::now();
			auto fullKMeansPP = staticKMeans.pickInitialCentersViaKMeansPlusPlus(staticPoints, staticWeights);
			auto fullKMeansPPCenters = staticKMeans.copyRows(staticPoints, *fullKMeansPP->getClusterIndices());
			fullResult = staticKMeans.runLloydsAlgorithm(staticPoints, fullKMeansPPCenters, staticWeights);
			auto stopFull = std::chrono::high_resolution_clock::now();

			fullTime = duration_cast<microseconds>(stopFull - startFull);

			CostWithoutCoreset = calcTotalSquaredCost(staticPoints, fullResult->getCentroids(), staticWeights) / blaze::sum(staticWeights);
			CostWithoutCoresetWithoutLoyd = calcTotalSquaredCost(staticPoints, fullKMeansPPCenters, staticWeights) / blaze::sum(staticWeights);

			
			/*
			double t = (dynamicKMeansTime+dynamicCoresetTime).count()/static_cast<double>(skip);
			double PTime = (t/1000000.)*PTimeFactor;
			ILPCost = getBaseLineCostWithCentersAsPoints(staticPoints,staticWeights,nrClusters, PTime, ILPTime) / blaze::sum(staticWeights);
			MIQPCost = getBaseLineCostWithArbitraryCenters(staticPoints,staticWeights,nrClusters, fullKMeansPP->getClusterIndices(), PTime, MIQPTime) / blaze::sum(staticWeights);
			*/
			



			randomTime = std::chrono::microseconds::zero();

			blaze::DynamicMatrix<double> randomCoreset(nrPointsFinalCoreset, dimensions);
			blaze::DynamicVector<double> randomCoresetWeights(nrPointsFinalCoreset);

			const int randomRepetitions{1};
			std::mt19937 randomEngine;
			std::array<int, 624> seedData;
			std::random_device randomDevice;
			std::generate_n(seedData.data(), seedData.size(), std::ref(randomDevice));
			std::seed_seq randomSeq(std::begin(seedData), std::end(seedData));
			randomEngine.seed(randomSeq);

			// Create a random Coresest, run KMeans with it and calculate the cost.
			// Do this "randomReptions" times and average

			for (int r = 0; r < randomRepetitions; r++)
			{
				auto startRand = std::chrono::high_resolution_clock::now();

				//std::discrete_distribution<size_t> weightedChoice(staticWeights.begin(), staticWeights.end());
				for (int e = 0; e < nrPointsFinalCoreset; e++)
				{
					//int randRow = weightedChoice(randomEngine);
					int randRow = rand() % staticPoints.rows();

					assert(randRow < staticPoints.rows());
					for (int c = 0; c < dimensions; c++)
						randomCoreset.at(e, c) = staticPoints.at(randRow, c);
					randomCoresetWeights.at(e) = staticWeights.at(randRow);
				}
				auto stopRand = std::chrono::high_resolution_clock::now();
				auto randomCoresetResult = staticKMeans.run(randomCoreset, randomCoresetWeights, false);
				randomTime += duration_cast<microseconds>(stopRand - startRand);


				RonT += calcTotalSquaredCost(staticPoints, randomCoresetResult->getCentroids(), staticWeights) / blaze::sum(staticWeights);
				RonR += calcTotalSquaredCost(randomCoreset, randomCoresetResult->getCentroids(), randomCoresetWeights) / blaze::sum(randomCoresetWeights);
				TonR += calcTotalSquaredCost(randomCoreset, fullResult->getCentroids(), randomCoresetWeights) / blaze::sum(randomCoresetWeights);
			}

			randomTime /= randomRepetitions;
			RonT /= randomRepetitions;
			RonR /= randomRepetitions;
			TonR /= randomRepetitions;

			double TonT {calcTotalSquaredCost(staticPoints, fullResult->getCentroids(), staticWeights) / blaze::sum(staticWeights)};

			double DonD {calcTotalSquaredCost(coresetPoints, dynamicResult->getCentroids(), coresetWeights) / blaze::sum(coresetWeights)};
			double DonT {calcTotalSquaredCost(staticPoints, dynamicResult->getCentroids(), staticWeights) / blaze::sum(staticWeights)};
			double TonD {calcTotalSquaredCost(coresetPoints, fullResult->getCentroids(), coresetWeights) / blaze::sum(coresetWeights)};

			double OonO {calcTotalSquaredCost(coresetPointsOpt, dynamicResultOpt->getCentroids(), coresetWeightsOpt) / blaze::sum(coresetWeightsOpt)};
			double OonT {calcTotalSquaredCost(staticPoints, dynamicResultOpt->getCentroids(), staticWeights) / blaze::sum(staticWeights)};
			double TonO {calcTotalSquaredCost(coresetPointsOpt, fullResult->getCentroids(), coresetWeightsOpt) / blaze::sum(coresetWeightsOpt)};

			double SonS {calcTotalSquaredCost(coresetPointsStat, staticResult->getCentroids(), coresetWeightsStat) / blaze::sum(coresetWeightsStat)};
			double SonT {calcTotalSquaredCost(staticPoints, staticResult->getCentroids(), staticWeights) / blaze::sum(staticWeights)};
			double TonS {calcTotalSquaredCost(coresetPointsStat, fullResult->getCentroids(), coresetWeightsStat) / blaze::sum(coresetWeightsStat)};






			// Write the statistics to the output file and the current times to the console
			std::cout << "Average Dynamic Time: " << (dynamicCoresetTime + dynamicKMeansTime).count() / skip << "ms\n";
			std::cout << "Average Static Time:  " << (staticCoresetTime + staticKMeansTime).count() << "ms\n";
			outFile << i << " ";

			outFile << dynamicCoresetTime.count() / skip << " ";
			outFile << dynamicKMeansTime.count() / skip << " ";
			outFile << DonT << " "; //Dyn cost
			outFile << std::max(std::max(TonT/TonD,TonD/TonT), std::max(DonD/DonT,DonT/DonD)) << " "; //Dyn coreset quality
			outFile << TonT/DonT << " "; //Dyn k-means quality

			outFile << dynamicCoresetTimeOpt.count() / skip << " ";
			outFile << dynamicKMeansTimeOpt.count() / skip << " ";
			outFile << OonT << " "; //Dyn Opt cost
			outFile << std::max(std::max(TonT/TonO,TonO/TonT), std::max(OonO/OonT,OonT/OonO)) << " "; //Dyn Opt coreset quality
			outFile << TonT/OonT << " "; //Dyn Opt k-means quality

			outFile << staticCoresetTime.count() << " ";
			outFile << staticKMeansTime.count() << " ";
			outFile << SonT << " "; //Stat cost
			outFile << std::max(std::max(TonT/TonS,TonS/TonT), std::max(SonS/SonT,SonT/SonS)) << " "; //Stat coreset quality
			outFile << TonT/SonT << " "; //Stat k-means quality

			outFile << randomTime.count() << " ";
			outFile << RonT << " "; //Rand cost
			outFile << std::max(std::max(TonT/TonR,TonR/TonT), std::max(RonR/RonT,RonT/RonR)) << " "; //Dyn coreset quality
			outFile << TonT/RonT << " "; //Rand k-means quality

			outFile << fullTime.count() << " ";
			outFile << CostWithoutCoreset << " ";

			outFile << ILPTime*1000000 << " ";
			outFile << ILPCost << " ";
			outFile << MIQPTime*1000000 << " ";
			outFile << MIQPCost << std::endl;


			// Reset the time measurements
		}

		if (i % skip == 0)
		{
			staticCoresetTime = std::chrono::microseconds::zero();
			dynamicCoresetTime = std::chrono::microseconds::zero();
			dynamicCoresetTimeOpt = std::chrono::microseconds::zero();
			staticKMeansTime = std::chrono::microseconds::zero();
			dynamicKMeansTime = std::chrono::microseconds::zero();
			dynamicKMeansTimeOpt = std::chrono::microseconds::zero();
		}

	}
	outFile.close();

	return 0;
}

int cardinalityOfCoreset(double epsilon, double lambda, int k, int d)
{
	if (::useFinalSize == true)
		return ::sizeOfFinalCoreset;
	if (epsilon <= 0 || lambda > 0.5)
		return -1;

	int c = 5;

	return (c / (epsilon * epsilon) * (k * d + blaze::log(1. / lambda)));
}
