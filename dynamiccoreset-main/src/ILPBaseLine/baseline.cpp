#include "baseline.h"
#include <memory>
#include <string>
#include <vector>

/**
 * @brief Calculate the squared distance between two points
 * @param points matrix containing at least the two points for which the distance should be calculated. Each row should contain the
 * coordinates of a point
 * @param index1 Row in the matrix "points" that contains the coordinates of the first point
 * @param index2 Row in the matrix "points" that contains the coordinates of the second point
 * @return Squared distance between the two points
 */
double calcSquDist(const blaze::DynamicMatrix<double> &points, int index1, int index2);
double getBaseLineCostWithCentersAsPoints(const blaze::DynamicMatrix<double> &points, const blaze::DynamicVector<double> &weights, const int nrCenters, const double timeCutOff, double& runtime)
{
	/**
	 * Uses the following Program:
	 *
	 * Minimize 	Sum (d_ij * x_ij)
	 *
	 * s.t.		1.	Sum_i (x_ii) 	<= k
	 * 			2.	Sum_j (x_ij)	>= 1 			for all i
	 * 			3.	Sum_i!=j (x_ij)	<= x_jj * n 	for all j
	 *		 	4.	x_ij 			in {0,1}
	 *
	 * d_ij: Distance between points i and j. Here stored in the matrix "distances" in a nxn matrix where n is the number of points
	 * x_ij: Is 1 iff point i is part of the cluster whose center is point j and 0 otherwise. Stored in a nxn matrix "isPartOfCluster"
	 * 		 When i == j, x_ij is 1 if the point i is a center of a cluster and 0 if the point is not a center of a cluster
	 *
	 * The variable "obj" stores the objective function.
	 * The constraint 1. is stored in the variable "constraint1"
	 * The constraint 2. is stored in the variable "constraint2"
	 * The constraint 3. is stored in the variable "constraint3"
	 * The constraint 4. is implicitly stored by setting the elements of "isPartOfCluster" to binary
	 *
	 */
	try
	{
		GRBEnv env = GRBEnv(true);
		env.set("TimeLimit",std::to_string(timeCutOff));
		env.set("LogFile", "ILPLog.log");
		env.start();

		GRBModel model = GRBModel(env);

		GRBVar *isPartOfCluster = new GRBVar[points.rows() * points.rows()];
		double distances[points.rows() * points.rows()];
		int index{0};
		GRBLinExpr obj = 0.0;
		GRBLinExpr constraint1 = 0.0;
		GRBLinExpr constraint2[points.rows()];
		GRBLinExpr constraint3[points.rows()];

		for (int row = 0; row < points.rows(); row++)
		{
			constraint2[row] = 0;
			for (int col = 0; col < points.rows(); col++)
			{
				isPartOfCluster[index] = model.addVar(0, 1, 0, GRB_BINARY);
				distances[index] = calcSquDist(points, row, col);

				if (col == row)
					constraint1 += isPartOfCluster[index];
				else
				{
					obj += (distances[index] * weights.at(row)) * isPartOfCluster[index];
					constraint3[col] += isPartOfCluster[index];
				}
				constraint2[row] += isPartOfCluster[index];

				index++;
			}
			model.addConstr(constraint2[row] >= 1);
		}
		for (int col = 0; col < points.rows(); col++)
		{
			model.addConstr(constraint3[col] <= isPartOfCluster[col * points.rows() + col] * points.rows());
		}
		model.addConstr(constraint1 <= nrCenters);
		assert(index == (points.rows() * points.rows()));
		model.setObjective(obj, GRB_MINIMIZE);

		model.optimize();


		//Print the resulting matrix "isPartOfCuster"
		/*
		auto vars = model.getVars();
		index = 0;
		for (int row = 0; row < points.rows(); row++)
		{
			for (int col = 0; col < points.rows(); col++)
			{
				std::cout << static_cast<unsigned int>(vars[index].get(GRB_DoubleAttr_X)) << " ";
				index++;
			}
			std::cout << "\n";
		}
		delete [] vars;
		*/

		delete [] isPartOfCluster;

		runtime = model.get(GRB_DoubleAttr_Runtime);
		return model.get(GRB_DoubleAttr_ObjVal);




	} catch ( GRBException e ) {
		std::cout << " Error code = " << e . getErrorCode () << std::endl ;
		std::cout << e . getMessage () << std::endl ;
	} catch (...) {
		std::cout << " Exception during optimization " << std::endl ;
	}

	return -1;

}

double getBaseLineCostWithArbitraryCenters (blaze::DynamicMatrix<double>& points, const blaze::DynamicVector<double>& weights, const int nrCenters, const std::shared_ptr<std::vector<size_t>> centerHints, const double timeCutOff, double& runtime)
{
	/**
	 * Uses the following Program:
	 *
	 * Minimize 	Sum_ij (d_ij)
	 *
	 * s.t.		1.	d_ij 			>= Sum_r (p_ir - c_ir)^2 - M * (1 - x_ij)	for all i and for all j
	 * 			2.	Sum_j (x_ij)	>= 1 										for all i
	 * 			3.	d_ij 			>= 0										for all i and for all j
	 *		 	4.	x_ij 			in {0,1}
	 *
	 * p_ir: Coordinate of the r'th dimension of point i
	 * c_ir: Coordinate of the r'th dimension of center i. Is a variable that will be determined and stored in "newCenters"
	 * M: A large number
	 * d_ij: Distance between the Point i and the center j if i is part of the cluster with center j. Else, 0
	 * x_ij: Is 1 iff point i is part of the cluster whose center is center j and 0 otherwise. Stored in a nxn matrix "isPartOfCluster"
	 * 		 When i == j, x_ij is 1 if the point i is a center of a cluster and 0 if the point is not a center of a cluster
	 *
	 * The variable "obj" stores the objective function.
	 * The constraint 1. is stored in the variable "constraint1". (a-b)^2 = a^2 - 2ab - b^2 was used
	 * The constraint 2. is stored in the variable "constraint2"
	 * The constraint 2. is set at line 179 
	 * The constraint 4. is implicitly stored by setting the elements of "isPartOfCluster" to binary
	 *
	 */
	try{
		GRBEnv env = GRBEnv(true);
		env.set("TimeLimit",std::to_string(timeCutOff));
		env.set("LogFile","ILPLog.log");
		env.start();

		int M = 100000;

		GRBModel model = GRBModel(env);
		GRBVar distances [points.rows()*nrCenters];
		GRBVar newCenters [points.columns()*nrCenters];
		GRBVar isPartOfCluster [points.rows()*nrCenters];
		GRBLinExpr obj = 0.0;
		GRBQuadExpr constraint1 [points.rows()*nrCenters];
		GRBLinExpr constraint2 [points.rows()];


		for (int p = 0; p < points.rows(); p++)
		{
			for (int c = 0; c < nrCenters; c++)
			{
				distances[p*nrCenters+c] = model.addVar(-GRB_INFINITY,GRB_INFINITY,0,GRB_CONTINUOUS,"d_"+std::to_string(p)+"_"+std::to_string(c));
				for (int d = 0; d < points.columns(); d++)
				{
					if (p == 0)
					{
						newCenters[c*points.columns()+d] = model.addVar(-GRB_INFINITY,GRB_INFINITY,0,GRB_CONTINUOUS,"c_"+std::to_string(c)+"_"+std::to_string(d));
						// Uncomment this to give the Model the centers of KMeans++ as hints. Acutally makes it worse whenever I tried it.
						//newCenters[c*points.columns()+d].set(GRB_DoubleAttr_VarHintVal,points.at(centerHints->at(c),d));
					}
					constraint1[p*nrCenters+c] += (points.at(p,d)*points.at(p,d));
					constraint1[p*nrCenters+c] -= (2*points.at(p,d)*newCenters[c*points.columns()+d]);
					constraint1[p*nrCenters+c] += (newCenters[c*points.columns()+d]*newCenters[c*points.columns()+d]);
				}
				isPartOfCluster[p*nrCenters+c] = model.addVar(0,1,0,GRB_BINARY,"x_"+std::to_string(p)+"_"+std::to_string(c));

				constraint1[p*nrCenters+c] -= (M * (1-isPartOfCluster[p*nrCenters+c]));
				constraint2[p] += isPartOfCluster[p*nrCenters+c];

				model.addQConstr(distances[p*nrCenters+c] - constraint1[p*nrCenters+c] >= 0);
				model.addConstr(distances[p*nrCenters+c] >= 0);
				obj += distances[p*nrCenters+c];
			}
			model.addConstr(constraint2[p] >= 1);
		}
		model.setObjective(obj,GRB_MINIMIZE);
		model.optimize();

		runtime = model.get(GRB_DoubleAttr_Runtime);
		double x = model.get(GRB_DoubleAttr_ObjVal);

		//Print the variables x_ij, then d_ij and finally c_ij
		/*
		for (int p = 0; p < points.rows(); p++)
		{
			for (int c = 0; c < nrCenters; c++)
			{
				std::cout <<"x_"+std::to_string(p)+"_"+std::to_string(c) << ": " << model.getVarByName("x_"+std::to_string(p)+"_"+std::to_string(c)).get(GRB_DoubleAttr_X) << " ";

			}
			std::cout << "\n";
		}
		std::cout << "\n";

		for (int p = 0; p < points.rows(); p++)
		{
			for (int c = 0; c < nrCenters; c++)
			{
				std::cout << "d_"+std::to_string(p)+"_"+std::to_string(c) << ": " << model.getVarByName("d_"+std::to_string(p)+"_"+std::to_string(c)).get(GRB_DoubleAttr_X) << " ";

			}
			std::cout << "\n";
		}
		std::cout << "\n";

		for (int c = 0; c < nrCenters; c++)
		{
			for (int d = 0; d < points.columns(); d++)
			{
				std::cout << "c_"+std::to_string(c)+"_"+std::to_string(d) << ": " << model.getVarByName("c_"+std::to_string(c)+"_"+std::to_string(d)).get(GRB_DoubleAttr_X) << " ";

			}
			std::cout << "\n";
		}
		std::cout << "\n";
		*/

		return x;
	}
	catch (GRBException e)
	{
		std::cout << " Error code = " << e.getErrorCode() << std::endl;
		std::cout << e.getMessage() << std::endl;
	}
	catch (...)
	{
		std::cout << " Exception during optimization " << std::endl;
	}

	return 0;
}

double calcSquDist(const blaze::DynamicMatrix<double> &points, int index1, int index2)
{
	if (index1 == index2)
		return 0;
	return blaze::sqrNorm(blaze::row(points, index2) - blaze::row(points, index1));
}
