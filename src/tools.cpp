#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth) 
{
	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;

	// Validating the inputs
	if (ground_truth.size() == 0 || estimations.size() != ground_truth.size())
	{
		cout << "Invalid inputs size. estimations= " << estimations.size() << ", ground_truth= " << ground_truth.size() << endl;
		return rmse;
	}

	//accumulate squared residuals
	for (unsigned int i = 0; i < estimations.size(); i++)
	{
		VectorXd residual = estimations[i] - ground_truth[i];
		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	//calculate the mean and squared root
	rmse = rmse / estimations.size();
	rmse = rmse.array().sqrt();

	// print the RMSE
	//cout << "RMSE = " << rmse << endl;
	return rmse;
}