#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rsme(4);
  rsme << 0,0,0,0;

  if ( estimations.size() == 0 || estimations.size() != ground_truth.size() ) {
		std::cout << "estimations or ground_truth vector size problem either 0 or not the same size" << std::endl;
		return rsme;
	}

	for (unsigned int i = 0; i < estimations.size(); i++) {
		VectorXd diff = estimations[i] - ground_truth[i];
		diff = diff.array()*diff.array();
		rsme += diff; 
	}

	rsme = rsme / estimations.size();

	rsme = rsme.array().sqrt();

  //DEBUGGING PURPOSE ONLY
	// cout << "----------------" << endl;
	// cout << "rmse_x:" << rsme(0) << endl;
	// cout << "rmse_y:" << rsme(1) << endl;
	// cout << "rmse_vx:" << rsme(2) << endl;
	// cout << "rmse_vy:" << rsme(3) << endl; 

  return rsme;
}