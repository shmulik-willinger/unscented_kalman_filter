#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  // radar measurement dimension
  n_zrad_ = 3;

  // radar measurement dimension
  n_zlas_ = 2;

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  int num_sigma_points = n_aug_ * 2 + 1;

  // Sigma point spreading parameter
  lambda_ = 0;

  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, num_sigma_points);

  // time when the state is true, in us
  time_us_ = 0;

  // the current NIS for laser and radar
  NIS_radar_ = 0;
  NIS_laser_ = 0;

  // Weights of sigma points
  weights_ = VectorXd::Zero(num_sigma_points);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i<2 * n_aug_ + 1; i++)
	  weights_(i) = 0.5 / (n_aug_ + lambda_);


  // Noise matrices
  R_radar = MatrixXd(3, 3);
  R_laser = MatrixXd(2, 2);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) 
{
	if (!is_initialized_)
	{

		if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
		{
			cout << "Initialize with first measurement- RADAR" << endl;
			double rho = meas_package.raw_measurements_(0);
			double phi = meas_package.raw_measurements_(1);
			double rhodot = meas_package.raw_measurements_(2);

			double x = rho * cos(phi);
			double y = rho * sin(phi);
			double vx = rhodot * cos(phi);
			double vy = rhodot * sin(phi);

			// polar to cartesian - r * cos(angle) for x and r * sin(angle) for y
			x_ << x, y, vx, vy;

			// create state covariance matrix
			P_ << std_radr_*std_radr_, 0, 0, 0, 0,
				0, std_radr_*std_radr_, 0, 0, 0,
				0, 0, 1, 0, 0,
				0, 0, 0, std_radphi_, 0,
				0, 0, 0, 0, std_radphi_;

			// create R for update noise later
			R_radar << std_radr_*std_radr_, 0, 0,
				0, std_radphi_*std_radphi_, 0,
				0, 0, std_radrd_*std_radrd_;
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
		{
			cout << "Initialize with first measurement- LASER" << endl;
			x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 4, 0.5, 0.0;

			// create state covariance matrix
			P_ << std_laspx_*std_laspx_, 0, 0, 0, 0,
				0, std_laspy_*std_laspy_, 0, 0, 0,
				0, 0, 1, 0, 0,
				0, 0, 0, 1, 0,
				0, 0, 0, 0, 1;

			// create R for update noise later
			R_laser << std_laspx_*std_laspx_, 0,
				0, std_laspy_*std_laspy_;
		}

		// save the first timestamp
		time_us_ = meas_package.timestamp_;

		is_initialized_ = true;
		return;
	}

	if (meas_package.timestamp_ != time_us_)
	{
		// computing the time in seconds between the current and previous measurements
		float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
		time_us_ = meas_package.timestamp_;

		Prediction(dt);
	}

	// Measurement updates
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
	{
		int n_zrad_ = 2;
		//mean predicted measurement
		VectorXd z_pred = VectorXd::Zero(n_zrad_);
		//measurement covariance matrix S
		MatrixXd S = MatrixXd::Zero(n_zrad_, n_zrad_);
		// cross-correlation matrix Tc
		MatrixXd Tc = MatrixXd::Zero(n_x_, n_zrad_);
		// get predictions for x,S and Tc in RADAR space
		PredictRadarMeasurement(z_pred, S, Tc);
		UpdateRadar(meas_package, z_pred, Tc, S);
	}
	else
	{
		int n_zlas_ = 3;
		//mean predicted measurement
		VectorXd z_pred = VectorXd::Zero(n_zlas_);
		//measurement covariance matrix S
		MatrixXd S = MatrixXd::Zero(n_zlas_, n_zlas_);
		// cross-correlation matrix Tc
		MatrixXd Tc = MatrixXd::Zero(n_x_, n_zlas_);
		// get predictions for x,S and Tc in Lidar space
		PredictLidarMeasurement(z_pred, S, Tc);
		UpdateLidar(meas_package, z_pred, Tc, S);
	}
}

void UKF::PredictRadarMeasurement(VectorXd &z_out, MatrixXd &S_out, MatrixXd &Tc_out) {

	//radar can measure r, phi, and r_dot

	//create matrix for sigma points in measurement space
	static MatrixXd Zsig = MatrixXd(n_zrad_, 2 * n_aug_ + 1);

	//transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points

												// extract values for better readibility
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double yaw = Xsig_pred_(3, i);

		// Avoid division by zero
		if (fabs(p_x) <= 0.0001) {
			p_x = 0.0001;
		}
		if (fabs(p_y) <= 0.0001) {
			p_y = 0.0001;
		}


		double v1 = cos(yaw)*v;
		double v2 = sin(yaw)*v;

		// measurement model
		Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);                        //r
		Zsig(1, i) = atan2(p_y, p_x);                                 //phi
		Zsig(2, i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
	}

	//mean predicted measurement
	static VectorXd z_pred = VectorXd(n_zrad_);
	z_pred.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}

	//measurement covariance matrix S
	static MatrixXd S = MatrixXd(n_zrad_, n_zrad_);
	S.fill(0.0);
	//create matrix for cross correlation Tc
	static MatrixXd Tc = MatrixXd(n_x_, n_zrad_);
	Tc.fill(0.0);

	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
												//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;

		//angle normalization
		while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

		S = S + weights_(i) * z_diff * z_diff.transpose();

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		//angle normalization
		while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();

	}

	//add measurement noise covariance matrix
	static MatrixXd R = MatrixXd(n_zrad_, n_zrad_);
	R << pow(std_radr_, 2), 0, 0,
		0, pow(std_radphi_, 2), 0,
		0, 0, pow(std_radrd_, 2);
	S = S + R;

	//write result
	z_out = z_pred;
	S_out = S;
	Tc_out = Tc;

	return;

}


void UKF::PredictLidarMeasurement(VectorXd &z_out, MatrixXd &S_out, MatrixXd &Tc_out) {

	//create matrix for sigma points in measurement space
	static MatrixXd Zsig = MatrixXd(n_zlas_, 2 * n_aug_ + 1);

	//transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points

												// measurement model
		Zsig(0, i) = Xsig_pred_(0, i);          //px
		Zsig(1, i) = Xsig_pred_(1, i);          //py

	}

	//mean predicted measurement
	static VectorXd z_pred = VectorXd(n_zlas_);
	z_pred.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}

	//measurement covariance matrix S
	static MatrixXd S = MatrixXd(n_zlas_, n_zlas_);
	S.fill(0.0);
	//create matrix for cross correlation Tc
	static MatrixXd Tc = MatrixXd(n_x_, n_zlas_);
	Tc.fill(0.0);

	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
												//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;

		S = S + weights_(i) * z_diff * z_diff.transpose();

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;

		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();

	}

	//add measurement noise covariance matrix
	static MatrixXd R = MatrixXd(n_zlas_, n_zlas_);
	R << pow(std_laspx_, 2), 0,
		0, pow(std_laspy_, 2);
	S = S + R;

	//write result
	z_out = z_pred;
	S_out = S;
	Tc_out = Tc;

	return;

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) 
{
	/////////////////////// Generate Sigma points ///////////////////////

	//define spreading parameter
	lambda_ = 3 - n_x_;

	//create sigma point matrix
	MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);

	//calculate square root of P
	MatrixXd A = P_.llt().matrixL();

	//set first column of sigma point matrix
	Xsig.col(0) = x_;

	//set remaining sigma points
	for (int i = 0; i < n_x_; i++)
	{
		Xsig.col(i + 1) = x_ + sqrt(lambda_ + n_x_) * A.col(i);
		Xsig.col(i + 1 + n_x_) = x_ - sqrt(lambda_ + n_x_) * A.col(i);
	}

	/////////////////////// Augmented Sigma points ///////////////////////

	//define spreading parameter
	double lambda = 3 - n_aug_;

	//create augmented mean vector
	VectorXd x_aug = VectorXd(7);

	//create augmented state covariance
	MatrixXd P_aug = MatrixXd(7, 7);

	//create sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

	//create augmented mean state
	x_aug.head(5) = x_;
	x_aug(5) = 0;
	x_aug(6) = 0;

	//create augmented covariance matrix
	P_aug.fill(0.0);
	P_aug.topLeftCorner(5, 5) = P_;
	MatrixXd Q = MatrixXd(2, 2);
	Q << std_a_*std_a_, 0,
		0, std_yawdd_*std_yawdd_;
	P_aug.bottomRightCorner(2, 2) = Q;
 
	//create square root matrix
	MatrixXd L = P_aug.llt().matrixL();

	//create augmented sigma points
	Xsig_aug.col(0) = x_aug;
	for (int i = 0; i< n_aug_; i++)
	{
		Xsig_aug.col(i + 1) = x_aug + sqrt(lambda + n_aug_) * L.col(i);
		Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda + n_aug_) * L.col(i);
	}


	/////////////////////// Sigma Point Prediction ///////////////////////

	//create matrix with predicted sigma points as columns
	MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

	//predict sigma points
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		//extract values for better readability
		double p_x = Xsig_aug(0, i);
		double p_y = Xsig_aug(1, i);
		double v = Xsig_aug(2, i);
		double yaw = Xsig_aug(3, i);
		double yawd = Xsig_aug(4, i);
		double nu_a = Xsig_aug(5, i);
		double nu_yawdd = Xsig_aug(6, i);

		//predicted state values
		double px_p, py_p;

		//avoid division by zero
		if (fabs(yawd) > 0.001) {
			px_p = p_x + v / yawd * (sin(yaw + yawd*delta_t) - sin(yaw));
			py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd*delta_t));
		}
		else {
			px_p = p_x + v*delta_t*cos(yaw);
			py_p = p_y + v*delta_t*sin(yaw);
		}

		double v_p = v;
		double yaw_p = yaw + yawd*delta_t;
		double yawd_p = yawd;

		//add noise
		px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
		py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
		v_p = v_p + nu_a*delta_t;

		yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
		yawd_p = yawd_p + nu_yawdd*delta_t;

		//write predicted sigma point into right column
		Xsig_pred(0, i) = px_p;
		Xsig_pred(1, i) = py_p;
		Xsig_pred(2, i) = v_p;
		Xsig_pred(3, i) = yaw_p;
		Xsig_pred(4, i) = yawd_p;
	}

	static VectorXd x_pred = VectorXd(n_x_);
	static MatrixXd P_pred = MatrixXd(n_x_, n_x_);
	PredictMeanAndCovariance(x_pred, P_pred);

	x_ = x_pred;
	P_ = P_pred;
	Xsig_pred_ = Xsig_pred;

	return;
}

void UKF::PredictMeanAndCovariance(VectorXd& x_out, MatrixXd& P_out)
{
	//create vector for predicted state
	VectorXd x = VectorXd(n_x_);

	//create covariance matrix for prediction
	MatrixXd P = MatrixXd(n_x_, n_x_);

	//predicted state mean
	x.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
		x = x + weights_(i) * Xsig_pred_.col(i);
	}

	//predicted state covariance matrix
	P.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x;
		//angle normalization
		while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

		P = P + weights_(i) * x_diff * x_diff.transpose();
	}
	x_out = x;
	P_out = P;

	return;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package, VectorXd &z_pred, MatrixXd &Tc, MatrixXd &S) 
{
	int n_zlas_ = 2;
	//mean predicted measurement
	static VectorXd z = VectorXd::Zero(2);
	z << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1);

	//Kalman gain K;
	static MatrixXd K = MatrixXd::Zero(n_x_, n_zlas_);
	K = Tc * S.inverse();

	//residual
	static VectorXd z_diff = VectorXd::Zero(n_zlas_);
	z_diff = z - z_pred;

	//update state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K*S*K.transpose();
	NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;

	return;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package, VectorXd &z_pred, MatrixXd &Tc, MatrixXd &S) 
{
	int n_zrad_ = 3;
	//mean predicted measurement  rho, phi, rho_dot
	static VectorXd z = VectorXd::Zero(n_zrad_);
	z << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), meas_package.raw_measurements_(2);

	//Kalman gain K;
	static MatrixXd K = MatrixXd::Zero(n_x_, n_zrad_);
	K = Tc * S.inverse();

	//residual
	static VectorXd z_diff = VectorXd::Zero(n_zrad_);
	z_diff = z - z_pred;

	//angle normalization
	while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
	while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

	//update state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K*S*K.transpose();
	NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

	return;
}
