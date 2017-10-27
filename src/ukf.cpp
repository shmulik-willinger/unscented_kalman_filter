#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
* Initializes Unscented Kalman filter
*/
UKF::UKF() 
{
	// Set initialization to false initially
	is_initialized_ = false;

	// State dimension
	n_x_ = 5;

	// Augmented state dimension
	n_aug_ = 7;

	// radar measurement dimension
	n_radar_ = 3;

	// radar measurement dimension
	n_laser_ = 2;

	// sigma points dimention
	n_sig_points_ = 2 * n_aug_ + 1;

	// Sigma point spreading parameter
	lambda_ = 3 - n_aug_;

	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;

	// initial state vector
	x_ = VectorXd(5);

	// initial covariance matrix
	P_ = MatrixXd::Zero(5, 5);

	// initialize the sigma points matrix with zeroes
	Xsig_pred_ = MatrixXd::Zero(n_x_, n_sig_points_);

	///* time when the state is true, in us
	time_us_ = 0;

	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 3;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 0.7;

	// Laser measurement noise standard deviation position1 in m
	std_laspx_ = 0.15;

	// Laser measurement noise standard deviation position2 in m
	std_laspy_ = 0.15;

	// Radar measurement noise standard deviation radius in m
	std_radr_ = 0.3;

	// Radar measurement noise standard deviation angle in rad
	std_radphi_ = 0.03;

	// Radar measurement noise standard deviation radius change in m/s
	std_radrd_ = 0.35;

	// Weights of sigma points
	weights_ = VectorXd::Zero(n_sig_points_);

	// the current NIS for radar
	NIS_radar_ = 0.0;

	// the current NIS for laser
	NIS_laser_ = 0.0;

	// Lidar noise covariance matrix
	R_lidar = MatrixXd(n_laser_, n_laser_);
	R_lidar << std_laspx_*std_laspx_, 0,
		0, std_laspy_*std_laspy_;
	
	// Radar noise covariance matrix
	R_radar = MatrixXd(n_radar_, n_radar_);
	R_radar << std_radr_*std_radr_, 0, 0,
		0, std_radphi_*std_radphi_, 0,
		0, 0, std_radrd_*std_radrd_;

	// set weights
	weights_(0) = lambda_ / (lambda_ + n_aug_);
	for (int i = 1; i< n_sig_points_; i++)
		weights_(i) = 0.5 / (n_aug_ + lambda_);

	return;

}

UKF::~UKF() {}

/**
* @param {MeasurementPackage} meas_package The latest measurement data of
* either radar or laser.
*/
void UKF::ProcessMeasurement(const MeasurementPackage meas_package) 
{
	if (!is_initialized_) 
	{
		// initialize the state covariance matrix
		P_ = MatrixXd::Identity(5, 5);

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
			x_ << x, y, sqrt(vx*vx + vy*vy), 0, 0;
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
		{
			cout << "Initialize with first measurement- LASER" << endl;
			x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 0.0, 0.0, 0.0;
		}

		// save the first timestamp
		time_us_ = meas_package.timestamp_;
		is_initialized_ = true;

		return;
	}

	// computing the time in seconds between the current and previous measurements
	float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
	Prediction(dt);

	/*****************************************************************************
	* Update the state and covariance matrices by the sensor type
	****************************************************************************/
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) 
	{
		// mean predicted measurement
		VectorXd z_pred = VectorXd::Zero(n_radar_);
		// cross-correlation matrix
		MatrixXd Tc = MatrixXd::Zero(n_x_, n_radar_);
		// measurement covariance matrix
		MatrixXd S = MatrixXd::Zero(n_radar_, n_radar_);

		// get predictions 
		PredictRadarMeasurement(z_pred, S, Tc);
		// update the state using the RADAR measurement
		UpdateRadar(meas_package, z_pred, Tc, S);
		time_us_ = meas_package.timestamp_;

	}
	else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
	{
		// mean predicted measurement
		VectorXd z_pred = VectorXd::Zero(n_laser_);
		// cross-correlation matrix
		MatrixXd Tc = MatrixXd::Zero(n_x_, n_laser_);
		//measurement covariance matrix
		MatrixXd S = MatrixXd::Zero(n_laser_, n_laser_);

		PredictLidarMeasurement(z_pred, S, Tc);
		// update the state using the LIDAR measurement
		UpdateLidar(meas_package, z_pred, Tc, S);
		time_us_ = meas_package.timestamp_;
	}

	return;
}

void UKF::AugmentedSigmaPoints(MatrixXd &Xsig_out) 
{
	//create augmented mean vector
	static VectorXd x_aug = VectorXd(n_aug_);

	//create augmented state covariance
	static MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

	//create sigma point matrix
	static MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_points_);

	//create augmented mean state
	x_aug.head(5) = x_;
	x_aug(5) = 0;
	x_aug(6) = 0;

	//create augmented covariance matrix
	P_aug.fill(0.0);
	P_aug.topLeftCorner(5,5) = P_;
	P_aug(5,5) = std_a_ * std_a_;
	P_aug(6,6) = std_yawdd_ * std_yawdd_;

	//create square root matrix
	MatrixXd L = P_aug.llt().matrixL();

	//create augmented sigma points
	Xsig_aug.col(0) = x_aug;
	double sqr = sqrt(lambda_ + n_aug_);
	for (int i = 0; i< n_aug_; i++)
	{
		Xsig_aug.col(i + 1) = x_aug + sqr * L.col(i);
		Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqr * L.col(i);
	}

	Xsig_out = Xsig_aug;

	return;

}

/**
* @param MatrixXd &Xsig_out. predicts augmented sigma points
*/
void UKF::SigmaPointPrediction(const MatrixXd &Xsig_aug, const double delta_t, MatrixXd  &Xsig_out) {

	//create matrix with predicted sigma points as columns
	static MatrixXd Xsig_pred = MatrixXd(n_x_, n_sig_points_);

	//predict sigma points
	for (int i = 0; i< n_sig_points_; i++)
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

	Xsig_out = Xsig_pred;
	return;
}

void UKF::PredictMeanAndCovariance(const MatrixXd &Xsig_pred) 
{
	static VectorXd x_pred = VectorXd(n_x_);
	static MatrixXd P_pred = MatrixXd(n_x_, n_x_);

	//create vector for predicted state
	static VectorXd x = VectorXd(n_x_);

	//create covariance matrix for prediction
	static MatrixXd P = MatrixXd(n_x_, n_x_);

	//predicted state mean
	x.fill(0.0);
	//iterate over sigma points
	for (int i = 0; i < n_sig_points_; i++) 
		x = x + weights_(i) * Xsig_pred.col(i);

	//predicted state covariance matrix
	P.fill(0.0);
	//iterate over sigma points
	for (int i = 0; i < n_sig_points_; i++) 
	{  
		VectorXd x_diff = Xsig_pred.col(i) - x;
		//angle normalization
		while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

		P = P + weights_(i) * x_diff * x_diff.transpose();
	}

	//set results
	x_ = x;
	P_ = P;
	Xsig_pred_ = Xsig_pred;
	return;
}

void UKF::PredictRadarMeasurement(VectorXd &z_out, MatrixXd &S_out, MatrixXd &Tc_out) 
{
	//create matrix for sigma points in measurement space
	static MatrixXd Zsig = MatrixXd(n_radar_, n_sig_points_);

	//transform sigma points into measurement space
	for (int i = 0; i < n_sig_points_; i++)
	{ 
		// extract values for better readibility
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double yaw = Xsig_pred_(3, i);

		// validate division by zero
		if (fabs(p_x) <= 0.0001) p_x = 0.0001;
		if (fabs(p_y) <= 0.0001) p_y = 0.0001;

		double v1 = cos(yaw)*v;
		double v2 = sin(yaw)*v;

		// measurement model
		Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);                       //r
		Zsig(1, i) = atan2(p_y, p_x);                               //phi
		Zsig(2, i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
	}

	//mean predicted measurement
	static VectorXd z_pred = VectorXd(n_radar_);
	z_pred.fill(0.0);
	for (int i = 0; i < n_sig_points_; i++) 
		z_pred = z_pred + weights_(i) * Zsig.col(i);

	//measurement covariance matrix S
	static MatrixXd S = MatrixXd(n_radar_, n_radar_);
	S.fill(0.0);
	//create matrix for cross correlation Tc
	static MatrixXd Tc = MatrixXd(n_x_, n_radar_);
	Tc.fill(0.0);

	for (int i = 0; i < n_sig_points_; i++)
	{
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
	S = S + R_radar;

	//write result
	z_out = z_pred;
	S_out = S;
	Tc_out = Tc;

	return;
}

void UKF::PredictLidarMeasurement(VectorXd &z_out, MatrixXd &S_out, MatrixXd &Tc_out) 
{
	// matrix for sigma points in measurement space
	static MatrixXd Zsig = MatrixXd(n_laser_, n_sig_points_);

	//transform sigma points into measurement space
	for (int i = 0; i < n_sig_points_; i++)
	{
		Zsig(0, i) = Xsig_pred_(0, i);          //px
		Zsig(1, i) = Xsig_pred_(1, i);          //py
	}

	//mean predicted measurement
	static VectorXd z_pred = VectorXd(n_laser_);
	z_pred.fill(0.0);
	for (int i = 0; i < n_sig_points_; i++)
		z_pred = z_pred + weights_(i) * Zsig.col(i);

	//measurement covariance matrix S
	static MatrixXd S = MatrixXd(n_laser_, n_laser_);
	S.fill(0.0);

	//create matrix for cross correlation Tc
	static MatrixXd Tc = MatrixXd(n_x_, n_laser_);
	Tc.fill(0.0);

	for (int i = 0; i < n_sig_points_; i++) 
	{
		VectorXd z_diff = Zsig.col(i) - z_pred;
		S = S + weights_(i) * z_diff * z_diff.transpose();

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}

	//add measurement noise covariance matrix
	S = S + R_lidar;

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
	static MatrixXd Xsig_aug = MatrixXd(n_sig_points_, n_aug_);

	//create example matrix with predicted sigma points
	static MatrixXd Xsig_pred = MatrixXd(n_x_, n_sig_points_);

	// compute augmented sigma points
	AugmentedSigmaPoints(Xsig_aug);

	// predict augmented sigma points
	SigmaPointPrediction(Xsig_aug, delta_t, Xsig_pred);

	// predicts the mean and covariance
	PredictMeanAndCovariance(Xsig_pred);

	return;
}

/**
* Updates the state and the state covariance matrix using a laser measurement.
* @param {MeasurementPackage} meas_package
*/
void UKF::UpdateLidar(MeasurementPackage meas_package, VectorXd &z_pred, MatrixXd &Tc, MatrixXd &S) 
{
	//mean predicted measurement
	static VectorXd z = VectorXd::Zero(n_laser_);
	z << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1);

	//Kalman gain K;
	static MatrixXd K = MatrixXd::Zero(n_x_, n_laser_);
	K = Tc * S.inverse();

	//residual
	static VectorXd z_diff = VectorXd::Zero(n_laser_);
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
	//mean predicted measurement  rho, phi, rho_dot
	static VectorXd z = VectorXd::Zero(n_radar_);
	z << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), meas_package.raw_measurements_(2);

	//Kalman gain K;
	static MatrixXd K = MatrixXd::Zero(n_x_, n_radar_);
	K = Tc * S.inverse();

	//residual
	static VectorXd z_diff = VectorXd::Zero(n_radar_);
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
