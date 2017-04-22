#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

ofstream ofs("debug.txt");

#define MIN_SENSOR_VALUE    0.000001

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
    if (!ofs)
        return;

    is_initialized_ = false;

    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    // initial state vector
    x_ = VectorXd(5);

    // initial covariance matrix
    P_ = MatrixXd(5, 5);

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 2.4;   // Need to tuning (ninety-five percent of the time)

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = .6;   // Need to tuning (ninety-five percent of the time)

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

    /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
    // Init some parameters
    n_x_ = 5;
    n_aug_= 7;
    lambda_ = 3 - n_aug_;
    // NIS_radar_? and NIS_laser_?
    // Init P_
    // [pos1 pos2 vel_abs yaw_angle yaw_rate]
//    P_ << 1., 0, 0, 0, 0,
//            0, 1, 0, 0, 0,
//            0, 0, 1, 0, 0,
//            0, 0, 0, 1, 0,
//            0, 0, 0, 0, 1;
    P_ << .5, 0, 0, 0, 0,
            0, .5, 0, 0, 0,
            0, 0, .25, 0, 0,
            0, 0, 0, .25, 0,
            0, 0, 0, 0, .25;
    // x_ will be initialized in ProcessMeasurement() based on reading sensor data
    // Init Xsig_pred_
    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
    Xsig_pred_.fill(0.0);
    // Init weights_
    weights_ = VectorXd(2*n_aug_+1);
    // set weights
    double weight_0 = lambda_/(lambda_+n_aug_);
    weights_(0) = weight_0;
    for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
        double weight = 0.5/(n_aug_+lambda_);
        weights_(i) = weight;
    }
}

UKF::~UKF() {
    ofs.close();
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

    if (use_radar_==false && (meas_package.sensor_type_ == MeasurementPackage::RADAR))
        return;
    if (use_laser_==false && (meas_package.sensor_type_ == MeasurementPackage::LASER))
        return;
    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) {
        /**
      TODO:
        * Initialize the state x_ with the first measurement.
        * Create the covariance matrix.
        * Remember: you'll need to convert radar from polar to cartesian coordinates.
      */
        // [pos1 pos2 vel_abs yaw_angle yaw_rate]
        x_ = VectorXd(5);

        time_us_ = meas_package.timestamp_;
        // Initializing x_
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            /**
            Convert radar from polar to CTRV coordinates.
            */
            float rho = meas_package.raw_measurements_(0);
            float theta = meas_package.raw_measurements_(1);
            theta=-theta;
            x_ << cos(theta)*rho, sin(theta)*rho, abs(meas_package.raw_measurements_(2)), theta+M_PI/2, 0;
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            /**
            Initialize state.
            */
            x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 1, M_PI, 0;
        }
        // Avoid initializing nonsence data
        if (x_[0] < MIN_SENSOR_VALUE && x_[1] < MIN_SENSOR_VALUE)
        {
            return;
        }

        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/
    double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
    time_us_ = meas_package.timestamp_;

    Prediction(dt);

    /*****************************************************************************
     *  Update
     ****************************************************************************/

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        // Radar updates
        UpdateRadar(meas_package);
    } else {
        // Laser updates
        UpdateLidar(meas_package);
    }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
//    std::cout << "Prediction: delta_t = " << delta_t << std::endl;
    // generating sigma points
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    AugmentedSigmaPoints(&Xsig_aug);

    Xsig_pred_.fill(0.0);

    //predict sigma points
    for (int i = 0; i< 2*n_aug_+1; i++)
    {
        //extract values for better readability
        double p_x = Xsig_aug(0,i);
        double p_y = Xsig_aug(1,i);
        double v = Xsig_aug(2,i);
        double yaw = Xsig_aug(3,i);

        //        // [Normalize]?
        //        yaw = angleNormalize(yaw);

        double yawd = Xsig_aug(4,i);
        double nu_a = Xsig_aug(5,i);
        double nu_yawdd = Xsig_aug(6,i);

        //predicted state values
        double px_p, py_p, v_p, yaw_p, yawd_p;

        if (delta_t<0.0001)
        {
            px_p = p_x;
            py_p = p_y;
            v_p = v;
            yaw_p = yaw;
            yawd_p = yawd;
        }
        else
        {
            //avoid division by zero
            if (fabs(yawd) > 0.001) {
                px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
                py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
            }
            else {
                px_p = p_x + v*delta_t*cos(yaw);
                py_p = p_y + v*delta_t*sin(yaw);
            }

            v_p = v;
            yaw_p = yaw + yawd*delta_t;
            yawd_p = yawd;

            //add noise
            px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
            py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
            v_p = v_p + nu_a*delta_t;

            yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
            yawd_p = yawd_p + nu_yawdd*delta_t;
        }

        //write predicted sigma point into right column
        Xsig_pred_(0,i) = px_p;
        Xsig_pred_(1,i) = py_p;
        Xsig_pred_(2,i) = v_p;
        Xsig_pred_(3,i) = yaw_p;
        //        // [Normalize]?
        //        Xsig_pred_(3,i) = angleNormalize(Xsig_pred_(3,i);
        Xsig_pred_(4,i) = yawd_p;

        if (delta_t>=0.0001)
        {
            // x_ and P_ become new predicted values
            PredictMeanAndCovarianceStateSpace();
        }
    }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
    float min_sensor_value_ = 0.0001;
    // processing nonsense data (both data are too small)
    if (x_(0)<=min_sensor_value_ && x_(1)<=min_sensor_value_)
    {
        x_(0) = min_sensor_value_;
        x_(1) = min_sensor_value_;
    }

    // measurement matrix
    int n_z = 2;
    Eigen::MatrixXd H_ = MatrixXd(n_z,n_x_);
    H_ << 1, 0, 0, 0, 0,
            0, 1, 0, 0, 0;
    VectorXd z_pred = H_ * x_;
    VectorXd z = meas_package.raw_measurements_;
    VectorXd y = z - z_pred;
    MatrixXd Ht = H_.transpose();
    MatrixXd PHt = P_ * Ht;
    //add measurement noise covariance matrix
    MatrixXd R_ = MatrixXd(n_z,n_z);
    R_ <<    std_laspx_*std_laspx_, 0,
            0, std_laspy_*std_laspy_;
    MatrixXd S = H_ * PHt + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K = PHt * Si;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;

    ofs<<"Lidar Update x_ = "<< x_ << "\nP_ = \n"<<P_ << std::endl; //[P_]
    //    std::cout<<"Lidar Update P_ = \n"<<P_ << std::endl; //[P_]

    NIS_laser_ = y.transpose()*Si*y;
    int tmp=0;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
    // Transform Xsig_pred_ to measurement space and calculate mean and covariance
    VectorXd z_pred_;    // mean
    MatrixXd S_;    // cov
    MatrixXd Zsig_; // z sigma points
    PredictRadarMeasurement(&z_pred_, &S_, &Zsig_);

    // Kamlman Filter update
    UpdateStateRadar(z_pred_, S_, Zsig_, meas_package.raw_measurements_);

}

void UKF::AugmentedSigmaPoints(Eigen::MatrixXd *Xsig_out)
{
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
    P_aug.topLeftCorner(5,5) = P_;
    P_aug(5,5) = std_a_*std_a_;
    P_aug(6,6) = std_yawdd_*std_yawdd_;

    //create square root matrix
    MatrixXd L = P_aug.llt().matrixL();

    //create augmented sigma points
    Xsig_aug.col(0)  = x_aug;
    for (int i = 0; i< n_aug_; i++)
    {
        Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
    }

    //write result
    *Xsig_out = Xsig_aug;
}

void UKF::PredictMeanAndCovarianceStateSpace()
{
    //create vector for predicted state
    VectorXd x = VectorXd(n_x_);

    //create covariance matrix for prediction
    MatrixXd P = MatrixXd(n_x_, n_x_);

    //predicted state mean
    x.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        x = x+ weights_(i) * Xsig_pred_.col(i);
    }

    //predicted state covariance matrix
    P.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x;
        //angle normalization
        x_diff(3) = angleNormalize(x_diff(3));

        P = P + weights_(i) * x_diff * x_diff.transpose() ;
    }
    //write result
    x_ = x;
    P_ = P;
    ofs<<"State space recurtion Prediction:\nx_="<<x<<"\nP_="<<P_ << std::endl; //[P_]
    //    std::cout<<P_ << std::endl; //[P_]
    int tmp=0;
}

void UKF::PredictRadarMeasurement(Eigen::VectorXd *z_mean_out, Eigen::MatrixXd *S_out, Eigen::MatrixXd *Zsig_out)
{
    //set measurement dimension, radar can measure r, phi, and r_dot
    int n_z = 3;

    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

        // extract values for better readibility
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);
        double v  = Xsig_pred_(2,i);
        double yaw = Xsig_pred_(3,i);
        // avoid numerical issue ??
        p_x = (fabs(p_x)<=MIN_SENSOR_VALUE) ? MIN_SENSOR_VALUE : p_x;
        p_y = (fabs(p_y)<=MIN_SENSOR_VALUE) ? MIN_SENSOR_VALUE : p_y;
        //        //angle normalization?
        //        yaw = angleNormalize(yaw);

        double v1 = cos(yaw)*v;
        double v2 = sin(yaw)*v;

        // measurement model
        Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
        Zsig(1,i) = atan2(p_y,p_x);                                 //phi
        Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
    }

    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i=0; i < 2*n_aug_+1; i++) {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }

    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        //angle normalization
        z_diff(1) = angleNormalize(z_diff(1));

        S = S + weights_(i) * z_diff * z_diff.transpose();
    }

    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z,n_z);
    R <<    std_radr_*std_radr_, 0, 0,
            0, std_radphi_*std_radphi_, 0,
            0, 0,std_radrd_*std_radrd_;
    S = S + R;

    //write result
    *z_mean_out = z_pred;
    *S_out = S;
    *Zsig_out = Zsig;
}

void UKF::UpdateStateRadar(const VectorXd &z_pred, const MatrixXd &S, const MatrixXd &Zsig, const VectorXd &z)
{
    //set measurement dimension, radar can measure r, phi, and r_dot
    int n_z = 3;

    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        //angle normalization
        z_diff(1) = angleNormalize(z_diff(1));

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        x_diff(3) = angleNormalize(x_diff(3));

        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    //Kalman gain K;
    MatrixXd Si = S.inverse();
    MatrixXd K = Tc * Si;

    //residual
    VectorXd z_diff = z - z_pred;

    //angle normalization
    z_diff(1) = angleNormalize(z_diff(1));

    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    x_(3) = angleNormalize(x_(3));

    P_ = P_ - K*S*K.transpose();
    ofs<<"Radar Update x_ = "<<x_<<"\nP_ = \n"<<P_ << std::endl; //[P_]
    //    std::cout<<"Radar Update P_ = \n"<<P_ << std::endl; //[P_]

    // Calculate the NIS
    NIS_radar_ = z_diff.transpose()*Si*z_diff;
    int tmp=0;
}

float UKF::angleNormalize(float x)
{
    //    while (x> M_PI) x-=2.*M_PI;
    //    while (x<-M_PI) x+=2.*M_PI;
    if (fabs(x) > M_PI)
    {
        x -= round(x / (2.0 * M_PI)) * (2.0 * M_PI);
    }
    return x;
}
