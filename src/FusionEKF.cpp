#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  ekf_.Q_ = MatrixXd(4,4);


  //the initial transition matrix F_
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
            0, 1, 0, 1,
            0, 0, 1, 0,
            0, 0, 0, 1;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


    /*****************************************************************************
   *  Initialization
   ****************************************************************************/
    if (!is_initialized_) {
        /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
        // first measurement
        cout << "EKF: " << endl;
        ekf_.x_ = VectorXd(4);
        ekf_.x_ << 1, 1, 5, 1;

        //Initialize the state covariance
        ekf_.P_ = MatrixXd(4, 4);
        ekf_.P_ << 1., 0, 0, 0,
                0, 1., 0, 0,
                0, 0, 1000., 0,
                0, 0, 0, 1000.;

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
        {
            /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */

            float rho = measurement_pack.raw_measurements_(0);
            float phi = measurement_pack.raw_measurements_(1);

            ekf_.x_(0) = rho * cos(phi);
            ekf_.x_(1) = rho * sin(phi);

        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
        {
            /**
      Initialize state.
      */

            ekf_.x_(0) = measurement_pack.raw_measurements_(0);
            ekf_.x_(1) = measurement_pack.raw_measurements_(1);
        }

        previous_timestamp_ = measurement_pack.timestamp_;

        // done initializing, no need to predict or update
        is_initialized_ = true;

        return;


    }

    /*****************************************************************************
   *  Prediction
   ****************************************************************************/

    /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
     * Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

    //compute the time elapsed between the current and previous measurements
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
    previous_timestamp_ = measurement_pack.timestamp_;

    ekf_.F_ << 1, 0, dt, 0,
            0, 1, 0, dt,
            0, 0, 1, 0,
            0, 0, 0, 1;

    //set the acceleration noise components
    float noise_ax = 9;
    float noise_ay = 9;
    float dt_2 = dt * dt;
    float dt_3 = dt * dt_2;
    float dt_4 = dt * dt_3;

    ekf_.Q_ << ((dt_4* noise_ax)/4 ), 0, ((dt_3* noise_ax)/2 ), 0,
            0, ((dt_4* noise_ay)/4 ), 0, ((dt_3* noise_ay)/2 ),
            ((dt_3* noise_ax)/2 ), 0, (dt_2* noise_ax) , 0,
            0, ((dt_3* noise_ay)/2 ), 0, (dt_2* noise_ay);

    ekf_.Predict();

    /*****************************************************************************
   *  Update
   ****************************************************************************/

    /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
        // Radar updates

        ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
        ekf_.R_ = R_radar_;


        //Calculate the measurement funcation
        VectorXd h_x(3);

        float rho = sqrt(ekf_.x_(0) * ekf_.x_(0) + ekf_.x_(1) * ekf_.x_(1));

        //check division by zero
        if(fabs(rho) < 0.0001)
        {

            cout << "ProcessMeasurement () - Error - Division by Zero" << endl;

            h_x(0) = rho;
            h_x(1) = atan2(ekf_.x_(1), ekf_.x_(0));
            h_x(2) = 0;
        }
        else
        {
            h_x(0) = rho;
            h_x(1) = atan2(ekf_.x_(1), ekf_.x_(0));
            h_x(2) = (ekf_.x_(0) * ekf_.x_(2) + ekf_.x_(1) * ekf_.x_(3)) / rho;
        }


        ekf_.UpdateEKF(measurement_pack.raw_measurements_, h_x);

    }
    else
    {
        // Laser updates

        ekf_.H_ = MatrixXd(2, 4);
        ekf_.H_ << 1, 0, 0, 0,
                0, 1, 0, 0;

        ekf_.R_ = R_laser_;
        ekf_.Update(measurement_pack.raw_measurements_);
    }

    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}
