#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
    x_ = F_ * x_;
    P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */

    VectorXd y;
    MatrixXd S;

    //Calculate the innovation estimate
    y = z - H_ * x_;

    //Calculate the innovation covariance
    S = H_ * P_ * H_.transpose() + R_;

    //Calculate Kalman Gain
    MatrixXd K = P_ * H_.transpose() * S.inverse();

    //Calculate the fused state expectation
    x_ = x_ + K * y;

    //Calculate the fused state covariance
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);

    MatrixXd temp = (I - K * H_)*P_;

    P_ = temp;
}

void KalmanFilter::UpdateEKF(const VectorXd &z, const VectorXd &h_x) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

    VectorXd y;
    MatrixXd S;

    //Calculate the innovation estimate
    y = z - h_x;

    // Normalize the angle difference to be between -PI and PI
    while(y(1) > M_PI) y(1) -= 2. * M_PI;
    while(y(1) < -M_PI) y(1) += 2. * M_PI;

    //Calculate the innovation covariance
    S = H_ * P_ * H_.transpose() + R_;

    //Calculate Kalman Gain
    MatrixXd K = P_ * H_.transpose() * S.inverse();

    //Calculate the fused state expectation
    x_ = x_ + K * y;

    //Calculate the fused state covariance
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);

    MatrixXd temp = (I - K * H_)*P_;

    P_ = temp;
}
