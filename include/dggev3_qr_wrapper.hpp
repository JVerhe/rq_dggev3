#pragma once
#include <Eigen/Dense>

struct GeneralizedEigenResult
{
    Eigen::VectorXd alphar;
    Eigen::VectorXd alphai;
    Eigen::VectorXd beta;
    Eigen::MatrixXd VL;
    Eigen::MatrixXd VR;
    int info;
};

GeneralizedEigenResult dggev3_qr_wrapper(
    bool computeVL,
    bool computeVR,
    Eigen::MatrixXd A,
    Eigen::MatrixXd B);
