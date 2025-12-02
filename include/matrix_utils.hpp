#pragma once
#include <Eigen/Dense>
#include <complex>
#include <vector>

double eigen_error_norm(const std::vector<std::complex<double>> &exact,
                        const Eigen::VectorXd &alphar,
                        const Eigen::VectorXd &alphai,
                        const Eigen::VectorXd &beta);
