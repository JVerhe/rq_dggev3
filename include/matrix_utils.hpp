#pragma once
#include <Eigen/Dense>
#include <complex>
#include <vector>

bool isComplexInfinite(const std::complex<double> &z);

bool isComplexFinite(const std::complex<double> &z);

double eigen_error_norm(const std::vector<std::complex<double>> &exact,
                        const Eigen::VectorXd &alphar,
                        const Eigen::VectorXd &alphai,
                        const Eigen::VectorXd &beta);
