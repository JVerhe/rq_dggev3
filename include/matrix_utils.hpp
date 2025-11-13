#pragma once
#include <Eigen/Dense>
#include <complex>
#include <vector>

// Compute Euclidean (L2) norm of the difference between computed and exact eigenvalues
double eigen_error_norm(const std::vector<std::complex<double>> &exact,
                        const Eigen::VectorXd &alphar,
                        const Eigen::VectorXd &alphai,
                        const Eigen::VectorXd &beta);
