#include "matrix_utils.hpp"
#include <limits>
#include <cmath>
#include <stdexcept>
#include <complex>
#include <vector>
#include <Eigen/Dense>

double eigen_error_norm(
    const std::vector<std::complex<double>> &reference,
    const Eigen::VectorXd &alphar,
    const Eigen::VectorXd &alphai,
    const Eigen::VectorXd &beta)
{
    const int N = static_cast<int>(reference.size());
    if (alphar.size() != N || alphai.size() != N || beta.size() != N)
        throw std::invalid_argument("Eigenvalue dimension mismatch in eigen_error_norm.");

    double error_sum = 0.0;
    int valid_count = 0;

    for (int i = 0; i < N; ++i)
    {
        const auto &ref = reference[i];

        // Ignore fake singular (NaN) eigenvalues
        if (std::isnan(ref.real()) || std::isnan(ref.imag()))
            continue;

        std::complex<double> comp;
        if (beta(i) == 0.0)
        {
            // True infinite eigenvalue
            comp = std::complex<double>(std::numeric_limits<double>::infinity(), 0.0);
        }
        else
        {
            comp = std::complex<double>(alphar(i), alphai(i)) / beta(i);
        }

        // Handle infinite cases
        if (std::isinf(ref.real()))
        {
            if (std::isinf(comp.real()))
            {
                continue; // both infinite → no error
            }
            else
            {
                error_sum += 1e6; // penalty for finite vs infinite mismatch
                valid_count++;
                continue;
            }
        }

        if (std::isinf(comp.real()) && !std::isinf(ref.real()))
        {
            error_sum += 1e6;
            valid_count++;
            continue;
        }

        // Both finite → normal difference
        double diff = std::abs(comp - ref);
        error_sum += diff * diff;
        valid_count++;
    }

    if (valid_count == 0)
        return std::numeric_limits<double>::quiet_NaN();

    return std::sqrt(error_sum / valid_count);
}