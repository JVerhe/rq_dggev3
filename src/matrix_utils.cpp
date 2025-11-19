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

    std::vector<std::complex<double>> computed;
    computed.reserve(N);
    for (int i = 0; i < N; ++i)
    {
        if (beta(i) == 0.0)
        {
            // True infinite eigenvalue
            computed.emplace_back(std::numeric_limits<double>::infinity(), 0.0);
        }
        else
        {
            computed.emplace_back(
                std::complex<double>(alphar(i), alphai(i)) / beta(i));
        }
    }

    std::vector<std::complex<double>> ref_sorted = reference;
    auto eigen_sort_key = [](const std::complex<double> &z)
    {
        double mag = std::abs(z);
        double ang = std::arg(z);

        // Fake singular (NaN) should go last
        if (std::isnan(z.real()) || std::isnan(z.imag()))
            return std::make_tuple(std::numeric_limits<double>::infinity(),
                                   std::numeric_limits<double>::infinity());

        // Infinite eigenvalues → magnitude = +∞, angle = ±π (keep consistent)
        if (std::isinf(z.real()) || std::isinf(z.imag()))
            return std::make_tuple(std::numeric_limits<double>::infinity(), ang);

        return std::make_tuple(mag, ang);
    };

    auto cmp = [&](const std::complex<double> &a, const std::complex<double> &b)
    {
        return eigen_sort_key(a) < eigen_sort_key(b);
    };

    std::sort(ref_sorted.begin(), ref_sorted.end(), cmp);
    std::sort(computed.begin(), computed.end(), cmp);

    double error_sum = 0.0;
    int valid_count = 0;

    for (int i = 0; i < N; ++i)
    {
        const auto &ref = ref_sorted[i];
        const auto &comp = computed[i];

        // Fake singular (NaN) on reference -> ignore
        if (std::isnan(ref.real()) || std::isnan(ref.imag()))
            continue;

        // Both infinite ->  zero error
        if (std::isinf(ref.real()) && std::isinf(comp.real()))
            continue;

        // One infinite, one finite -> heavy penalty
        if (std::isinf(ref.real()) != std::isinf(comp.real()))
        {
            error_sum += 1e6;
            valid_count++;
            continue;
        }

        // Both finite → normal error
        double diff = std::abs(comp - ref);
        error_sum += diff * diff;
        valid_count++;
    }

    if (valid_count == 0)
        return std::numeric_limits<double>::quiet_NaN();

    return std::sqrt(error_sum / valid_count);
}