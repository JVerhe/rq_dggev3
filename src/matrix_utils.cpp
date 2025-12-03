#include <vector>
#include <complex>
#include <limits>
#include <cmath>
#include <stdexcept>
#include <Eigen/Dense>
#include <iostream>

static constexpr double INF_MISMATCH_PENALTY = 1e16;

bool isComplexInf(const std::complex<double> &z)
{
    return std::isinf(z.real()) || std::isinf(z.imag());
}

bool isComplexFinite(const std::complex<double> &z)
{
    return std::isfinite(z.real()) || std::isfinite(z.imag());
}

double eigErrorNorm(
    const std::vector<std::complex<double>> &reference,
    const Eigen::VectorXd &alphar,
    const Eigen::VectorXd &alphai,
    const Eigen::VectorXd &beta)
{
    const int N = static_cast<int>(reference.size());
    if (alphar.size() != N || alphai.size() != N || beta.size() != N)
        throw std::invalid_argument("Eigenvalue dimension mismatch in eigen_error_norm.");

    // Construct computed eigenvalues
    std::vector<std::complex<double>> computed(N);
    for (int i = 0; i < N; ++i)
    {
        if (beta(i) == 0.0) // Tolerance for inifinite vectors?
        {
            computed[i] = std::complex<double>(std::numeric_limits<double>::infinity(), 0.0);
        }
        else
        {
            computed[i] = std::complex<double>(alphar(i), alphai(i)) / beta(i);
        }
    }

    // List of indices for reference entries that are valid (not NaN)
    std::vector<int> ref_idx;
    ref_idx.reserve(N);
    for (int i = 0; i < N; ++i)
    {
        const auto &r = reference[i];
        if (std::isnan(r.real()) || std::isnan(r.imag()))
        {
            continue;
        }
        ref_idx.push_back(i);
    }
    if (ref_idx.empty()) // Only fake eigenvalues
    {
        return std::numeric_limits<double>::quiet_NaN(); // nothing to compare
    }

    auto mag_ang = [&](const std::complex<double> &z)
    {
        double m = std::abs(z);
        double a = std::arg(z);
        return std::make_pair(m, a);
    };

    std::sort(ref_idx.begin(), ref_idx.end(), [&](int i, int j)
              { return mag_ang(reference[i]) < mag_ang(reference[j]); });

    std::vector<char> used(N, 0); // flags for computed entries
    double error_sum = 0.0;
    int matched_count = 0;

    for (int ridx : ref_idx) // loop over valid references
    {
        const std::complex<double> &r = reference[ridx];

        // Special-case: reference infinite -> try to match an infinite computed first
        if (std::isinf(r.real()) || std::isinf(r.imag()))
        {
            int found_inf = -1;
            for (int j = 0; j < N; ++j) // loop over computed values
            {
                if (used[j])
                    continue;
                const auto &c = computed[j];
                if (std::isinf(c.real()) || std::isinf(c.imag()))
                {
                    found_inf = j;
                    break;
                }
            }
            if (found_inf != -1)
            {
                used[found_inf] = 1;
                matched_count++; // No penalty is added
                continue;
            }
        }

        // General case: find best unused computed index minimizing distance (with infinite handling)
        int best_j = -1;
        double best_cost = std::numeric_limits<double>::infinity();

        for (int j = 0; j < N; ++j)
        {
            if (used[j])
                continue;
            const auto &c = computed[j];

            if (std::isnan(c.real()) || std::isnan(c.imag())) // Safety check, nan normally already filtered
                continue;

            // If reference finite but computed infinite -> large penalty cost
            if ((isComplexInf(c)) && (isComplexFinite(r)))
            {
                // express as large constant cost (but allow selection if no better)
                if (INF_MISMATCH_PENALTY < best_cost)
                {
                    best_cost = INF_MISMATCH_PENALTY;
                    best_j = j;
                }
                continue;
            }

            // If reference infinite but computed finite -> large penalty
            if (isComplexInf(r) && isComplexFinite(c))
            {
                if (INF_MISMATCH_PENALTY < best_cost)
                {
                    best_cost = INF_MISMATCH_PENALTY;
                    best_j = j;
                }
                continue;
            }

            double dist = std::abs((c - r) / r);
            double cost = dist;
            if (cost < best_cost)
            {
                best_cost = cost;
                best_j = j;
            }
        }
        if (best_j == -1)
        {
            throw std::runtime_error("Couldn't find a match");
        }

        used[best_j] = 1;
        const auto &best_c = computed[best_j];

        // infinite-infinite -> no error
        if (isComplexInf(r) && isComplexInf(best_c))
        {
            matched_count++;
            continue;
        }

        // Reference Finite, computed is infinite
        if (isComplexFinite(r) && isComplexInf(best_c))
        {
            double squared_penalty = 0.0;
            if (alphar[best_j] != 0.0)
            {
                squared_penalty = (r * r).real() / (alphar[best_j] * alphar[best_j]);
            }
            else
            {
                squared_penalty = (r * r).real();
            }
            error_sum += squared_penalty;
            matched_count++;
            continue;
        }

        // Reference Infintie, computed finite
        if (isComplexInf(r) && isComplexFinite(best_c))
        {
            double squared_penalty = beta[best_j] * beta[best_j] / (alphar[best_j] * alphar[best_j]);
            error_sum += squared_penalty;
            matched_count++;
            continue;
        }

        // Both finite:
        double penalty = std::abs((best_c - r) / r);
        error_sum += penalty * penalty;
        matched_count++;
    }

    if (matched_count == 0)
        return std::numeric_limits<double>::quiet_NaN();

    return std::sqrt(error_sum);
}

double invEigenErrorNorm(const std::vector<std::complex<double>> &reference,
                         const Eigen::VectorXd &alphar,
                         const Eigen::VectorXd &alphai,
                         const Eigen::VectorXd &beta)
{
    // Calculate inverse of eigs
    std::vector<std::complex<double>> inv_ref;
    inv_ref.reserve(reference.size());
    for (const auto &z : reference)
    {
        inv_ref.push_back(1.0 / z);
    }

    // for now, all eigenvalues are real -> switch beta and alpha_r
    return eigErrorNorm(inv_ref, beta, alphai, alphar);
}
