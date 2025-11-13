#include "pencil_generator.hpp"
#include <Eigen/Dense>
#include <random>
#include <limits>
#include <cmath>
#include <stdexcept>

using namespace Eigen;

Pencil generate_regular_pencil(int N)
{
    if (N <= 0)
        throw std::invalid_argument("Matrix dimension N must be positive");

    std::mt19937 rng(std::random_device{}());
    std::normal_distribution<> dist(0.0, 1.0);

    MatrixXd X = MatrixXd::NullaryExpr(N, N, [&]()
                                       { return dist(rng); });
    MatrixXd Y = MatrixXd::NullaryExpr(N, N, [&]()
                                       { return dist(rng); });

    HouseholderQR<MatrixXd> qrX(X);
    HouseholderQR<MatrixXd> qrY(Y);

    MatrixXd Qx = qrX.householderQ() * MatrixXd::Identity(N, N);
    MatrixXd Qy = qrY.householderQ() * MatrixXd::Identity(N, N);

    // Diagonal matrix D = diag(1, 2, ..., N)
    MatrixXd D = MatrixXd::Zero(N, N);
    for (int i = 0; i < N; ++i)
        D(i, i) = static_cast<double>(i + 1);

    // Build matrices A and B
    MatrixXd A = Qx * D * Qy;
    MatrixXd B = Qx * Qy;

    // Exact eigenvalues
    std::vector<std::complex<double>> eigvals;
    eigvals.reserve(N);
    for (int i = 0; i < N; ++i)
        eigvals.emplace_back(static_cast<double>(i + 1), 0.0);

    return {A, B, eigvals};
}

Pencil generate_singular_pencil(int N)
{
    if (N <= 1)
        throw std::invalid_argument("Matrix dimension N must be at least 2 for a singular pencil.");

    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<> dist(0.5, 2.0); // avoid near-zero random numbers

    MatrixXd A = MatrixXd::Zero(N, N);
    MatrixXd B = MatrixXd::Zero(N, N);

    // Generate random upper-triangular matrices
    for (int i = 0; i < N; ++i)
    {
        for (int j = i; j < N; ++j)
        {
            A(i, j) = dist(rng);
            B(i, j) = dist(rng);
        }
    }

    // Zero out some diagonal entries in both A and B
    int num_zeros = std::max(1, N / 10);
    std::uniform_int_distribution<> diag_dist(0, N - 1);

    for (int k = 0; k < num_zeros; ++k)
    {
        int idx = diag_dist(rng);
        A(idx, idx) = 0.0;
        B(idx, idx) = 0.0;
    }

    // Compute “exact” eigenvalues
    std::vector<std::complex<double>> eigvals;
    eigvals.reserve(N);

    for (int i = 0; i < N; ++i)
    {
        if (A(i, i) == 0.0 && B(i, i) == 0.0)
        {
            // undefined (fake singular)
            eigvals.emplace_back(std::numeric_limits<double>::quiet_NaN(), 0.0);
        }
        else if (B(i, i) == 0.0)
        {
            // infinite eigenvalue (well-defined)
            eigvals.emplace_back(std::numeric_limits<double>::infinity(), 0.0);
        }
        else
        {
            eigvals.emplace_back(A(i, i) / B(i, i), 0.0);
        }
    }

    return {A, B, eigvals};
}
