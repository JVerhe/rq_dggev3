#include "pencil_generator.hpp"
#include <Eigen/Dense>
#include <random>
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
