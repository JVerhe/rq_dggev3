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

Pencil generate_singular_triangular_pencil(int N)
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

Pencil generate_singular_pencil(int N)
{
    if (N < 2)
        throw std::invalid_argument("N must be >= 2");

    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<> dist_val(0.5, 2.0);
    std::uniform_int_distribution<> diag_pick(0, N - 1);

    // --- Step 1: Create diagonal DA, DB with singular structure ---
    Eigen::VectorXd diagA(N), diagB(N);

    for (int i = 0; i < N; i++)
    {
        diagA(i) = dist_val(rng);
        diagB(i) = dist_val(rng);
    }

    // Force some singularities
    int k = std::max(1, N / 10);

    for (int i = 0; i < k; i++)
    {
        int idx = diag_pick(rng);

        if (rng() & 1)
        {
            // produce fake singular (undefined eigenvalue)
            diagA(idx) = 0.0;
            diagB(idx) = 0.0;
        }
        else
        {
            // produce true infinite eigenvalue
            diagB(idx) = 0.0;
            // diagA stays random and nonzero
        }
    }

    Eigen::MatrixXd DA = diagA.asDiagonal();
    Eigen::MatrixXd DB = diagB.asDiagonal();

    // --- Step 2: Generate two random orthogonal matrices U, V ---
    Eigen::MatrixXd R1 = Eigen::MatrixXd::Random(N, N);
    Eigen::MatrixXd R2 = Eigen::MatrixXd::Random(N, N);

    Eigen::HouseholderQR<Eigen::MatrixXd> qr1(R1);
    Eigen::HouseholderQR<Eigen::MatrixXd> qr2(R2);

    Eigen::MatrixXd U = qr1.householderQ() * Eigen::MatrixXd::Identity(N, N);
    Eigen::MatrixXd V = qr2.householderQ() * Eigen::MatrixXd::Identity(N, N);

    // --- Step 3: Construct dense singular A and B ---
    Eigen::MatrixXd A = U * DA * V.transpose();
    Eigen::MatrixXd B = U * DB * V.transpose();

    // --- Step 4: Compute reference eigenvalues ---
    std::vector<std::complex<double>> eigvals;
    eigvals.reserve(N);

    for (int i = 0; i < N; i++)
    {
        double a = diagA(i);
        double b = diagB(i);

        if (a == 0.0 && b == 0.0)
        {
            eigvals.emplace_back(std::numeric_limits<double>::quiet_NaN(), 0.0);
        }
        else if (b == 0.0)
        {
            eigvals.emplace_back(std::numeric_limits<double>::infinity(), 0.0);
        }
        else
        {
            eigvals.emplace_back(a / b, 0.0);
        }
    }

    return {A, B, eigvals};
}

Pencil generate_illconditioned_B_pencil(int N,
                                        bool use_integer_DA)
{
    if (N <= 0)
        throw std::invalid_argument("N must be positive");

    // Random generator
    std::mt19937_64 rng(std::random_device{}());
    std::normal_distribution<double> nd(0.0, 1.0);

    // --- 1) Build diagonal DB = [10^(base_exponent), 10^(base_exponent+1), ...] ---
    VectorXd diagB(N);
    for (int i = 0; i < N; ++i)
    {
        double exp = -16 * i / (N - 1); // e.g. -16, -15, ...
        diagB(i) = std::pow(10.0, exp);
    }

    // --- 2) Build diagonal DA ---
    VectorXd diagA(N);
    if (use_integer_DA)
    {
        for (int i = 0; i < N; ++i)
            diagA(i) = static_cast<double>(i + 1);
    }
    else
    {
        // Optionally you could use random positive values
        for (int i = 0; i < N; ++i)
            diagA(i) = std::abs(nd(rng)) + 0.5;
    }

    MatrixXd DA = diagA.asDiagonal();
    MatrixXd DB = diagB.asDiagonal();

    MatrixXd R1(N, N), R2(N, N);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
        {
            R1(i, j) = nd(rng);
            R2(i, j) = nd(rng);
        }

    HouseholderQR<MatrixXd> qr1(R1);
    HouseholderQR<MatrixXd> qr2(R2);
    MatrixXd U = qr1.householderQ() * MatrixXd::Identity(N, N);
    MatrixXd V = qr2.householderQ() * MatrixXd::Identity(N, N);

    MatrixXd A = U * DA * V.transpose();
    MatrixXd B = U * DB * V.transpose();

    std::vector<std::complex<double>> eigs;
    eigs.reserve(N);
    for (int i = 0; i < N; ++i)
    {
        double a = diagA(i);
        double b = diagB(i);
        // b should be > 0 here (powers of 10), but guard anyway:
        if (b == 0.0)
        {
            if (a == 0.0)
                eigs.emplace_back(std::numeric_limits<double>::quiet_NaN(), 0.0);
            else
                eigs.emplace_back(std::numeric_limits<double>::infinity(), 0.0);
        }
        else
        {
            eigs.emplace_back(a / b, 0.0);
        }
    }

    return {A, B, eigs};
}
