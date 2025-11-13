#pragma once
#include <Eigen/Dense>
#include <complex>
#include <vector>

struct Pencil
{
    Eigen::MatrixXd A;
    Eigen::MatrixXd B;
    std::vector<std::complex<double>> eigenvalues;
};

// Generate a well-conditioned regular pencil
// A = P * D * Q,  B = P * Q
// where P and Q are orthogonal (from QR of random matrices)
// and D = diag(1, 2, ..., N)
Pencil generate_regular_pencil(int N);

Pencil generate_singular_pencil(int N);

// (Optional extensions for future)
// Pencil generate_singular_pencil(int N);
// Pencil generate_defective_pencil(int N);
// Pencil generate_complex_pencil(int N);
