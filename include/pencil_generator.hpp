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
// where P and Q are orthogonal
// and D = diag(1, 2, ..., N)
Pencil generate_regular_pencil(int N);

// Generates a singular pencil
// where A, B are upper-triangular matrices
// Some random elements on the diagonals have been set to 0
// Resulting in fake eigs for (aii, bii == 0)
// And infinite eigs (only bii == 0)
Pencil generate_singular_triangular_pencil(int N);

Pencil generate_singular_pencil(int N, bool infEigvals = false);

Pencil generate_illconditioned_B_pencil(int N,
                                        bool use_integer_DA = true);
