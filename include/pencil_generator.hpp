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
Pencil generateRegularPencil(int N);

// Generates a singular pencil
// where A, B are upper-triangular matrices
// Some random elements on the diagonals have been set to 0
// Resulting in fake eigs for (aii, bii == 0)
// And infinite eigs (only bii == 0)
Pencil generateSingularTriangularPencil(int N);

Pencil generateRandomSingularPencil(int N);

Pencil generateLogspaceSingularPencil(int N);
