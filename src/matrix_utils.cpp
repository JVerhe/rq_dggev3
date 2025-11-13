#include "matrix_utils.hpp"
#include <algorithm>
#include <cmath>

using namespace Eigen;

double eigen_error_norm(const std::vector<std::complex<double>> &exact,
                        const VectorXd &alphar,
                        const VectorXd &alphai,
                        const VectorXd &beta)
{
    int n = static_cast<int>(exact.size());
    std::vector<std::complex<double>> computed(n);
    for (int i = 0; i < n; ++i)
        computed[i] = std::complex<double>(alphar[i], alphai[i]) / beta[i];

    // Sort both sets by magnitude for stable pairing
    auto mag_cmp = [](auto &a, auto &b)
    { return std::abs(a) < std::abs(b); };
    std::sort(computed.begin(), computed.end(), mag_cmp);

    auto sorted_exact = exact;
    std::sort(sorted_exact.begin(), sorted_exact.end(), mag_cmp);

    double err = 0.0;
    for (int i = 0; i < n; ++i)
        err += std::norm(computed[i] - sorted_exact[i]);

    return std::sqrt(err);
}
