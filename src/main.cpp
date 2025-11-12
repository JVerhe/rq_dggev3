#include "dggev3_qr_wrapper.hpp"
#include <iostream>

int main()
{
    Eigen::Matrix3d A;
    Eigen::Matrix3d B;

    A << 1.0, 1.0, 5.0,
        1.0, 2.0, 6.0,
        1.0, 1.0, 3.0;

    B << 2.0, 4.0, 3.0,
        1.0, 1.0, 2.0,
        1.0, 1.0, 1.0;

    auto result = dggev3_qr_wrapper(true, true, A, B);

    if (result.info == 0)
    {
        std::cout << "Generalized eigenvalues:\n";
        for (int i = 0; i < A.rows(); ++i)
        {
            std::cout << "λ" << i << " = ("
                      << result.alphar[i] << " + "
                      << result.alphai[i] << "i) / "
                      << result.beta[i] << "\n";
        }

        if (result.VR.size() > 1)
        {
            std::cout << "\nRight eigenvectors (VR):\n"
                      << result.VR << "\n";
        }
    }
    else
    {
        std::cerr << "LAPACK dggev3 failed with info = " << result.info << "\n";
    }

    return 0;
}
