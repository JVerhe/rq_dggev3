#include "dggev3_wrapper.hpp"
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

    std::cout << "=== Using DGGEV3_QR ===\n";
    auto qr_result = dggev3_qr_wrapper(false, true, A, B);

    if (qr_result.info == 0)
    {
        for (int i = 0; i < A.rows(); ++i)
            std::cout << "λ" << i << " = (" << qr_result.alphar[i]
                      << " + " << qr_result.alphai[i]
                      << "i) / " << qr_result.beta[i] << "\n";
    }

    std::cout << "\n=== Using DGGEV3_RQ ===\n";
    auto rq_result = dggev3_rq_wrapper(false, true, A, B);

    if (rq_result.info == 0)
    {
        for (int i = 0; i < A.rows(); ++i)
            std::cout << "λ" << i << " = (" << rq_result.alphar[i]
                      << " + " << rq_result.alphai[i]
                      << "i) / " << rq_result.beta[i] << "\n";
    }
}
