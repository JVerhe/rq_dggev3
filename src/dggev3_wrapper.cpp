#include "dggev3_wrapper.hpp"
#include <stdexcept>

extern "C"
{
    void dggev3_qr_(const char *jobvl, const char *jobvr, const int *n,
                    double *A, const int *lda, double *B, const int *ldb,
                    double *alphar, double *alphai, double *beta,
                    double *VL, const int *ldvl, double *VR, const int *ldvr,
                    double *work, const int *lwork, int *info);

    void dggev3_rq_(const char *jobvl, const char *jobvr, const int *n,
                    double *A, const int *lda, double *B, const int *ldb,
                    double *alphar, double *alphai, double *beta,
                    double *VL, const int *ldvl, double *VR, const int *ldvr,
                    double *work, const int *lwork, int *info);
}

// Generic helper to call either Fortran function
static GeneralizedEigenResult dggev3_common(
    void (*fortran_func)(const char *, const char *, const int *, double *, const int *, double *, const int *,
                         double *, double *, double *, double *, const int *, double *, const int *,
                         double *, const int *, int *),
    bool computeVL,
    bool computeVR,
    Eigen::MatrixXd A,
    Eigen::MatrixXd B)
{
    if (A.rows() != A.cols() || B.rows() != B.cols() || A.rows() != B.rows())
        throw std::invalid_argument("A and B must be square and of same size");

    int n = static_cast<int>(A.rows());
    char jobvl = computeVL ? 'V' : 'N';
    char jobvr = computeVR ? 'V' : 'N';

    Eigen::VectorXd alphar(n), alphai(n), beta(n);
    Eigen::MatrixXd VL = Eigen::MatrixXd::Zero(computeVL ? n : 1, computeVL ? n : 1);
    Eigen::MatrixXd VR = Eigen::MatrixXd::Zero(computeVR ? n : 1, computeVR ? n : 1);

    int lda = n, ldb = n, ldvl = n, ldvr = n, info = 0;
    int lwork = -1;
    double work_query;

    // Workspace query
    fortran_func(&jobvl, &jobvr, &n,
                 A.data(), &lda, B.data(), &ldb,
                 alphar.data(), alphai.data(), beta.data(),
                 VL.data(), &ldvl, VR.data(), &ldvr,
                 &work_query, &lwork, &info);

    lwork = static_cast<int>(work_query);
    Eigen::VectorXd work(lwork);

    // Actual computation
    fortran_func(&jobvl, &jobvr, &n,
                 A.data(), &lda, B.data(), &ldb,
                 alphar.data(), alphai.data(), beta.data(),
                 VL.data(), &ldvl, VR.data(), &ldvr,
                 work.data(), &lwork, &info);

    return {alphar, alphai, beta, VL, VR, info};
}

GeneralizedEigenResult dggev3_qr_wrapper(
    bool computeVL,
    bool computeVR,
    Eigen::MatrixXd A,
    Eigen::MatrixXd B)
{
    return dggev3_common(dggev3_qr_, computeVL, computeVR, std::move(A), std::move(B));
}

GeneralizedEigenResult dggev3_rq_wrapper(
    bool computeVL,
    bool computeVR,
    Eigen::MatrixXd A,
    Eigen::MatrixXd B)
{
    return dggev3_common(dggev3_rq_, computeVL, computeVR, std::move(A), std::move(B));
}
