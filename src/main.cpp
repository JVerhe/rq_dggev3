#include "dggev3_wrapper.hpp"
#include "pencil_generator.hpp"
#include "matrix_utils.hpp"

#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <chrono>
#include <cmath>
#include <vector>
#include <complex>
#include <functional>

using namespace Eigen;
using namespace std;

// Define function pointer type for pencil generators:
using PencilGenerator = std::function<Pencil(int)>;

// Forward declaration
int runBenchmark(PencilGenerator pencilGenerate);

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        cerr << "Usage: ./main <pencil_type>\n"
             << " 1 = random singular\n"
             << " 2 = triangular singular\n"
             << " 3 = logspace singular\n"
             << " 4 = regular\n";
        return 1;
    }

    std::string arg = argv[1];
    PencilGenerator pencilGenerate;

    if (arg == "1")
        pencilGenerate = generateRandomSingularPencil;
    else if (arg == "2")
        pencilGenerate = generateSingularTriangularPencil;
    else if (arg == "3")
        pencilGenerate = generateLogspaceSingularPencil;
    else if (arg == "4")
        pencilGenerate = generateRegularPencil;
    else
    {
        cerr << "Unknown pencil type '" << arg << "'.\n";
        return 1;
    }

    return runBenchmark(pencilGenerate);
}

int runBenchmark(PencilGenerator pencilGenerate)
{
    namespace fs = std::filesystem;
    fs::create_directories("results");

    // Timestamped output file
    auto now = chrono::system_clock::now();
    time_t now_time = chrono::system_clock::to_time_t(now);
    tm *local_tm = localtime(&now_time);

    char timestamp[64];
    strftime(timestamp, sizeof(timestamp), "%Y%m%d_%H%M%S", local_tm);
    string filename = string("results/results_") + timestamp + ".txt";

    ofstream fout(filename);
    if (!fout.is_open())
    {
        cerr << "Failed to open output file: " << filename << "\n";
        return 1;
    }

    fout << "# Generalized Eigenvalue Comparison (QR vs RQ)\n";
    fout << "# Timestamp: " << timestamp << "\n";
    fout << "# Columns: N, QR_error, RQ_error, QR_time_ms, RQ_time_ms\n\n";
    fout << std::scientific << std::setprecision(6);

    cout << "Benchmarking dggev3_qr vs dggev3_rq ...\n";
    cout << "Writing to " << filename << "\n\n";

    // Loop over matrix sizes
    for (int N = 8; N <= 210; N += 6)
    {
        int trials = 1e5;
        if (N >= 16)
            trials = 1e4;
        if (N >= 64)
            trials = 500;
        if (N >= 128)
            trials = 50;
        if (N >= 256)
            trials = 18;
        if (N >= 1024)
            trials = 2;
        if (N >= 2048)
            trials = 1;

        double qr_err_sum = 0.0, rq_err_sum = 0.0;
        double qr_time_sum = 0.0, rq_time_sum = 0.0;

        for (int t = 0; t < trials; ++t)
        {
            Pencil pencil = pencilGenerate(N);

            // QR method
            auto t1 = chrono::high_resolution_clock::now();
            auto qr_res = dggev3_qr_wrapper(false, true, pencil.A, pencil.B);
            auto t2 = chrono::high_resolution_clock::now();
            double qr_time = chrono::duration<double, milli>(t2 - t1).count();

            double qr_err = eigErrorNorm(pencil.eigenvalues,
                                         qr_res.alphar,
                                         qr_res.alphai,
                                         qr_res.beta);

            // RQ method
            t1 = chrono::high_resolution_clock::now();
            auto rq_res = dggev3_rq_wrapper(false, true, pencil.A, pencil.B);
            t2 = chrono::high_resolution_clock::now();
            double rq_time = chrono::duration<double, milli>(t2 - t1).count();

            double rq_err = eigErrorNorm(pencil.eigenvalues,
                                         rq_res.alphar,
                                         rq_res.alphai,
                                         rq_res.beta);

            qr_err_sum += qr_err;
            rq_err_sum += rq_err;
            qr_time_sum += qr_time;
            rq_time_sum += rq_time;
        }

        double qr_err_avg = qr_err_sum / trials;
        double rq_err_avg = rq_err_sum / trials;
        double qr_time_avg = qr_time_sum / trials;
        double rq_time_avg = rq_time_sum / trials;

        fout << setw(5) << N << "  "
             << setw(14) << qr_err_avg << "  "
             << setw(14) << rq_err_avg << "  "
             << setw(14) << qr_time_avg << "  "
             << setw(14) << rq_time_avg << "\n";

        cout << "N=" << setw(4) << N
             << " | trials=" << setw(5) << trials
             << " | QR_err=" << setw(10) << qr_err_avg
             << " | RQ_err=" << setw(10) << rq_err_avg
             << " | QR_t=" << setw(8) << qr_time_avg << " ms"
             << " | RQ_t=" << setw(8) << rq_time_avg << " ms\n";
    }

    fout.close();
    cout << "\nFinished. Results saved to " << filename << "\n";
    return 0;
}
