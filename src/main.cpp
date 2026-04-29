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
#include <string>
#include <sstream>
#include <algorithm>
#include <map>

using namespace Eigen;
using namespace std;

using PencilGenerator = std::function<Pencil(int)>;

class ConfigParser {
private:
    static string trim(const string& s) {
        auto start = s.find_first_not_of(" \t\r\n");
        if (start == string::npos) return "";
        auto end = s.find_last_not_of(" \t\r\n");
        return s.substr(start, end - start + 1);
    }

public:
    map<string, string> values;
    map<string, vector<vector<string>>> blocks;

    bool load(const string& filename) {
        ifstream fin(filename);
        if (!fin.is_open()) return false;

        string line, current_section = "";
        while (getline(fin, line)) {
            line = trim(line);
            if (line.empty() || line[0] == '#' || line[0] == ';') continue; 

            if (line.front() == '[' && line.back() == ']') {
                current_section = line.substr(1, line.size() - 2);
                continue;
            }

            if (current_section.empty()) {
                auto pos = line.find('=');
                if (pos != string::npos) {
                    string key = trim(line.substr(0, pos));
                    string val = trim(line.substr(pos + 1));
                    values[key] = val;
                }
            } else {
                istringstream iss(line);
                vector<string> row;
                string token;
                while (iss >> token) {
                    row.push_back(token);
                }
                if (!row.empty()) {
                    blocks[current_section].push_back(row);
                }
            }
        }
        return true;
    }

    string getString(const string& key, const string& default_value = "") const {
        auto it = values.find(key);
        return (it != values.end()) ? it->second : default_value;
    }

    template<typename T>
    T get(const string& key, const T& default_value = T()) const {
        auto it = values.find(key);
        if (it != values.end()) {
            istringstream iss(it->second);
            T result;
            iss >> result;
            return result;
        }
        return default_value;
    }
};

int runBenchmark(PencilGenerator pencilGenerate, const ConfigParser& cfg);
int errorVariance(PencilGenerator pencilGenerate, const ConfigParser& cfg);

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        cerr << "Usage: ./main <config_file.txt>\n";
        return 1;
    }

    ConfigParser cfg;
    if (!cfg.load(argv[1])) {
        cerr << "Failed to parse configuration file: " << argv[1] << "\n";
        return 1;
    }

    string task_name = cfg.getString("task");
    string pencil_type = cfg.getString("pencil_type");

    PencilGenerator pencilGenerate;

    if (pencil_type == "1")
        pencilGenerate = generateRandomSingularPencil;
    else if (pencil_type == "2")
        pencilGenerate = generateSingularTriangularPencil;
    else if (pencil_type == "3")
        pencilGenerate = generateLogspaceSingularPencil;
    else if (pencil_type == "4")
        pencilGenerate = generateRegularPencil;
    else
    {
        cerr << "Unknown pencil type '" << pencil_type << "'.\n";
        return 1;
    }

    // Route task based on config file
    if (task_name == "runBenchmark") {
        return runBenchmark(pencilGenerate, cfg);
    } else if (task_name == "error_variance") {
        return errorVariance(pencilGenerate, cfg);
    } else {
        cerr << "Unknown task '" << task_name << "'.\n";
        return 1;
    }
}


int errorVariance(PencilGenerator pencilGenerate, const ConfigParser& cfg)
{
    namespace fs = std::filesystem;
    fs::create_directories("results");

    int size = cfg.get<int>("size", 50);
    int trials = cfg.get<int>("trials", 12);

    auto now = chrono::system_clock::now();
    time_t now_time = chrono::system_clock::to_time_t(now);
    tm *local_tm = localtime(&now_time);

    char timestamp[64];
    strftime(timestamp, sizeof(timestamp), "%Y%m%d_%H%M%S", local_tm);
    string filename = string("results/variance_") + timestamp + ".txt";

    ofstream fout(filename);
    if (!fout.is_open())
    {
        cerr << "Failed to open output file: " << filename << "\n";
        return 1;
    }

    fout << "# Generalized Eigenvalue Comparison (QR vs RQ)\n";
    fout << "# Timestamp: " << timestamp << "\n";
    fout << "# Task: " << cfg.getString("task") 
         << " | Pencil Type: " << cfg.getString("pencil_type") 
         << " | size = " << size << "\n";
    fout << "# Columns: QR_error, RQ_error, QR_time_ms, RQ_time_ms\n\n";
    
    fout << scientific << setprecision(6);

    cout << "Running error variance analysis ...\n";
    cout << "Matrix Size: " << size << " | Trials: " << trials << "\n";
    cout << "Writing to " << filename << "\n\n";

    for (int t = 0; t < trials; ++t)
    {
        Pencil pencil = pencilGenerate(size);

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

        // Log exact trial result to file
        fout << qr_err << "    "
             << rq_err << "    "
             << qr_time << "    "
             << rq_time << "\n";

        // print progress to console for long runs
        if ((t + 1) % 10 == 0 || t == trials - 1) {
            cout << "Completed " << (t + 1) << " / " << trials << " trials...\r" << flush;
        }
    }

    fout.close();
    cout << "\nFinished. Results saved to " << filename << "\n";
    return 0;
}

int runBenchmark(PencilGenerator pencilGenerate, const ConfigParser& cfg)
{
    namespace fs = std::filesystem;
    fs::create_directories("results");

    int start_N = cfg.get<int>("range_start", 8);
    int end_N   = cfg.get<int>("range_end", 210);
    int step_N  = cfg.get<int>("range_step", 6);

    vector<pair<int, int>> trial_limits;
    if (cfg.blocks.count("trial_limits")) {
        for (const auto& row : cfg.blocks.at("trial_limits")) {
            if (row.size() >= 2) {
                int threshold = stoi(row[0]);
                int trials = static_cast<int>(stod(row[1]));
                trial_limits.push_back({threshold, trials});
            }
        }
        sort(trial_limits.begin(), trial_limits.end());
    }

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
    fout << "# Task: " << cfg.getString("task") << " | Pencil Type: " << cfg.getString("pencil_type") << "\n";
    fout << "# Columns: N, QR_error, RQ_error, QR_time_ms, RQ_time_ms\n\n";
    fout << scientific << setprecision(6);

    cout << "Benchmarking dggev3_qr vs dggev3_rq ...\n";
    cout << "Writing to " << filename << "\n\n";

    for (int N = start_N; N <= end_N; N += step_N)
    {
        int trials = 1; 
        for (const auto& limit : trial_limits) {
            if (N >= limit.first) {
                trials = limit.second;
            }
        }

        double qr_err_sum = 0.0, rq_err_sum = 0.0;
        double qr_time_sum = 0.0, rq_time_sum = 0.0;

        for (int t = 0; t < trials; ++t)
        {
            Pencil pencil = pencilGenerate(N);

            auto t1 = chrono::high_resolution_clock::now();
            auto qr_res = dggev3_qr_wrapper(false, true, pencil.A, pencil.B);
            auto t2 = chrono::high_resolution_clock::now();
            double qr_time = chrono::duration<double, milli>(t2 - t1).count();

            double qr_err = eigErrorNorm(pencil.eigenvalues,
                                         qr_res.alphar,
                                         qr_res.alphai,
                                         qr_res.beta);

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