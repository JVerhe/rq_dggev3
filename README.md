# QZ Benchmark Tool

This project is a C++ benchmarking tool designed to compare the accuracy and performance (execution time) of the QR and RQ methods for solving Generalized Eigenvalue Problems using custom wrappers for `dggev3`.

## Prerequisites (Debian)

* The GNU C++ and Fortran compilers

`sudo apt install g++ gfortran`

* The Eigen c++ library for linear algebra from link

`sudo apt-get install libeigen3-dev`

* The BLAS and LLAPACK libraries

`sudo apt install libblas-dev liblapack-dev`

* The Python modules matplotlib, numpy for visualization

`python3 -m pip install numpy matplotlib`

## Compiling the project

Update the **CXXFLAGS** field in the Makefile so the include path points to where eigen3 is installed locally. For example:

`CXXFLAGS = -O3 -Iinclude -I/path/to/eigen3/`

To compile the project, from the main folder simply type:

`make`

## Running the project

The program uses a configuration file in the subdirectory `/config/` to read input arguments. 2 different tasks can be executed as of now.

1) **runBenchmark**: generates random pencils of varying sizes and computes the average errors for both the QR and RQ method

2) **errorVariance**: generates random pencils for one fixed size and measures how the error compares between runs.

**To run the program with a specific config file, execute the command:**

`./main config/<nameOfConfigFile>.txt`

### 1. runBenchmark

**Config file layout:**
```
# Core Settings
task = runBenchmark
pencil_type = 1

# Benchmark Loop Parameters
range_start = 8
range_end = 210
range_step = 6

# Trial Limits: [Matrix Size Threshold] [Number of Trials]
[trial_limits]
0 1e5
16 1e4
64 500
128 50
256 18
1024 2
2048 1
```

Lines that start with "#" are ignored by the program.

* `task`: determines the current task
* `pencil_type`: selects the type of generated matrix pencil
  * `1` Random Singular
  * `2` Triangular Singular
  * `3` Logspace Singular (ill-posed)
  * `4` Random Regular
* `range_start, range_end, range_step`: define the step conditions for the size of the tested pencils
* `[trial limits]`: determines the number of trials for each size of matrix
  * `[Matrix Size Threshold]`: minimum matrix size
  * `[Number of Trials]`: corresponding total number of trials

### 2. errorVariance
**Config file layout:**
```
# Core Settings
task = error_variance
pencil_type = 1

# Variance Parameters
size = 50
trials = 12
```

Lines that start with "#" are ignored by the program.

* `task`: determines the current task
* `pencil_type`: selects the type of generated matrix pencil
  * `1` Random Singular
  * `2` Triangular Singular
  * `3` Logspace Singular (ill-posed)
  * `4` Random Regular
* `size`: determines the fixed size of the pencil
* `trials`: determines the number of individual trials ran

## Visualizing Results

The result of the program will be written into a .txt file which contains the timestap of the beginning of execution. The output file is located in the directory `/results/`.

* **Benchmark output format**: `results_YYYYMMDD_HHMMSS.txt`
  * Columns: `N`, `QR_error_avg`, `RQ_error_avg`, `QR_time_avg_ms`, `RQ_time_avg_ms`
* **Variance output format**: `variance_YYYYMMDD_HHMMSS.txt`
  * Columns: `QR_error_avg`, `RQ_error_avg`, `QR_time_avg_ms`, `RQ_time_avg_ms`

The output of the runBenchmark task can be visualized using a Python script by issueing the following command from the main directory:

`python3 plots/make_convergence_plot.py`

## Deleting results

Deleting all results in `/results/` and plots in `/plots/` can be done via:

`make clean`