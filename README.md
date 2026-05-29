# QREV Voltammetry Simulation & Parameter Optimization Pipeline

> **This repository contains the complete computational and automation pipeline developed as part of the engineering thesis.** It connects every stage of the research workflow: from loading experimental cyclic voltammetry data to high-speed kinetic modeling and automated graphic exports.

For permission requests or academic collaboration, please contact: **your.email@gmail.com**

---

## 📊 Pipeline Workflow

```text
Experimental Data (.txt/.csv)
             │
             ▼
    qrev_find_alpha.py 
             │
             ├──► _qrev_sim()     (Numba-JIT Accelerated Fortran Engine)
             ├──► _qrev_rmse()    (Automated Data Fitting Calculation)
             │
             ▼
        output_dir/  ───► qrev_rmse_results.csv  (Complete Grid-Search Grid)
                     ───► qrev_best_fit.csv      (Aligned Simulation vs. Exp)
                     ───► qrev_best_curve.csv    (Full Optimized CV Curve)
                     ───► qrev_best_fit_plot.png (Validation Figure)
```
Key Numerical Capabilities & Settings:
Grid-Search Matrix: Tests $99 \times 11$ ($1,089$) total combinations per execution.
$\alpha$ Optimization Resolution: Scans the full range of $\alpha \in [0.01, 0.99]$ with a step of $0.01$.
Dimensionless Frequency Range: Covers $\omega \in [10^{-6}, 10^{3}]$ including intermediate fractional values ($10^{-1.5}$).
Execution Speed: Leverages Numba JIT (Just-In-Time) compilation to achieve native execution speeds (equivalent to Fortran).
📁 Project Structure
```plaintext
qrev-voltammetry/
│
├── qrev_find_alpha.py        # Main execution script (UI Dialog + Grid-Search)
├── mcmc_voltammetry.py       # Bayesian parameter estimation program
├── requirements.txt          # Python dependencies (numpy, numba, matplotlib)
├── .gitignore                # Prevents caching and output directory commits
└── README.md                 # Project documentation & structure

```
⚙️ Core Scripts & Functions
`qrev_find_alpha.py` — Main Pipeline
An automated optimization environment that runs the main parameter-search loops:
Section / Function	Objective	Operational Description
`@njit _qrev_sim()`	Core Numerical Solver	Fast Soluble-Insoluble cyclic voltammetry model translated from original Fortran code.
`@njit _qrev_rmse()`	Objective Function	Computes Root Mean Square Error, matching each experimental point to its nearest simulated potential.
`load_experimental_data()`	Data Parser	Cleans, filters, and loads space/tab-separated files containing `POT Y Q` data.
`find_best_alpha()`	Parameter Optimization	Grid-search solver isolating minimum RMSE across all combinations of $\alpha$ and $\omega$.
`save_results()`	File Automation	Generates data tables, full curve files, and logs a Top 10 Best Fits summary table.
`plot_results()`	Data Visualization	Generates and exports the final thesis-ready figure comparing simulation and experimental data.
Running the execution script:
```bash
# Opens a graphical file dialog to select data manually
python qrev_find_alpha.py

# Pass data path directly via CLI with instant plotting enabled
python qrev_find_alpha.py data.txt --plot

```

