# QREV Voltammetry Simulation & Parameter Optimization Pipeline

> **This repository contains the complete computational and automation pipeline developed as part of the engineering thesis.** It connects every stage of the research workflow: from loading experimental cyclic voltammetry data to high-speed kinetic modeling and automated graphic exports.

contact: **abderraouf.haddad@ensh.dz**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20451254.svg)](https://doi.org/10.5281/zenodo.20451254)

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
Dimensionless Frequency Range: Covers $\omega \in [10^{-6}, 10^{2}]$ including intermediate fractional values ($10^{-1.5}$).
Execution Speed: Leverages Numba JIT (Just-In-Time) compilation to achieve native execution speeds (equivalent to Fortran).
---
## 📁 Project Structure
```plaintext
qrev-voltammetry/
│
├── qrev_find_alpha.py        # Main execution script (UI Dialog + Grid-Search)
├── mcmc_voltammetry.py       # Bayesian parameter estimation program
├── requirements.txt          # Python dependencies (numpy, numba, matplotlib)
├── .gitignore                # Prevents caching and output directory commits
└── README.md                 # Project documentation & structure

```
---
⚙️ Core Scripts & Functions
### `qrev_find_alpha.py` — Main Pipeline
An automated optimization environment that runs the main parameter-search loops:

| Section / Function | Objective | Operational Description |
| :--- | :--- | :--- |
| `@njit _qrev_sim()` | Core Numerical Solver | Fast Soluble-Insoluble cyclic voltammetry model translated from original Fortran code. |
| `@njit _qrev_rmse()` | Objective Function | Computes Root Mean Square Error, matching each experimental point to its nearest simulated potential. |
| `load_experimental_data()` | Data Parser | Cleans, filters, and loads space/tab-separated files containing `POT Y ` data. |
| `find_best_alpha()` | Parameter Optimization | Grid-search solver isolating minimum RMSE across all combinations of $\alpha$ and $\omega$. |
| `save_results()` | File Automation | Generates data tables, full curve files, and logs a **Top 10 Best Fits** summary table. |
| `plot_results()` | Data Visualization | Generates and exports the final thesis-ready figure comparing simulation and experimental data. |

Running the execution script:
```bash
# Opens a graphical file dialog to select data manually
python qrev_find_alpha.py

# Pass data path directly via CLI with instant plotting enabled
python qrev_find_alpha.py data.txt --plot

```
🛠️ Setup & Installation
1. Clone the repository
```bash
git clone https://github.com/Saymyname-ops/qrev-voltammetry.git
cd qrev-voltammetry

```
2. Install Python dependencies
```bash
pip install -r requirements.txt

```
3. Data Format Expected
The custom parser reads plain text or CSV data formatted into 3 aligned columns (`POT` = Potential, `Y` = Dimensionless Current):
```plaintext
# POT          Y         
0.0000e+00   0.0000e+00   
-1.0000e-02  -1.5430e-03   
...

```
📈 Generated Outputs & Artifacts
All data files and visuals are formatted automatically for direct copy-pasting or inclusion into your text layout:

| Exported File | Artifact Class | Description for Thesis Integration |
| :--- | :--- | :--- |
| `qrev_rmse_results.csv` | Data Table | Complete error matrix for all analyzed $\omega$ and $\alpha$ sets. |
| `qrev_best_fit.csv` | Target Match | Paired potentials comparing matching indices of experimental and simulated models. |
| `qrev_best_curve.csv` | Curve Profile | High-resolution simulated points of the calculated optimal voltammogram. |
| `qrev_best_fit_plot.png` | Validation Figure | Finished plot ($I_{adm}$ vs $E_{adm}$) with an automatically generated title containing optimized variables. |

 

## 🔬 Research & Thesis Context

* **Degree Framework:** Engineering Degree in Process Engineering (*Ingénieur d'État en Génie des Procédés*).
* **Thesis Title:** *“Numerical simulation and critical analysis of Tafel curves obtained from cyclic voltammetry for the optimised determination of the charge transfer coefficient.”*


### Purpose of Public Code Tracking:
This publication guarantees complete transparency, absolute reproducibility, and verifiable technical documentation for the computational modeling component of the defense. It serves as a secure, **timestamped academic record** to prevent any false positives or conflicts with plagiarism detection algorithms during manuscript submission.

 

## 📜 License & Intellectual Property

All rights reserved © 2026 — Raouf haddad  
Deposited at **ONDA** (National Office of Copyright and Related Rights), Algeria.

---

## ✍️ Citation

```
Haddad,A. R. (2026). Numerical simulation and critical analysis of Tafel curves obtained from cyclic voltammetry for the optimised determination of the charge transfer coefficient [Engineering Degree Thesis]. 
Constantine National Polytechnic School,
Department of Process Engineering.
```
