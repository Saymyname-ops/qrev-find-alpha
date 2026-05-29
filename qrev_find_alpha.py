"""
QREV - Python implementation of Fortran electrochemistry simulation
SOLUBLE-INSOLUBLE mechanism with cyclic voltammetry

Finds the best alpha (and omega) by comparing simulation with experimental data.

Usage:
    python qrev_find_alpha.py                     # Opens file dialog to choose data
    python qrev_find_alpha.py data.txt             # Pass data file directly
    python qrev_find_alpha.py data.txt --plot       # Also show plots

Experimental data format: Tab/space-separated file with columns POT Y Q
"""

import numpy as np
import numba
from numba import njit
import sys
import csv
import time
from pathlib import Path


# ---------------------------------------------------------------------------
# Numba-JIT core simulation (Fortran-speed)
# ---------------------------------------------------------------------------

@njit(cache=True)
def _qrev_sim(init, limit, omega, alpha, beta, delta):
    """
    Core QREV simulation compiled with Numba.
    Returns (pot_arr, y_arr, q_arr) trimmed to actual computed length.
    Stops backward scan when q returns to >= 0 (oxidation complete).
    """
    it = int(round((init - limit) / delta))
    int1 = it
    int2 = 7 * int1

    f1 = 2.0 * (limit - init)
    c_val = 1.0 / (2.0 * np.sqrt(delta))
    sqrt_pi = np.sqrt(np.pi)
    alpha_beta = alpha + beta

    d = np.empty(int2 + 1)
    for i in range(int2 + 1):
        d[i] = np.sqrt(float(i + 1)) - np.sqrt(float(i))

    pot = np.empty(int2 + 1)
    y_arr = np.zeros(int2 + 1)
    q_arr = np.zeros(int2 + 1)
    x_arr = np.zeros(int2 + 1)

    pot[0] = init
    y_arr[0] = 0.0
    q_arr[0] = 0.0
    x_arr[0] = 0.0

    actual_len = int2 + 1

    for k in range(1, int2 + 1):
        if k <= int1:
            pot[k] = init - delta * k
            exp_factor = np.exp(-delta * k)
        else:
            pot[k] = 2.0 * limit - init + delta * (k - int1)
            exp_factor = np.exp(delta * (k - int1) + f1)

        a = c_val * (exp_factor ** alpha_beta)
        b = (c_val / omega) * (exp_factor ** alpha)

        sigma = 0.0
        for j in range(1, k):
            sigma += d[k - j] * x_arr[j]

        x_arr[k] = (c_val - a - sigma) / (1.0 + b)
        y_arr[k] = -x_arr[k] * sqrt_pi
        q_arr[k] = q_arr[k - 1] + y_arr[k] * delta

        if k > int1 and q_arr[k] >= 0.0:
            actual_len = k + 1
            break

        if abs(y_arr[k]) > 1e30:
            actual_len = k + 1
            break

    return pot[:actual_len], y_arr[:actual_len], q_arr[:actual_len]


@njit(cache=True)
def _qrev_rmse(sim_pot, sim_y, exp_pot, exp_y):
    """
    Calculate RMSE between simulation and experimental Y,
    matching each experimental point to the nearest simulation potential.
    """
    n_exp = len(exp_pot)
    n_sim = len(sim_pot)
    total_sq = 0.0

    for i in range(n_exp):
        best_idx = 0
        best_dist = abs(sim_pot[0] - exp_pot[i])
        for j in range(1, n_sim):
            dist = abs(sim_pot[j] - exp_pot[i])
            if dist < best_dist:
                best_dist = dist
                best_idx = j
        diff = sim_y[best_idx] - exp_y[i]
        total_sq += diff * diff

    return np.sqrt(total_sq / float(n_exp))


# ---------------------------------------------------------------------------
# Warm up JIT (runs once on import)
# ---------------------------------------------------------------------------

def _warmup_jit():
    _ = _qrev_sim(0.0, -1.0, 1.0, 0.5, 0.5, 0.1)
    _ = _qrev_rmse(
        np.array([0.0, -0.5]), np.array([1.0, 2.0]),
        np.array([0.0, -0.5]), np.array([1.1, 2.1])
    )

print("Warming up Numba JIT compiler...", end=' ', flush=True)
_warmup_jit()
print("done.")


# ---------------------------------------------------------------------------
# Experimental data loader
# ---------------------------------------------------------------------------

def load_experimental_data(filepath):
    """Load experimental data from tab/space-separated file (POT Y Q columns)."""
    data = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('"'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    pot = float(parts[0])
                    y = float(parts[1])
                    q = float(parts[2]) if len(parts) >= 3 else 0.0
                    data.append({'pot': pot, 'y': y, 'q': q})
                except ValueError:
                    continue
    return data


# ---------------------------------------------------------------------------
# Find best alpha/omega
# ---------------------------------------------------------------------------

def find_best_alpha(experimental_data, alphas, omegas, init=0.0, limit=-30.0, delta=0.01):
    """Search over all (omega, alpha) combinations and find the minimum RMSE."""

    exp_pot = np.array([d['pot'] for d in experimental_data], dtype=np.float64)
    exp_y = np.array([d['y'] for d in experimental_data], dtype=np.float64)

    it = int(round((init - limit) / delta))
    int2 = 7 * it

    total = len(omegas) * len(alphas)
    print(f"Experimental data points: {len(exp_pot)}")
    print(f"Simulation steps per run: {int2}")
    print(f"Total combinations to test: {total}")

    best = {'omega': None, 'alpha': None, 'beta': None, 'rmse': float('inf')}
    all_results = []
    count = 0
    t0 = time.time()

    for omega_idx, omega in enumerate(omegas):
        t_omega = time.time()
        print(f"\n  Omega={omega:.2e} [{omega_idx+1}/{len(omegas)}]...", end='', flush=True)

        omega_best_rmse = float('inf')

        for alpha in alphas:
            beta = 1.0 - alpha
            count += 1

            try:
                sim_pot, sim_y, sim_q = _qrev_sim(init, limit, omega, alpha, beta, delta)
                rmse = _qrev_rmse(sim_pot, sim_y, exp_pot, exp_y)
            except Exception as e:
                rmse = float('inf')

            result = {'omega': omega, 'alpha': alpha, 'beta': beta, 'rmse': rmse}
            all_results.append(result)

            if rmse < omega_best_rmse:
                omega_best_rmse = rmse

            if rmse < best['rmse']:
                best = result.copy()
                print(f"\n    ** New best: alpha={alpha:.2f}, omega={omega:.2e}, RMSE={rmse:.6e}")

        elapsed = time.time() - t_omega
        print(f" done in {elapsed:.1f}s (best RMSE this omega: {omega_best_rmse:.6e})")

        if omega_idx == 0:
            est_total = elapsed * len(omegas)
            print(f"    Estimated total time: {est_total:.0f}s ({est_total/60:.1f}min)")

    total_time = time.time() - t0
    print(f"\nSearch completed in {total_time:.1f}s ({total_time/60:.1f}min)")

    return all_results, best


# ---------------------------------------------------------------------------
# Save results
# ---------------------------------------------------------------------------

def save_results(all_results, best, experimental_data, init, limit, delta, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    # 1) RMSE table
    rmse_file = output_dir / "qrev_rmse_results.csv"
    with open(rmse_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['omega', 'alpha', 'beta', 'rmse'])
        writer.writeheader()
        writer.writerows(all_results)
    print(f"RMSE results saved: {rmse_file}")

    # 2) Best-fit simulation data
    sim_pot, sim_y, sim_q = _qrev_sim(init, limit, best['omega'], best['alpha'], best['beta'], delta)
    best_file = output_dir / "qrev_best_fit.csv"

    exp_pot = np.array([d['pot'] for d in experimental_data], dtype=np.float64)
    exp_y = np.array([d['y'] for d in experimental_data], dtype=np.float64)
    exp_q = np.array([d['q'] for d in experimental_data], dtype=np.float64)

    with open(best_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['pot_sim', 'y_sim', 'q_sim', 'pot_exp', 'y_exp', 'q_exp'])
        for i in range(len(exp_pot)):
            idx = np.argmin(np.abs(sim_pot - exp_pot[i]))
            writer.writerow([
                f"{sim_pot[idx]:.8e}", f"{sim_y[idx]:.8e}", f"{sim_q[idx]:.8e}",
                f"{exp_pot[i]:.8e}", f"{exp_y[i]:.8e}", f"{exp_q[i]:.8e}"
            ])
    print(f"Best fit data saved: {best_file}")

    # 3) Full simulation curve for best fit
    curve_file = output_dir / "qrev_best_curve.csv"
    with open(curve_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['pot', 'y', 'q'])
        for i in range(len(sim_pot)):
            writer.writerow([f"{sim_pot[i]:.8e}", f"{sim_y[i]:.8e}", f"{sim_q[i]:.8e}"])
    print(f"Best fit curve saved: {curve_file}")

    # 4) Top-10 results
    sorted_results = sorted(all_results, key=lambda r: r['rmse'])
    print("\n" + "=" * 65)
    print("TOP 10 BEST FITS")
    print("=" * 65)
    print(f"{'Rank':<5} {'Alpha':<8} {'Beta':<8} {'Omega':<14} {'RMSE':<14}")
    print("-" * 65)
    for i, r in enumerate(sorted_results[:10]):
        print(f"{i+1:<5} {r['alpha']:<8.3f} {r['beta']:<8.3f} {r['omega']:<14.4e} {r['rmse']:<14.6e}")

    return rmse_file, best_file


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_results(best, experimental_data, init, limit, delta, output_dir):
    try:
        import matplotlib
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed, skipping plots.")
        return

    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    sim_pot, sim_y, sim_q = _qrev_sim(init, limit, best['omega'], best['alpha'], best['beta'], delta)

    exp_pot = [d['pot'] for d in experimental_data]
    exp_y = [d['y'] for d in experimental_data]
    exp_q = [d['q'] for d in experimental_data]

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    ax.plot(sim_pot, sim_y, 'b-', label='Simulation', linewidth=1.5)
    ax.plot(exp_pot, exp_y, 'ro', label='Experimental', markersize=3)
    ax.set_xlabel('E adm')
    ax.set_ylabel('I adm')
    ax.set_title(f'alpha={best["alpha"]:.2f}, omega={best["omega"]:.2e}, RMSE={best["rmse"]:.4e}')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plot_file = output_dir / "qrev_best_fit_plot.png"
    plt.savefig(plot_file, dpi=150)
    print(f"Plot saved: {plot_file}")
    plt.show()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    INIT = 0.0
    LIMIT = -30.0
    DELTA = 0.01

    alphas = [round(i * 0.01, 4) for i in range(1, 100)]
    # Focus on sensitive range where alpha matters most: 1e-3 to 1e3
    # Based on debugging: alpha sensitivity is LOW at high omega (>1e4) 
    # and HIGH at low omega (<1e-2)
    # Include 10^-1.5 as specifically requested
    # Also include lower omega values 10^-4, 10^-5, 10^-6 as requested
    omegas = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 10**-1.5, 1e-1, 1e0, 1e1, 1e2, 1e3]  # 11 values

    output_dir = Path(r"C:\Users\ranin\qrev_output")
    output_dir.mkdir(exist_ok=True)

    show_plot = '--plot' in sys.argv

    # --- Load experimental data ---
    if len(sys.argv) > 1 and not sys.argv[1].startswith('--'):
        exp_file = sys.argv[1]
    else:
        try:
            import tkinter as tk
            from tkinter import filedialog
            root = tk.Tk()
            root.withdraw()
            root.attributes('-topmost', True)
            exp_file = filedialog.askopenfilename(
                title="Choose Experimental Data File (Tab-separated: POT Y Q)",
                filetypes=[("Text files", "*.txt"), ("CSV files", "*.csv"), ("All files", "*.*")],
                initialdir=r"C:\Users\ranin"
            )
            root.destroy()
        except Exception:
            print("ERROR: Could not open file dialog. Pass the file path as argument:")
            print("  python qrev_find_alpha.py <data_file.txt>")
            return

    if not exp_file:
        print("ERROR: No data file selected.")
        return

    print("=" * 65)
    print("QREV - FIND BEST ALPHA")
    print("=" * 65)
    print(f"Data file:      {exp_file}")
    print(f"INIT =          {INIT}")
    print(f"LIMIT =         {LIMIT}")
    print(f"DELTA =         {DELTA}")
    print(f"Alpha range:    0.01 to 0.99 ({len(alphas)} values)")
    print(f"Omega values:   {len(omegas)}")
    print(f"Combinations:   {len(alphas) * len(omegas)}")
    print("=" * 65)

    experimental_data = load_experimental_data(exp_file)
    if not experimental_data:
        print("ERROR: No data loaded from file.")
        return

    print(f"Loaded {len(experimental_data)} data points")
    print(f"  POT range: {experimental_data[0]['pot']:.2f} to {experimental_data[-1]['pot']:.2f}")
    print(f"  Y range:   {min(d['y'] for d in experimental_data):.4e} to "
          f"{max(d['y'] for d in experimental_data):.4e}")

    # --- Search ---
    all_results, best = find_best_alpha(
        experimental_data, alphas, omegas, INIT, LIMIT, DELTA
    )

    # --- Results ---
    print("\n" + "=" * 65)
    print("BEST FIT FOUND")
    print("=" * 65)
    print(f"  Alpha:  {best['alpha']}")
    print(f"  Beta:   {best['beta']}")
    print(f"  Omega:  {best['omega']}")
    print(f"  RMSE:   {best['rmse']:.6e}")
    print("=" * 65)

    # --- Save ---
    save_results(all_results, best, experimental_data, INIT, LIMIT, DELTA, output_dir)

    # --- Plot ---
    plot_results(best, experimental_data, INIT, LIMIT, DELTA, output_dir)

    try:
        input("\nPress ENTER to exit...")
    except EOFError:
        # Handle case when running without stdin (e.g., from script)
        pass


if __name__ == "__main__":
    main()
