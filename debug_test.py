import numpy as np
from qrev_find_alpha import _qrev_sim

# Test data from known parameters (alpha=0.5, omega=10000)
print("Testing with alpha=0.5, omega=10000")
sim_pot, sim_y, sim_q = _qrev_sim(0.0, -30.0, 10000, 0.5, 0.5, 0.01)
print(f"Simulation: {len(sim_pot)} points")
print(f"Peak Y: {np.min(sim_y):.6f} at pot={sim_pot[np.argmin(sim_y)]:.2f}")

# Create "experimental" data with noise (should recover alpha=0.5)
np.random.seed(42)
idx = np.arange(0, len(sim_pot), max(1, len(sim_pot)//50))  # ~50 points
noise_y = np.random.normal(0, 0.01, len(idx))
noise_q = np.random.normal(0, 0.01, len(idx))

exp_pot = sim_pot[idx]
exp_y = sim_y[idx] + noise_y
exp_q = sim_q[idx] + noise_q

print(f"Experimental: {len(exp_pot)} points")

# Now test different alphas with omega fixed at 10000
print("\nTesting different alphas with omega=10000:")
best_rmse = float('inf')
best_alpha = None

for alpha_test in [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]:
    beta_test = 1.0 - alpha_test
    try:
        sim_pot_t, sim_y_t, sim_q_t = _qrev_sim(0.0, -30.0, 10000, alpha_test, beta_test, 0.01)
        
        # Simple RMSE (match by index for now)
        min_len = min(len(sim_y_t), len(exp_y))
        if min_len > 10:
            diff = sim_y_t[:min_len] - exp_y[:min_len]
            rmse = np.sqrt(np.mean(diff**2))
            print(f"  alpha={alpha_test:.2f}: RMSE={rmse:.6f}")
            if rmse < best_rmse:
                best_rmse = rmse
                best_alpha = alpha_test
    except Exception as e:
        print(f"  alpha={alpha_test:.2f}: ERROR - {e}")

print(f"\nBest alpha found: {best_alpha} (RMSE={best_rmse:.6f})")
print(f"Expected alpha: 0.5")