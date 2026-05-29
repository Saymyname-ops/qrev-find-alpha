import numpy as np
from qrev_find_alpha import _qrev_sim, _qrev_rmse

print("Testing sensitivity at LOW vs HIGH omega")
print("="*50)

# Generate "true" simulation with alpha=0.5, omega=0.01 (LOW omega - should be sensitive)
print("\nLOW omega case (omega=0.01):")
sim_pot_low, sim_y_low, sim_q_low = _qrev_sim(0.0, -30.0, 0.01, 0.5, 0.5, 0.01)
print(f"  Simulation points: {len(sim_pot_low)}")
print(f"  Y range: [{np.min(sim_y_low):.4e}, {np.max(sim_y_low):.4e}]")

# Create noisy experimental copy
np.random.seed(42)
idx = np.arange(0, len(sim_pot_low), 100)
exp_pot = sim_pot_low[idx].copy()
exp_y = sim_y_low[idx].copy() + np.random.normal(0, 0.005, len(idx))

print(f"  Experimental points: {len(exp_pot)}")
print(f"  Y range: [{np.min(exp_y):.4e}, {np.max(exp_y):.4e}]")

# Test different alphas at this LOW omega
print("\nRMSE vs alpha (omega=0.01):")
best_rmse = float('inf')
best_alpha = None
for alpha_test in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
    beta_test = 1.0 - alpha_test
    sim_pot_t, sim_y_t, sim_q_t = _qrev_sim(0.0, -30.0, 0.01, alpha_test, beta_test, 0.01)
    rmse_t = _qrev_rmse(sim_pot_t, sim_y_t, exp_pot, exp_y)
    print(f"  alpha={alpha_test:.1f}: RMSE={rmse_t:.6f}")
    if rmse_t < best_rmse:
        best_rmse = rmse_t
        best_alpha = alpha_test
print(f"  -> Best alpha: {best_alpha} (RMSE={best_rmse:.6f})")

print("\n" + "="*50)
print("HIGH omega case (omega=10000):")
# Now test with HIGH omega
sim_pot_high, sim_y_high, sim_q_high = _qrev_sim(0.0, -30.0, 10000, 0.5, 0.5, 0.01)
print(f"  Simulation points: {len(sim_pot_high)}")
print(f"  Y range: [{np.min(sim_y_high):.4e}, {np.max(sim_y_high):.4e}]")

# Same experimental data (from LOW omega sim)
print(f"  Experimental points: {len(exp_pot)} (from LOW omega simulation)")

# Test different alphas at HIGH omega
print("\nRMSE vs alpha (omega=10000):")
best_rmse_high = float('inf')
best_alpha_high = None
for alpha_test in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
    beta_test = 1.0 - alpha_test
    sim_pot_t, sim_y_t, sim_q_t = _qrev_sim(0.0, -30.0, 10000, alpha_test, beta_test, 0.01)
    rmse_t = _qrev_rmse(sim_pot_t, sim_y_t, exp_pot, exp_y)  # Compare to LOW omega exp data
    print(f"  alpha={alpha_test:.1f}: RMSE={rmse_t:.6f}")
    if rmse_t < best_rmse_high:
        best_rmse_high = rmse_t
        best_alpha_high = alpha_test
print(f"  -> Best alpha: {best_alpha_high} (RMSE={best_rmse_high:.6f})")

print("\n" + "="*50)
print("CONCLUSION:")
print("- At LOW omega (0.01): alpha sensitivity is HIGH (different alpha -> different shapes)")  
print("- At HIGH omega (10000): alpha sensitivity is LOW (all alpha give similar shapes)")
print("- Your experimental data likely came from a LOW omega process")
print("- So you need to search LOW omega values (1e-3 to 1e0) to find correct alpha")