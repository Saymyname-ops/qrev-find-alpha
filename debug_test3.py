import numpy as np
from qrev_find_alpha import _qrev_sim, _qrev_rmse

# Test the RMSE function directly
print("Testing RMSE function:")

# Generate simulation data
sim_pot, sim_y, sim_q = _qrev_sim(0.0, -30.0, 10000, 0.5, 0.5, 0.01)
print(f"Simulation: {len(sim_pot)} points")

# Create experimental data as a NOISY copy of simulation (should give low RMSE)
np.random.seed(42)
# Sample every 200 points to reduce data size
idx = np.arange(0, len(sim_pot), 200)
exp_pot = sim_pot[idx].copy()
exp_y = sim_y[idx].copy() + np.random.normal(0, 0.01, len(idx))  # 1% noise
exp_q = sim_q[idx].copy() + np.random.normal(0, 0.01, len(idx))

print(f"Experimental (noisy copy): {len(exp_pot)} points")
print(f"Y range: exp=[{np.min(exp_y):.4e}, {np.max(exp_y):.4e}], sim=[{np.min(sim_y):.4e}, {np.max(sim_y):.4e}]")

# Calculate RMSE
rmse = _qrev_rmse(sim_pot, sim_y, exp_pot, exp_y)
print(f"RMSE (should be small): {rmse:.6f}")

# Now test with shifted experimental data (add constant offset)
exp_y_shifted = exp_y + 0.5  # Add 0.5 to all Y values
rmse_shifted = _qrev_rmse(sim_pot, sim_y, exp_pot, exp_y_shifted)
print(f"RMSE with +0.5 offset: {rmse_shifted:.6f}")

# Now test with scaled experimental data (multiply by 2)
exp_y_scaled = exp_y * 2.0
rmse_scaled = _qrev_rmse(sim_pot, sim_y, exp_pot, exp_y_scaled)
print(f"RMSE with *2.0 scale: {rmse_scaled:.6f}")

# Test different alphas with fixed omega=10000
print("\n" + "="*50)
print("Testing RMSE vs alpha (omega=10000 fixed, noisy exp data):")
for alpha_test in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
    beta_test = 1.0 - alpha_test
    try:
        sim_pot_t, sim_y_t, sim_q_t = _qrev_sim(0.0, -30.0, 10000, alpha_test, beta_test, 0.01)
        rmse_t = _qrev_rmse(sim_pot_t, sim_y_t, exp_pot, exp_y)
        print(f"  alpha={alpha_test:.1f}: RMSE={rmse_t:.6f}")
    except Exception as e:
        print(f"  alpha={alpha_test:.1f}: ERROR")