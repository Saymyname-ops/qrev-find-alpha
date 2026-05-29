import numpy as np
from qrev_find_alpha import _qrev_sim

# Check simulation output scale
print("Simulation output scale check:")
sim_pot, sim_y, sim_q = _qrev_sim(0.0, -30.0, 10000, 0.5, 0.5, 0.01)
print(f"Pot range: [{np.min(sim_pot):.2f}, {np.max(sim_pot):.2f}]")
print(f"Y range: [{np.min(sim_y):.4e}, {np.max(sim_y):.4e}]")
print(f"Q range: [{np.min(sim_q):.4e}, {np.max(sim_q):.4e}]")

print("\nFirst 5 points:")
for i in range(5):
    print(f"  pot={sim_pot[i]:.2f}, y={sim_y[i]:.4e}, q={sim_q[i]:.4e}")

print("\nLast 5 points:")
for i in range(-5, 0):
    print(f"  pot={sim_pot[i]:.2f}, y={sim_y[i]:.4e}, q={sim_q[i]:.4e}")

# Now test with different omega values
print("\n" + "="*50)
print("Testing different omega values with alpha=0.5:")
for omega in [1e6, 1e5, 1e4, 1e3, 1e2, 1e1, 1e0]:
    sim_pot, sim_y, sim_q = _qrev_sim(0.0, -30.0, omega, 0.5, 0.5, 0.01)
    print(f"omega={omega:.0e}: Y_range=[{np.min(sim_y):.4e}, {np.max(sim_y):.4e}], points={len(sim_pot)}")