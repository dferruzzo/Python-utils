"""Compare the fixed-step RK4 and RK45 solvers on the Lorenz system.

The script integrates the system with both methods, prints the resulting
trajectories, checks that their state arrays agree within a tolerance, and
displays phase-space plots. Running this module performs the full simulation
and opens Matplotlib windows.
"""
 
import numpy as np
from runge_kutta import rk4, rk45
import matplotlib.pyplot as plt


def lorenz(y_state, params):
    """Return the time derivatives for the Lorenz system.

    Parameters:
    y_state : array-like
        The current state [x, y, z].
    params : tuple
        The Lorenz parameters (sigma, rho, beta).

    Returns:
    numpy.ndarray
        The derivatives [dx/dt, dy/dt, dz/dt].
    """
    x, y, z = y_state
    sigma, rho, beta = params
    dxdt = sigma * (y - x)
    dydt = x * (rho - z) - y
    dzdt = x * y - beta * z
    return np.array([dxdt, dydt, dzdt])


# Configure the Lorenz system and its fixed-step integration interval.
params = (10.0, 28.0, 8.0 / 3.0)  # sigma, rho, beta
rhs = lambda t, y: lorenz(y, params)
y0 = np.array([1.0, 1.0, 1.0])
t0 = 0.0
t_end = 5000.0
h = 0.01
N = int((t_end - t0) / h)

# Allocate and initialize the RK4 trajectory.
t_values_rk4 = np.linspace(t0, t_end, N + 1)
y_values_rk4 = np.zeros((N + 1, len(y0)))
y_values_rk4[0] = y0
t_values_rk4[0] = t0

# Allocate and initialize the RK45 trajectory and error estimates.
t_values_rk45 = np.linspace(t0, t_end, N + 1)
y_values_rk45 = np.zeros((N + 1, len(y0)))
error_values_rk45 = np.zeros((N + 1, len(y0)))
y_values_rk45[0] = y0
t_values_rk45[0] = t0

# Integrate both methods over the same time grid.
for i in range(N):
    t_values_rk4[i + 1], y_values_rk4[i + 1] = rk4(rhs, t_values_rk4[i], y_values_rk4[i], h)
    t_values_rk45[i + 1], y_values_rk45[i + 1], error_values_rk45[i + 1] = rk45(rhs, t_values_rk45[i], y_values_rk45[i], h)

# Report the computed trajectories before comparing their states.
print("Resultados do método RK4:")
print("t_values_rk4:", t_values_rk4)
print("y_values_rk4:", y_values_rk4)
print("\nResultados do método RK45:")
print("t_values_rk45:", t_values_rk45)
print("y_values_rk45:", y_values_rk45)

# Verificação de consistência entre os métodos RK4 e RK45
tolerance = 1e-5
if np.allclose(y_values_rk4, y_values_rk45, atol=tolerance
):
    print("\nOs resultados dos métodos RK4 e RK45 são consistentes dentro da tolerância especificada.")
else:
    print("\nOs resultados dos métodos RK4 e RK45 não são consistentes dentro da tolerância especificada.")

# Plot the x-y and x-z phase-space projections for both methods.
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.plot(y_values_rk4[:, 0], y_values_rk4[:,1], label='RK4_x vs RK4_y', color='blue')
plt.plot(y_values_rk45[:, 0], y_values_rk45[:,1], label='RK45_x vs RK45_y', color='red', linestyle='dashed')
plt.legend()
plt.subplot(1, 2, 2)
plt.plot(y_values_rk4[:, 0], y_values_rk4[:,2], label='RK4_x vs RK4_z', color='blue')
plt.plot(y_values_rk45[:, 0], y_values_rk45[:,2], label='RK45_x vs RK45_y', color='red', linestyle='dashed')
plt.legend()

# Display the RK45 error estimates over time.
plt.figure(figsize=(12, 6))
plt.subplot(1, 3, 1)
plt.plot(t_values_rk45, error_values_rk45[:, 0], label='Error in x', color='blue')
plt.subplot(1, 3, 2)
plt.plot(t_values_rk45, error_values_rk45[:, 1], label='Error in y', color='green')
plt.subplot(1, 3, 3)
plt.plot(t_values_rk45, error_values_rk45[:, 2], label='Error in z', color='red')
plt.legend()
plt.show()

