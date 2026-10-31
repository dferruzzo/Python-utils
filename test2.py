import numpy as np
from runge_kutta import rk4, rk45
import matplotlib.pyplot as plt

def second_ordem_system(y_state, params):
    """Return the time derivatives for a second-order system.

    Parameters:
    y_state : array-like
        The current state [x, v].
    params : tuple
        The system parameters (k, m, c).

    Returns:
    numpy.ndarray
        The derivatives [dx/dt, dv/dt].
    """
    x, v = y_state
    k, m, c = params
    dxdt = v
    dvdt = -(k/m) * x - (c/m) * v
    return np.array([dxdt, dvdt])

# Configure the second-order system and its fixed-step integration interval.
params = (1.0, 1.0, 0.1)  #
k, m, c = params
rhs = lambda t, y: second_ordem_system(y, params)
y0 = np.array([1.0, 0.0])  # Initial conditions
t0 = 0.0
t_end = 20.0
h = 0.001
N = int((t_end - t0) / h)

# Allocate and initialize the RK4 trajectory.
t_values_rk4 = np.linspace(t0, t_end, N + 1)
y_values_rk4 = np.zeros((N + 1, len(y0)))
y_values_rk4[0] = y0

# Allocate and initialize the RK45 trajectory and error estimates.
t_values_rk45 = np.linspace(t0, t_end, N + 1)
y_values_rk45 = np.zeros((N + 1, len(y0)))
error_values_rk45 = np.zeros((N + 1, len(y0)))
y_values_rk45[0] = y0

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

# Gráficos dos resultados
plt.figure(figsize=(12, 6))
# Plot RK4 results
plt.subplot(1, 2, 1)
plt.plot(t_values_rk4, y_values_rk4[:, 0], label='RK4 - Posição (x)')
plt.plot(t_values_rk45, y_values_rk45[:, 0], label='RK45 - Posição (x)', linestyle='--')
plt.title('Posição ao longo do tempo')
plt.xlabel('Tempo (s)')
plt.ylabel('Posição (x)')
plt.legend()

plt.subplot(1, 2, 2)
# plot the error estimates for RK45
plt.plot(t_values_rk45, error_values_rk45[:, 0], label='Erro estimado RK45 - Posição (x)', color='orange')
plt.title('Erro estimado do método RK45 ao longo do tempo')
plt.xlabel('Tempo (s)')
plt.ylabel('Erro estimado')
plt.legend()

plt.tight_layout()
plt.show()
