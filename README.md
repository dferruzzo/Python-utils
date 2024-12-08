# Python Utilities

This repository contains the implementation of 4th-order Runge-Kutta algorithms and other utilities related to control systems in Python.

## Description

The 4th-order Runge-Kutta method is a numerical technique used to solve ordinary differential equations (ODEs). This repository includes several implementations of the Runge-Kutta algorithm, including the standard method, a method with control inputs, and the Runge-Kutta-Fehlberg method for adaptive step size.

## Algorithms

### 1. Standard Runge-Kutta 4th Order (RK4)

- **Function:** `rk4(f, x0, t0, tf, h)`
- **Description:** Solves ODEs using the 4th-order Runge-Kutta method.
- **Parameters:**
  - `f`: Function representing the ODE.
  - `x0`: Initial state vector.
  - `t0`: Initial time.
  - `tf`: Final time.
  - `h`: Integration step size.
- **Returns:** Time vector `t` and state vector `x`.

### 2. Runge-Kutta 4th Order with Control Inputs (RK4U)

- **Function:** `rk4u(f, x0, t0, tf, h, gfunc, m)`
- **Description:** Solves ODEs with control inputs using the 4th-order Runge-Kutta method.
- **Parameters:**
  - `f`: Function representing the ODE with control inputs.
  - `x0`: Initial state vector.
  - `t0`: Initial time.
  - `tf`: Final time.
  - `h`: Integration step size.
  - `gfunc`: Control input function.
  - `m`: Number of control inputs.
- **Returns:** Time vector `t`, state vector `x`, and control input vector `u`.

### 3. Runge-Kutta-Fehlberg 4th-5th Order (RKF45)

- **Function:** `rkf45(f, x0, t0, tf, h)`
- **Description:** Solves ODEs using the adaptive step size 4th-5th order Runge-Kutta-Fehlberg method.
- **Parameters:**
  - `f`: Function representing the ODE.
  - `x0`: Initial state vector.
  - `t0`: Initial time.
  - `tf`: Final time.
  - `h`: Initial integration step size.
- **Returns:** Time vector `t`, state vector `x`, error estimates `erro`, and step sizes `hs`.

## Additional Functions

### Control Systems

- **Function:** `controlavel(A, B)`
- **Description:** Checks the controllability of the system defined by matrices `A` and `B`.

### Solutions to Riccati Equation

- **Function:** `P_analitico(tf, t, A, B, F, R, Q)`
- **Description:** Calculates the analytical solution to the Riccati equation.
- **Function:** `P_tf_inf(A, B, R, Q)`
- **Description:** Calculates the analytical solution to the Riccati equation as time approaches infinity.

### Sub-Optimal Control

- **Function:** `SDRE(A, B, R, Q, t0, tf, h, x0)`
- **Description:** Computes the sub-optimal trajectory for a given system using State-Dependent Riccati Equation (SDRE) method.

### Linear Quadratic Regulator

- **Function:** `lqr_regulador(A, B, Q, R, h, t0, tf, dt, x0)`
- **Description:** Computes the optimal control trajectory for a Linear Quadratic Regulator (LQR).

## Installation

To use the functions and algorithms in this repository, clone the repository and install the required dependencies.

```bash
git clone https://github.com/dferruzzo/Runge-Kutta.git
cd Runge-Kutta
pip install -r requirements.txt
```

## Usage

Import the required functions from the repository and use them as needed in your Python scripts.

```python
from myfunctions import rk4, rk4u, rkf45

# Example usage
t, x = rk4(f, x0, t0, tf, h)
```

## License

This project is licensed under the MIT License.
