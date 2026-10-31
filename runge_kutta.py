
from numpy.linalg import norm
from numpy import inf

def rk4(f, t0, y0, h):
    """Advance an ordinary differential equation by one fourth-order RK step.

    Solves ``dy/dt = f(t, y)`` from ``(t0, y0)`` using the classical
    four-stage Runge-Kutta method and a fixed step size.

    Parameters:
    f : callable
        Function defining the differential equation ``dy/dt = f(t, y)``.
    t0 : float
        The initial time.
    y0 : float or numpy.ndarray
        State at ``t0``.
    h : float
        Fixed integration step size.

    Returns:
    tuple[float, float or numpy.ndarray]
        The time and state after one step.
    """
    k1 = f(t0, y0)
    k2 = f(t0 + (h/2), y0 + (h*k1/2))
    k3 = f(t0 + (h/2), y0 + (h*k2/2))
    k4 = f(t0 + h, y0 + (h*k3))

    t = t0 + h
    y = y0 + (h/6)*(k1 + 2*k2 + 2*k3 + k4)

    return t, y

def rk45(f, t0, y0, h):
    """Advance an ODE by one fixed Runge-Kutta-Fehlberg (RK45) step.

    Computes the fourth-order state estimate and the difference between the
    embedded fourth- and fifth-order estimates, which callers can use for
    local-error control.

    Parameters:
    f : callable
        Function defining the differential equation ``dy/dt = f(t, y)``.
    t0 : float
        The initial time.
    y0 : float or numpy.ndarray
        State at ``t0``.
    h : float
        Fixed integration step size.

    Returns:
    tuple[float, float or numpy.ndarray, float or numpy.ndarray]
        The time and fourth-order state estimate after one step, followed by
        the embedded local-error estimate.
    """
    k1 = f(t0, y0)
    k2 = f(t0 + (h/4), y0 + (h*k1/4))
    k3 = f(t0 + (3*h/8), y0 + (3*h*k1/32) + (9*h*k2/32))
    k4 = f(t0 + (12*h/13), y0 + (1932*h*k1/2197) - (7200*h*k2/2197) + (7296*h*k3/2197))
    k5 = f(t0 + h, y0 + (439*h*k1/216) - (8*h*k2) + (3680*h*k3/513) - (845*h*k4/4104))
    k6 = f(t0 + (h/2), y0 - (8*h*k1/27) + (2*h*k2) - (3544*h*k3/2565) + (1859*h*k4/4104) - (11*h*k5/40))

    t = t0 + h
    y = y0 + h*((25*k1/216) + (1408*k3/2565) + (2197*k4/4104) - (1*k5/5))
    error = h*((1*k1/360) - (128*k3/4275) - (2197*k4/75240) + (1*k5/50) + (2*k6/55))

    return t, y, error

def rk45_h_adj(f, t0, y0, h, epsilon):
    """Perform an RK45 step and calculate an adjusted next step size.

    The returned state is the fourth-order Runge-Kutta-Fehlberg estimate.
    The proposed next step uses the infinity norm of the difference between
    the embedded fourth- and fifth-order estimates and a safety factor of
    0.84.

    Parameters:
    f : callable
        Function defining the differential equation ``dy/dt = f(t, y)``.
    t0 : float
        The initial time.
    y0 : float or numpy.ndarray
        State at ``t0``.
    h : float
        Step size used for this integration step.
    epsilon : float
        Desired local-error tolerance used to calculate the next step size.

    Returns:
    tuple[float, float or numpy.ndarray, float or numpy.ndarray, float]
        The time and fourth-order state estimate after one step, the embedded
        local-error estimate, and the proposed step size for the next step.
    """
    k1 = f( t0, y0 );
    k2 = f( t0 + h/4, y0 + k1*h/4 );
    k3 = f( t0 + 3*h/8, y0 + 3*k1*h/32 + 9*k2*h/32 );
    k4 = f( t0 + 12*h/13, y0 + 1932*k1*h/2197 - 7200*k2*h/2197 + 7296*k3*h/2197 );
    k5 = f( t0 + h, y0 + 439*k1*h/216 - 8*k2*h + 3680*k3*h/513 - 845*k4*h/4104 );
    k6 = f( t0 + h/2, y0 - 8*k1*h/27 + 2*k2*h - 3544*k3*h/2565 + 1859*k4*h/4104 - 11*k5*h/40 );
    
    # RK4
    y_rk4 = y0 + h*( 25*k1/216 + 1408*k3/2565 + 2197*k4/4104 - k5/5 );
    y_rk5 = y0 + h*( 16*k1/135 + 6656*k3/12825 + 28561*k4/56430 - 9*k5/50 + 2*k6/55 );
    y = y_rk4;
    t = t0 + h;
    # estimativa do erro
    error = h*( k1/360 - 128*k3/4275 - 2197*k4/75240 + k5/50 + 2*k6/55 );
    # h ajustado
    alpha = 0.84; # Segundo a ref.
    p = 4; # RK order
    h_adj = alpha*h*(epsilon*h/norm(y_rk5 - y_rk4, inf))**(1/p);
    #
    return t, y, error, h_adj
