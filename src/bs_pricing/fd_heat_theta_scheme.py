import numpy as np
from .transforms import h_initial, g_left, g_right, reconstruct_V_from_u
from .linalg import solve_tridiagonal

def solve_heat_theta(
    K=100.0, r=0.02, sigma=0.2, T_days=180, option="call",
    x_min=-8.0, x_max=8.0, Nx=400, Nt=4000,
    theta=0.0, cfl_check=True
):
    """
    Solve u_tau = u_xx on tau in [0, Theta], x in [x_min,x_max] with Dirichlet boundaries.
    theta = 0 (explicit), 1 (implicit), 0.5 (Crank–Nicolson).
    Returns dict containing S-grid and V(0,S).
    """
    T = T_days / 365.0
    q = 2.0 * r / (sigma**2)
    Theta = 0.5 * sigma**2 * T

    x = np.linspace(x_min, x_max, Nx + 1)
    dx = x[1] - x[0]
    tau = np.linspace(0.0, Theta, Nt + 1)
    dtau = tau[1] - tau[0]
    alpha = dtau / (dx**2)

    if cfl_check and abs(theta) < 1e-15 and alpha > 0.5:
        raise ValueError(f"Explicit unstable: alpha={alpha:.4g} > 1/2. Increase Nt or decrease Nx.")

    u = h_initial(x, K=K, q=q, option=option)
    u[0]  = g_left(0.0, x_min, q, option)
    u[-1] = g_right(0.0, x_max, q, option)

    n_int = Nx - 1  # i=1..Nx-1

    if theta > 0:
        a = np.zeros(n_int)
        b = np.zeros(n_int)
        c = np.zeros(n_int)
        a[1:]  = -theta * alpha
        b[:]   = 1.0 + 2.0 * theta * alpha
        c[:-1] = -theta * alpha

    for m in range(Nt):
        tau_m   = tau[m]
        tau_mp1 = tau[m+1]

        Lm = g_left(tau_m,   x_min, q, option)
        Rm = g_right(tau_m,  x_max, q, option)
        Lp = g_left(tau_mp1, x_min, q, option)
        Rp = g_right(tau_mp1,x_max, q, option)

        if abs(theta) < 1e-15:
            u_new = u.copy()
            u_new[0], u_new[-1] = Lp, Rp
            u_new[1:-1] = alpha*u[2:] + (1.0-2.0*alpha)*u[1:-1] + alpha*u[:-2]
            u = u_new
        else:
            lap_u_m = u[2:] - 2.0*u[1:-1] + u[:-2]
            rhs = u[1:-1] + (1.0-theta)*alpha*lap_u_m
            rhs[0]  += theta*alpha*Lp + (1.0-theta)*alpha*Lm
            rhs[-1] += theta*alpha*Rp + (1.0-theta)*alpha*Rm

            u_int_next = solve_tridiagonal(a, b, c, rhs)
            u_new = u.copy()
            u_new[0], u_new[-1] = Lp, Rp
            u_new[1:-1] = u_int_next
            u = u_new

    S_grid = K * np.exp(x)
    V0_grid = reconstruct_V_from_u(u, x, Theta, K, q)

    return {
        "x": x, "S": S_grid, "u_final": u,
        "V0": V0_grid, "Theta": Theta,
        "params": {"q": q, "dx": dx, "dtau": dtau, "alpha": alpha, "theta": theta}
    }
