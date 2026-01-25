import numpy as np

def payoff(S, K, option = 'call'):
    return np.maximum(S-K,0.0) if option == 'call' else np.maxmum(K-S,0.0)


#Initial condition function
def h_initial(x, K, q, option = 'call'):
    S = K * np.exp(x)
    return np.exp(0.5 * (q - 1.0) * x) * payoff(S, K, option=option) / K

#left asymptotic condition approximation
def g_left(tau, x_min, q, option = 'call'):
    if option =='call':
       return 0.0
    return np.exp(0.5 * (q - 1.0) * x_min + 0.25 * (q - 1.0)**2 * tau)

#right asymptotic condition approximation
def g_right(tau, x_max, q, option="call"):
    if option == "put":
        return 0.0
    return np.exp(0.5 * (q + 1.0) * x_max + 0.25 * (q + 1.0)**2 * tau)

def reconstruct_V_from_u(u, x_grid, tau, K, q):
    # recall: V(Ke^x, T-2tau/sigma^2) = K exp(-1/2(q-1)x - 1/4(q+1)^2 tau) u(tau,x)
    return K * np.exp(-0.5 * (q - 1.0) * x_grid - 0.25 * (q + 1.0)**2 * tau) * u
