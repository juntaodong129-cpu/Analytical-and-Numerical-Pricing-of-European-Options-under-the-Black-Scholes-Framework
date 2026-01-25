import numpy as np
from scipy.stats import norm

def bs_price(S, K, r, sigma, T, option = 'call'):
    """
    Black–Scholes closed-form price for European options.

    Parameters
    ----------
    S : float or array-like
    K : float
    r : float
    sigma : float
    T : float (year fraction)
    option : {"call","put"}
    """
    if K <= 0:
        raise ValueError("K must be > 0")
    if T < 0:
        raise ValueError("T must be >=0")
    if sigma < 0:
       raise ValueError("sigma must be >= 0")
    if option not in ("call", "put"):
        raise ValueError("option must be 'call' or 'put'")
    
    
    
    S_arr = np.asarray(S, dtype=float)
    if np.any(S_arr < 0):
        raise ValueError("S must be >= 0")
    
    # two special cases
    if T == 0:
        out = np.maximum(S_arr - K, 0.0) if option == "call" else np.maximum(K - S_arr, 0.0)
        return float(out) if np.isscalar(S) else out
    
    if sigma == 0:
       ST = S_arr * np.exp(r*T)
       disc = np.exp(-r*T)
       out = disc * (np.maximum(ST - K, 0.0) if option == "call" else np.maximum(K - ST, 0.0))
       return float(out) if np.isscalar(S) else out

    # normal cases
    sqrtT = np.sqrt(T)

    d1 = (np.log(S_arr/K)+ (r + 0.5 * sigma**2) * T)/(sigma * sqrtT)
    d2 = d1 - sigma * sqrtT

    if option =='call':
        out = S_arr * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)
    else:
        out = K * np.exp(-r*T) * norm.cdf(-d2) - S_arr * norm.cdf(-d1)
    
    return float(out) if np.isscalar(S) else out