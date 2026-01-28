import numpy as np

def european_payoff(x, K, payoff = 'call'):
    if payoff == 'call':
        return max(x - K, 0)
    elif payoff == 'put':
        return max(K-x, 0)
    else:
        raise ValueError("payoff must be 'call' or 'put'")

def binomial_tree(S_0, K, T, r, N, sigma, payoff = 'call'):
    """
    Price a European option using the Cox–Ross–Rubinstein (CRR) binomial tree model.

    This function builds a recombining binomial tree for the underlying asset price,
    computes the terminal payoff of a European option, and performs backward induction
    under the risk-neutral measure to obtain the option price at time 0.

    Parameters
    ----------
    S_0 : float
        Initial price of the underlying asset.
    K : float
        Strike price of the option.
    T : float
        Time to maturity (in years).
    r : float
        Constant risk-free interest rate.
    N : int
        Number of time steps in the binomial tree.
    sigma : float
        Volatility of the underlying asset.
    payoff : str, optional
        Type of option payoff ('call' or 'put'). Default is 'call'.

    Returns
    -------
    float
        Option price at time t = 0.
    """
    delta_t = T/N
    disc = np.exp(-r * delta_t)
    u = np.exp(sigma * np.sqrt(delta_t))
    d = np.exp(-sigma * np.sqrt(delta_t))
    p = (np.exp(r * delta_t)-d)/(u-d)
    final_stock_prices = [S_0 * (u**j) * (d**(N-j)) for j in range(N+1)]
    Values = [european_payoff(s, K, payoff=payoff) for s in final_stock_prices]
    for i in range(N-1,-1,-1):
        last_layer_values = [disc * (p*Values[j+1] + (1-p) * Values[j]) for j in range(i+1)]
        Values = last_layer_values
    return Values[0]