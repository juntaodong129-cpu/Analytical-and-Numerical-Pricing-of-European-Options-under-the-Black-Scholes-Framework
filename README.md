# Numerical Pricing of European Options  
### under the Black–Scholes Framework

**Author:** Juntao DONG  
**Date:** February 2026  

This project studies the pricing of European call and put options under the Black–Scholes model using:

- Closed-form analytical formulas  
- Finite-difference methods for the Black–Scholes PDE  
- Cox–Ross–Rubinstein (CRR) binomial tree method  

The complete mathematical derivations and numerical experiments are presented in:

> *Analytical and Numerical Pricing of European Options under the Black–Scholes Framework*


---

## 1. Model Setup

We consider a financial market where the risky asset follows:

dS_t = S_t (μ_t dt + σ dB_t)

with:

- Constant volatility σ > 0
- Constant interest rate r
- European payoff:
  - Call: (S_T − K)^+
  - Put:  (K − S_T)^+

Under the risk-neutral measure Q:

V_t = E^Q [ e^{-r(T-t)} f(S_T) | F_t ]


---

## 2. Closed-Form Black–Scholes Formula

For a European call:

C(t,S) = S N(d1) − K e^{-r(T−t)} N(d2)

For a European put:

P(t,S) = K e^{-r(T−t)} N(−d2) − S N(−d1)

where:

d1 = [ ln(S/K) + (r + 0.5σ²)(T−t) ] / (σ√(T−t))  
d2 = d1 − σ√(T−t)

We analyze price sensitivity with respect to:

- Spot S
- Volatility σ
- Interest rate r
- Time to maturity T


---

## 3. Finite-Difference Method

The Black–Scholes PDE:

∂_t v + ½ σ² S² ∂²_SS v + rS ∂_S v − r v = 0

After change of variables:

S = K e^x  
t = T − 2τ/σ²  

it becomes the heat equation:

∂_τ u = ∂²_xx u


### Numerical Schemes (θ-scheme)

u_i^{m+1} − u_i^m
------------------ = (1−θ)Δ_xx u^m + θ Δ_xx u^{m+1}
       Δτ

- θ = 0 → Explicit
- θ = 1 → Implicit Euler
- θ = 1/2 → Crank–Nicolson

Explicit stability condition:

Δτ / (Δx)² ≤ 1/2


### Numerical Experiments

- Comparison with analytical solution
- Error concentration near strike S = K
- Grid refinement study
- L² convergence analysis
- Log–log convergence plots

Main observations:

- Explicit scheme has largest error
- Implicit and Crank–Nicolson are significantly more accurate
- Convergence observed under grid refinement


---

## 4. Binomial Tree Method (CRR)

Time discretisation:

Δt = T / N  

Tree parameters:

u = e^{σ√Δt}  
d = e^{−σ√Δt}  

Risk-neutral probability:

p = (e^{rΔt} − d) / (u − d)

Backward recursion:

V_{i,j} = e^{−rΔt} [ p V_{i+1,j+1} + (1−p)V_{i+1,j} ]


### Experiments

- Comparison with Black–Scholes formula
- Error maximum near at-the-money region
- Convergence with N = 2^k, k = 4,...,10
- Log–log L² convergence plots

Observed convergence rate ≈ first order.


---

## 5. Numerical Parameters

Typical configuration:

- K = 100
- r = 2%
- σ = 20%
- T ∈ {5, 180, 360} days
- S ∈ {50, 70, 90, 100, 110, 130, 150}


---


