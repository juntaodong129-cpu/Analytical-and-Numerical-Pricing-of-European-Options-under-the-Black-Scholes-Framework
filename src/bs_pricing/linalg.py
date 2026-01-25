import numpy as np

def solve_tridiagonal(a, b, c, d):
    """
    Thomas algorithm for tridiagonal systems.
    a,b,c,d length n; a[0] unused, c[-1] unused.
    """
    n = len(b)
    ac, bc, cc, dc = map(np.array, (a, b, c, d))

    for i in range(1, n):
        w = ac[i] / bc[i-1]
        bc[i] -= w * cc[i-1]
        dc[i] -= w * dc[i-1]

    x = np.zeros(n, dtype=float)
    x[-1] = dc[-1] / bc[-1]
    for i in range(n-2, -1, -1):
        x[i] = (dc[i] - cc[i] * x[i+1]) / bc[i]
    return x
