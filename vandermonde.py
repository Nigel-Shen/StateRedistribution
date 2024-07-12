import numpy as np
from scipy.special import legendre

def get_lagrange_legendre(Q: int, N: int) -> tuple:
    ts = np.linspace(0, 1, Q + 1, endpoint=True)
    xs, ws = np.polynomial.legendre.leggauss(N)
    xs = (xs + 1) / 2
    ws = ws / 2
    M = np.zeros((N, Q + 1))
    dM = np.zeros((N, Q + 1))
    for j in range(Q + 1):
        y = 1
        z = 0
        for i in range(len(ts)):
            if i != j:
                y = y * (xs - ts[i]) / (ts[j] - ts[i])
                z = z + 1 / (xs - ts[i])
        M[:, j] = y
        dM[:, j] = y * z
    return M, dM, ws

def get_legendre_vandermonde(xs, ys, P: int, **kwargs) -> tuple | np.ndarray:
    get_gradient = kwargs.get("get_gradient", False)
    
    V = np.zeros((len(xs), (P + 1) * (P + 2) // 2))
    if get_gradient:
        dxV = np.zeros((len(xs), (P + 1) * (P + 2) // 2))
        dyV = np.zeros((len(xs), (P + 1) * (P + 2) // 2))
    col = 0
    for j in range(P + 1):
        for k in range(P + 1 - j):
            V[:, col] = legendre(j)(xs) * legendre(k)(ys) * np.sqrt((2 * j + 1)*(2 * k + 1))
            if get_gradient:
                dxV[:, col] = np.polyder(legendre(j))(xs) * legendre(k)(ys) * np.sqrt((2 * j + 1)*(2 * k + 1))
                dyV[:, col] = legendre(j)(xs) * np.polyder(legendre(k))(ys) * np.sqrt((2 * j + 1)*(2 * k + 1))
            col += 1
    
    if get_gradient:
        return V, dxV, dyV
    else:
        return V