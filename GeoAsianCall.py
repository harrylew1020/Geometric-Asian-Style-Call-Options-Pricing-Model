import numpy as np
from scipy.integrate import quad as quad


def GeomAsianCall(S0, v0, theta, sigma, kappa, rho, r, n, T, K):
    call = np.exp(-r * T) * ((psi(1, 0, S0, v0, theta, sigma, kappa, rho, r, n, T) - K) * 0.5 + 1 / np.pi *
                             GeoIntegral(S0, v0, theta, sigma, kappa, rho, r, n, T, K))
    return call


def GeoIntegral(S0, v0, theta, sigma, kappa, rho, r, n, T, K):
    ret = quad(integrand, 0, 1e+5, args=(S0, v0, theta, sigma, kappa, rho, r, n, T, K))
    return ret


def integrand(x, S0, v0, theta, sigma, kappa, rho, r, n, T, K):
    A = psi(1 + 1j * x, 0, S0, v0, theta, sigma, kappa, rho, r, n, T)
    B = psi(x * 1j, 0, S0, v0, theta, sigma, kappa, rho, r, n, T)
    C = np.exp(-1j * x * np.log(K)) / (1j * x)
    ret = np.real((A - K * B) * C)
    return ret


def psi(s, w, S0, v0, theta, sigma, kappa, rho, r, n, T):
    a1 = 2 * v0 / sigma ** 2
    a2 = 2 * kappa * theta / sigma ** 2
    a3 = np.log(S0) + ((r * sigma - kappa * theta * rho) * T) / (2 * sigma) - (rho * v0) / sigma
    a4 = np.log(S0) - (rho * v0 / sigma) + (r - rho * kappa * theta / sigma) * T
    a5 = (kappa * v0 + kappa ** 2 * theta * T) / (sigma ** 2)

    hmat = np.zeros((n + 3, s.shape[1]))
    hmat[2] = 1
    hmat[3] = T * (kappa - w * rho * sigma) / 2
    nmat = np.linspace(1, n, n)
    A1 = 1 / (4 * nmat[1:, 0] * (nmat[1:, 0] - 1))
    A2 = -s ** 2 * sigma ** 2 * (1 - rho ** 2) * T ** 2
    A3 = (s * sigma * T * (sigma - 2 * rho * kappa) - 2 * s * w * sigma ** 2 * T * (1 - rho ** 2))
    A4 = T * (kappa ** 2 * T - 2 * s * rho * sigma - w * (2 * rho * kappa - sigma) * sigma * T
              - w ** 2 * (1 - rho ** 2) * sigma ** 2 * T)

    for i in range(4, n + 2):
        hmat[i] = A1[i - 3, 0] * (A2 * hmat[i - 3] + A3 * (T * hmat[i - 2])
                                  + A4 * hmat[i - 1])
        print(hmat)

    H = sum(hmat[2])
    h_tilde = (nmat / T) * hmat[3]
    H_tilde = sum(h_tilde)
    ret = np.exp(-a1 * (H_tilde / H) - a2 * np.log(H) + a3 * s + a4 * w + a5)
    return ret

