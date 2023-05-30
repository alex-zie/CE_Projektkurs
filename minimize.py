from crane import crane
from fem import FEM
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

E = 210e9  # E-Modul in Pa
rho = 7850  # Dichte in kg/m^
price = 10


def tension(x: np.ndarray):
    myCrane = crane(1, x[0], x[1], x[2], x[3], x[4], rho, E)
    print(x)
    for i in range(-1, -5, -1):
        myCrane.addExternalForce(i, 0, 0, -500e3 / 4)
    fem = FEM(myCrane)
    # N, R, U = fem.TrussAnalysis()
    # t = np.max(N[np.newaxis]) / x[4]
    t = np.max(fem.getTension())
    print(t / 1e9)
    return -t


def cost(x: np.ndarray):
    print(x)
    myCrane = crane(1, x[0], x[1], x[2], x[3], x[4], rho, E)
    return np.sum(myCrane.mass)


if __name__ == "__main__":
    print(tension(np.array([7.5, 7.5, 1, 1, 0.0225])))
    cons = (
        {'type': 'ineq', 'fun': lambda x: tension(x) + 0.2e9},  # tension < 0.2e9
        {'type': 'ineq', 'fun': lambda x: x[0] - 2 * x[2]},
        {'type': 'ineq', 'fun': lambda x: x[1] - 2 * x[3]})
    res = minimize(cost, np.array([10, 10, 3, 3, 0.0225]),
                   bounds=((5, 10), (5, 10), (0.5, 2), (0.5, 2), (2.5e-3, 6.25e-2)),
                   method='SLSQP', constraints=cons, tol=1e-12)
    print(res)
