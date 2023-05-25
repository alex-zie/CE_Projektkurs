import numpy as np
import scipy.optimize as sp


def objective(x):
    return x[0] ** 2.0 + x[1] ** 2.0


def cons(x):
    return x[1] - 1


print(sp.minimize(objective, np.array([1, 1]), bounds=[],constraints={'type': 'ineq', 'fun': cons}))
#
