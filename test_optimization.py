import numpy as np
import scipy.optimize as sp


def objective(x):
    return x[0] ** 2.0 + x[1] ** 2.0


def cons(x):
    return x[1] - 1


print(sp.minimize(objective, np.array([2, 2]),bounds=[[1, 2], [-1, -0.5]], constraints={'type': 'ineq', 'fun': cons}))
#
