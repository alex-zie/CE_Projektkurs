from crane import crane
from fem import FEM
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp


def optimize_area(E, rho):
    t = np.array([])
    a = np.array([])
    # for h_l in np.linspace(0.7, 1.5 / np.sqrt(2), 10):
    for h_b in range(5, 26):
        crane1 = crane(10, 10, 1, 1, h_b ** 2 * 10 ** (-4), rho, E)
        fem = FEM(crane1)
        N, R, U = fem.TrussAnalysis()
        t = np.append(t, np.max(N[np.newaxis]) / (h_b ** 2 * 10 ** (-4)))
        a = np.append(a, np.array([h_b ** 2 * 10 ** (-4), h_b]))
        a = a.reshape(-1, 2)
    ind = np.argsort(t)
    t = np.sort(t)
    a = a[ind, :]
    plt.plot(a[:, 0], t, 'b', label="Maximale Spannung")
    plt.xlabel("A in m^2")
    plt.ylabel("Maximale Spannung in N/m^2")
    plt.plot(a[:, 0], np.ones(len(a)) * 0.2 * 10 ** 9, 'r', label="Streckgrenze")
    plt.xlim(25 * 10 ** (-4), 25 * 25 * 10 ** (-4))
    plt.legend()
    plt.show()
    # print(np.where(t < 0.2e9)[0])


def optimize_lengths(E, rho, A):
    t = np.array([])
    a = np.array([])

    for l in np.linspace(0.6, 1.5 / np.sqrt(2), 100):
        crane1 = crane(10, 10, l, l, A, rho, E)
        fem = FEM(crane1)
        N, R, U = fem.TrussAnalysis()
        t = np.append(t, np.max(N[np.newaxis]) / A)
        a = np.append(a, [l])
    plt.plot(a, t, 'b', label="Maximale Spannung")
    plt.xlabel("l in m")
    plt.ylabel("Maximale Spannung in N/m^2")
    plt.plot(a, np.ones(len(a)) * 0.2 * 10 ** 9, 'r', label="Streckgrenze")
    plt.xlim((0.7, 1.5 / np.sqrt(2)))
    plt.legend()
    plt.show()
    # print(np.where(t < 0.2e9)[0])


def f(x):
    E = 210e9  # E-Modul in Pa
    rho = 7850  # Dichte in kg/m^
    crane1 = crane(10, 10, x[0], x[0], x[1], rho, E)
    fem = FEM(crane1)
    weight = np.sum(fem.computeWeight() / 9.81)
    return rho * x[1] * weight


def cons(x):
    E = 210e9  # E-Modul in Pa
    rho = 7850  # Dichte in kg/m^
    crane1 = crane(10, 10, x[0], x[0], x[1], rho, E)
    fem = FEM(crane1)
    return 2e8 - fem.TrussAnalysis()[0]


l_min = 0.6
l_max = 1.5
A_min = 0.01
A_max = 0.75
o_y = 2e8
A = 0.018
E = 210e9  # E-Modul in Pa
rho = 7850  # Dichte in kg/m^
g = 9.81  # m/s^2
print(sp.minimize(f, np.array([0.8, 0.02]), bounds=[[A_min, A_max], [l_min, l_max]],
                  constraints={'type': 'ineq', 'fun': cons}))
