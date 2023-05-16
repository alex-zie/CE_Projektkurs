from crane import crane
from fem import FEM
import numpy as np
import matplotlib.pyplot as plt


def optimize_area(E, rho):
    t = np.array([])
    a = np.array([])
    for h_b in range(5, 26):
        for b_b in range(h_b, 26):
            crane1 = crane(10, 10, 1, 1, (h_b * b_b) * 10 ** (-4), rho, E)
            fem = FEM(crane1)
            N, R, U = fem.TrussAnalysis()
            t = np.append(t, np.max(N[np.newaxis]) / ((h_b * b_b) * 10 ** (-4)))
            a = np.append(a, np.array([(h_b * b_b) * 10 ** (-4), h_b, b_b]))
            a = a.reshape(-1, 3)
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
    print(np.where(t < 0.2e9)[0])
E = 210e9  # E-Modul in Pa
rho = 7850  # Dichte in kg/m^
g = 9.81  # m/s^2
optimize_area(E, g)
