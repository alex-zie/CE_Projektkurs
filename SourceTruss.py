import numpy as np
import matplotlib.pyplot as plt

# %% Input truss structure data
E = 1e4
A = 0.111

nodes = []
bars = []

nodes.append([0, 0, 0])
nodes.append([1, 0, 0])
nodes.append([0.5, 1, 0])
nodes.append([0.5, 0.5, 1])
bars.append([0, 1])
bars.append([1, 2])
bars.append([2, 1])
bars.append([0, 3])
bars.append([1, 3])
bars.append([2, 3])
nodes = np.array(nodes).astype(float)
bars = np.array(bars)

# Applied forces
P = np.zeros_like(nodes)
P[3, 2] = -100

# Support Displacement
Ur = [0, 0, 0, 0, 0, 0, 0, 0, 0]
# Condition of DOF (1 = free, 0 = fixed)
DOFCON = np.ones_like(nodes).astype(int)
DOFCON[0, :] = 0
DOFCON[1, :] = 0
DOFCON[2, :] = 0
# %% Truss structural analysis
def TrussAnalysis():
    NN = len(nodes)
    NE = len(bars)
    DOF = 3
    NDOF = DOF * NN

    # structural analysis
    d = nodes[bars[:, 1], :] - nodes[bars[:, 0], :]
    L = np.sqrt((d ** 2).sum(axis=1))
    angle = d.T / L
    a = np.concatenate((-angle.T, angle.T), axis=1)
    K = np.zeros([NDOF, NDOF])
    for k in range(NE):
        aux = 3 * bars[k, :]
        index = np.r_[aux[0]:aux[0] + 3, aux[1]:aux[1] + 3]

        ES = np.dot(a[k][np.newaxis].T * E * A, a[k][np.newaxis]) / L[k]
        K[np.ix_(index, index)] = K[np.ix_(index, index)] + ES

    freeDOF = DOFCON.flatten().nonzero()[0]
    supportDOF = (DOFCON.flatten() == 0).nonzero()[0]
    Kff = K[np.ix_(freeDOF, freeDOF)]
    print("Kff-determinate = " + str(np.linalg.det(Kff)))
    Kfr = K[np.ix_(freeDOF, supportDOF)]
    Krf = Kfr.T
    Krr = K[np.ix_(supportDOF, supportDOF)]
    Pf = P.flatten()[freeDOF]
    Uf = np.linalg.solve(Kff, Pf)
    U = DOFCON.astype(float).flatten()
    U[freeDOF] = Uf
    U[supportDOF] = Ur
    U = U.reshape(NN, DOF)
    u = np.concatenate((U[bars[:, 0]], U[bars[:, 1]]), axis=1)
    N = E * A / L[:] * (a[:] * u[:]).sum(axis=1)
    R = (Krf[:] * Uf).sum(axis=1) + (Krr[:] * Ur).sum(axis=1)
    R = R.reshape(3, DOF)
    return np.array(N), np.array(R), U


def Plot(nodes, c, lt, lw, lg):
    plt.subplot(projection='3d')
    plt.gca().set_aspect('equal')
    for i in range(len(bars)):
        xi, xf = nodes[bars[i, 0], 0], nodes[bars[i, 1], 0]
        yi, yf = nodes[bars[i, 0], 1], nodes[bars[i, 1], 1]
        zi, zf = nodes[bars[i, 0], 2], nodes[bars[i, 1], 2]
        line, = plt.plot([xi, xf], [yi, yf], [zi, zf], color=c, linestyle=lt, linewidth=lw)
    line.set_label(lg)
    plt.legend(prop={'size': 14})


# %% Result
N, R, U = TrussAnalysis()
print('Axial Forces (positive = tension, negative = compression)')
print(N[np.newaxis].T)
print('Reaction Forces (positive = upward, negative = downward')
print(R)
print('Deformation at nodes')
print(U)
Plot(nodes, 'gray', '--', 1, 'Undeformed')
scale = 5
Dnodes = U * scale + nodes
Plot(Dnodes, 'red', '-', 2, 'Deformed')
plt.show()
# plt.savefig('fig-1.png', dpi=300)
