import numpy as np
import matplotlib.pyplot as plt

E = 1e4  # E-Modul
A = 0.111  # Querschnittsfläche der Balken

nodes = []
bars = []

nodes.append([-37.5,0,200])
nodes.append([37.5,0,200])
nodes.append([-37.5,37.5,100])
nodes.append([37.5,37.5,100])
nodes.append([37.5,-37.5,100])
nodes.append([-37.5,-37.5,100])
nodes.append([-100,100,0])
nodes.append([100,100,0])
nodes.append([100,-100,0])
nodes.append([-100,-100,0])

bars.append([0,1])
bars.append([3,0])
bars.append([2,1])
bars.append([4,0])
bars.append([5,1])
bars.append([3,1])
bars.append([4,1])
bars.append([2,0])
bars.append([5,0])
bars.append([5,2])
bars.append([4,3])
bars.append([2,3])
bars.append([5,4])
bars.append([9,2])
bars.append([6,5])
bars.append([8,3])
bars.append([7,4])
bars.append([6,3])
bars.append([7,2])
bars.append([9,4])
bars.append([8,5])
bars.append([9,5])
bars.append([6,2])
bars.append([7,3])
bars.append([8,4])

nodes = np.array(nodes).astype(float)
bars = np.array(bars)

#  Äußere Kräfte
P = np.zeros_like(nodes)
P[0,2] = -10
P[1,2] = -10
# P[-3,2] = -1e9
# P[-4,2] = -1e9

# Lager Verschiebung
Ur = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # 4 Festlager = 4*3 blockierte Freiheitsgrade

#Freiheitsgrade (1 = beweglich, 0 = fest)
DOFCON = np.ones_like(nodes).astype(int)
# Festlager
DOFCON[6,:] = 0
DOFCON[7,:] = 0
DOFCON[8,:] = 0
DOFCON[9,:] = 0

#%% Truss structural analysis
def TrussAnalysis():
    NN = len(nodes)
    NE = len(bars)
    DOF = 3  # weil wir uns in 3D befinden
    NDOF = DOF * NN

    # Geometrie
    d = nodes[bars[:, 1], :] - nodes[bars[:, 0], :]  # Endknoten - Anfangsknoten
    L = np.sqrt((d ** 2).sum(axis=1))  # Länge der Balken (Euklidische Norm)
    orientations = d.T / L  # Richtungsvektoren der Balken
    trans = np.concatenate((-orientations.T, orientations.T), axis=1)  # Transformationsvektor

    # Steifigkeitsmatrix
    K = np.zeros([NDOF, NDOF])
    for k in range(NE):
        aux = 3 * bars[k, :]
        index = np.r_[aux[0]:aux[0] + 3, aux[1]:aux[1] + 3]
        ES = np.dot(trans[k][np.newaxis].T * E * A, trans[k][np.newaxis]) / L[k]  # lokale Steifigkeiten
        K[np.ix_(index, index)] = K[np.ix_(index, index)] + ES

    freeDOF = DOFCON.flatten().nonzero()[0]
    supportDOF = (DOFCON.flatten() == 0).nonzero()[0]
    Kff = K[np.ix_(freeDOF, freeDOF)]
    Kfr = K[np.ix_(freeDOF, supportDOF)]
    Krf = Kfr.T
    Krr = K[np.ix_(supportDOF, supportDOF)]  # für die Lagerkräfte
    Pf = P.flatten()[freeDOF]
    # Uf = np.linalg.solve(Kff,Pf) # Deformation an jedem Freiheitsgrad
    Uf = np.linalg.lstsq(Kff, Pf)[0]  # Deformation an jedem Freiheitsgrad
    U = DOFCON.astype(float).flatten()
    U[freeDOF] = Uf
    U[supportDOF] = Ur
    U = U.reshape(NN, DOF)
    u = np.concatenate((U[bars[:, 0]], U[bars[:, 1]]), axis=1)
    N = E * A / L[:] * (trans[:] * u[:]).sum(axis=1)  # externe Kräfte
    R = (Krf[:] * Uf).sum(axis=1) + (Krr[:] * Ur).sum(axis=1)  # Reaktionskräfte
    R = R.reshape(4, DOF)  # 4 ist die Anzahl der Lager
    return np.array(N), np.array(R), U


def Plot(nodes, c, lt, lw, lg):
    plt.subplot(projection='3d')
    #plt.gca().set_aspect('equal')
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
