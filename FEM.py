import numpy as np
import matplotlib.pyplot as plt

# Material
h_b = 10e-2  # Höhe des Querschnitts der Balken
b_b = 10e-2  # Breite des Querschnitts der Balken
E = 210e9  # E-Modul in Pa
A = h_b * b_b  # Querschnittsfläche der Balken
rho = 7850  # Dichte in kg/m^3
g = 9.81

# Nodes und Bars
nodes = []
bars = []
# Turm
height = 10
numSegmentsT = 10
hs = height / numSegmentsT  # Höhe jedes Segments

# Ausleger
offset = 4 * numSegmentsT - 5
numSegmentsA = 10
length = 10
ls = length / numSegmentsA

# Erstelle gestapelte Nodes Segmente des Turms
for i in range(numSegmentsT):
    nodes.append([0, 0, i * hs])  # Left Top
    nodes.append([hs, 0, i * hs])  # Right Top
    nodes.append([0, hs, i * hs])  # Left Bottom
    nodes.append([hs, hs, i * hs])  # Right Bottom

# Erstelle die Nodes des Auslegers in positive x Richtung
# for i in range(1, numSegmentsA):
#     nodes.append([hs + i * ls, 0, (numSegmentsT - 1) * hs])  # Left Top
#     nodes.append([hs + i * ls, hs, (numSegmentsT - 1) * hs])  # Right Top
#     nodes.append([hs + i * ls, 0, (numSegmentsT - 2) * hs])  # Left Bottom
#     nodes.append([hs + i * ls, hs, (numSegmentsT - 2) * hs])  # Right Bottom

# Turm
# x- und y-Richtung (LT für Left Top usw.)
for i in range(numSegmentsT):
    bars.append([4 * i, 4 * i + 1])  # LT -> RT
    bars.append([4 * i + 2, 4 * i + 3])  # LB -> RB
    bars.append([4 * i, 4 * i + 2])  # LT -> LB
    bars.append([4 * i + 1, 4 * i + 3])  # RT -> RB

# z-Richtung
for i in range(numSegmentsT - 1):
    bars.append([4 * i, 4 * i + 4])  # LT
    bars.append([4 * i + 1, 4 * i + 5])  # RT
    bars.append([4 * i + 2, 4 * i + 6])  # LB
    bars.append([4 * i + 3, 4 * i + 7])  # RB

# Kreuzstreben (+1 jeweils für nächste Stufe)
for i in range(numSegmentsT - 1):
    bars.append([4 * i, 4 * i + 5])  # LT -> RT+1
    bars.append([4 * i + 1, 4 * i + 4])  # RT -> RT+1
    bars.append([4 * i + 2, 4 * i + 7])  # LB -> RB+1
    bars.append([4 * i + 3, 4 * i + 6])  # RB -> LB+1
    bars.append([4 * i + 1, 4 * i + 7])  # RT -> RB+1
    bars.append([4 * i + 3, 4 * i + 5])  # RB -> RT+1
    bars.append([4 * i + 2, 4 * i + 4])  # LB -> LT+1
    bars.append([4 * i, 4 * i + 6])  # LT -> LB+1

# Ausleger

# python list zu np.array
nodes = np.array(nodes).astype(float)
bars = np.array(bars)

# Äußere Kräfte
P = np.zeros_like(nodes)
P[16, 0] = 1
P[17, 0] = 1
P[18, 0] = 1
P[19, 0] = 1

# Lager Verschiebung
Ur = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # 4 Festlager = 4*3 blockierte Freiheitsgrade

# Freiheitsgrade (1 = beweglich, 0 = fest) # evtl. booleans benutzen?
DOFCON = np.ones_like(nodes).astype(int)
# Festlager
DOFCON[0, :] = 0
DOFCON[1, :] = 0
DOFCON[2, :] = 0
DOFCON[3, :] = 0


# %% Truss structural analysis
def TrussAnalysis():
    NN = len(nodes)
    NE = len(bars)
    DOF = 3  # weil wir uns in 3D befinden
    NDOF = DOF * NN  # Gesamtanzahl der Freihetsgrade

    # Geometrie
    d = nodes[bars[:, 1], :] - nodes[bars[:, 0], :]  # Endknoten - Anfangsknoten
    L = np.sqrt((d ** 2).sum(axis=1))  # Länge der Balken (Euklidische Norm)
    orientations = d.T / L  # Richtungsvektoren der Balken (Transponieren für Dimensionen)
    trans = np.concatenate((-orientations.T, orientations.T), axis=1)  # Transformationsvektor lokal -> global

    # Steifigkeitsmatrix
    K = np.zeros([NDOF, NDOF])
    for k in range(NE):
        aux = DOF * bars[k, :]
        # Indices der vom Element betroffenen Knoten für Position der Summierung
        index = np.r_[aux[0]:aux[0] + DOF, aux[1]:aux[1] + DOF]
        # lokale Steifigkeiten, np.newaxis für 0 Zeile(4x4)
        ES = np.dot(trans[k][np.newaxis].T * E * A, trans[k][np.newaxis]) / L[k]
        # Globale Steifigkeiten durch Summierung der Einzelsteifigkeiten, Position !
        K[np.ix_(index, index)] = K[np.ix_(index, index)] + ES

    freeDOF = DOFCON.flatten().nonzero()[0]  # Prüfe, welche Knoten FG > 0 haben
    supportDOF = (DOFCON.flatten() == 0).nonzero()[0]  # Knoten mit Lagern
    Kff = K[np.ix_(freeDOF, freeDOF)]  # Vollkommen Bewegliche knoten
    Kfr = K[np.ix_(freeDOF, supportDOF)]  # Teilweise bewegliche Knoten
    Krf = Kfr.T
    Krr = K[np.ix_(supportDOF, supportDOF)]  # für die Lagerkräfte
    # Kraftmatrix passend zu K mit nicht null Einträgen, wie oben definiert
    Pf = P.flatten()[freeDOF]
    # Uf = np.linalg.solve(Kff,Pf)
    Uf = np.linalg.lstsq(Kff, Pf)[
        0]  # Deformation an jedem Freiheitsgrad # least squares damit auch überbestimmte Systeme fkt.
    U = DOFCON.astype(float).flatten()
    U[freeDOF] = Uf
    U[supportDOF] = Ur
    U = U.reshape(NN, DOF)
    # Verschiebungsvektor für die einzelnen Elemente
    u = np.concatenate((U[bars[:, 0]], U[bars[:, 1]]), axis=1)
    N = E * A / L[:] * (trans[:] * u[:]).sum(axis=1)  # interne Kräfte
    R = (Krf[:] * Uf).sum(axis=1) + (Krr[:] * Ur).sum(axis=1)  # Reaktionskräfte
    R = R.reshape(4, DOF)  # 4 ist die Anzahl der Lager
    return np.array(N), np.array(R), U


def Plot(nodes, c, lt, lw, lg):
    plt.subplot(projection='3d')
    plt.gca().set_aspect('auto')
    # plt.gca(projection='3d')
    for i in range(len(bars)):
        # Jeweilige Start und Endkoordiante
        xi, xf = nodes[bars[i, 0], 0], nodes[bars[i, 1], 0]
        yi, yf = nodes[bars[i, 0], 1], nodes[bars[i, 1], 1]
        zi, zf = nodes[bars[i, 0], 2], nodes[bars[i, 1], 2]
        # Plotte die Elemente
        line, = plt.plot([xi, xf], [yi, yf], [zi, zf], color=c, linestyle=lt, linewidth=lw)
    line.set_label(lg)
    plt.legend(prop={'size': 7})


# %% Result
N, R, U = TrussAnalysis()
print('Axial Forces (positive = tension, negative = compression)')
# Anschaulichkeit
print(N[np.newaxis].T)
print('Reaction Forces (positive = upward, negative = downward')
print(R)
print('Deformation at nodes')
print(U)
Plot(nodes, 'gray', '--', 1, 'Undeformed')
scale = 5
# Berechne die neue Position der Knoten
Dnodes = U * scale + nodes
Plot(Dnodes, 'red', '-', 2, 'Deformed')
plt.show()
# plt.savefig('fig-1.png', dpi=300)
