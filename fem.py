import numpy as np
import crane
import matplotlib.pyplot as plt



def Plot(nodes,bars, c, lt, lw, lg):
    plt.subplot(projection='3d')
    plt.gca().set_aspect('equal')
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

def TrussAnalysis(crane: crane.Crane , P: np.ndarray, DOFCON: np.ndarray, Ur: np.ndarray, A, rho,E) -> (np.ndarray, np.ndarray, np.ndarray) :
    nodes,bars = crane.getNB()
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