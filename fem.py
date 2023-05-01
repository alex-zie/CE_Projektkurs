import numpy as np
from truss import Truss
import matplotlib.pyplot as plt

class FEM:

    def __init__(self, truss: Truss):
        self.truss = truss
        self.NN = len(truss.nodes)
        self.NE = len(truss.bars)
        self.DOF = 3  # weil wir uns in 3D befinden
        self.NDOF = self.DOF * self.NN  # Gesamtanzahl der Freihetsgrade

    def TrussAnalysis(self):
        """
        returns: axial forces, reactional forces, displacements
        """
        E = self.truss.E 
        A = self.truss.A
        L = self.truss.lengths
        trans = np.concatenate((-self.truss.orientations.T, self.truss.orientations.T), axis=1)  # Transformationsvektor lokal -> global

        K = self.computeStiffnessMatrix(E, A, L, trans)

        freeDOF = self.truss.supports.flatten().nonzero()[0]  # Prüfe, welche Knoten FG > 0 haben
        supportDOF = (self.truss.supports.flatten() == 0).nonzero()[0]  # Knoten mit Lagern
        Kff = K[np.ix_(freeDOF, freeDOF)]  # Vollkommen Bewegliche knoten
        Kfr = K[np.ix_(freeDOF, supportDOF)]  # Teilweise bewegliche Knoten # sicher?
        Krf = Kfr.T
        Krr = K[np.ix_(supportDOF, supportDOF)]  # für die Lagerkräfte

        weights = np.zeros_like(self.truss.F)
        weights[:, 2] = -self.computeWeight()
        self.truss.addExternalForces(weights)
        F = self.truss.F.flatten()[freeDOF] # Kraftmatrix passend zu K mit nicht null Einträgen, wie oben definiert

        Uf = np.linalg.lstsq(Kff, F)[0]  # Deformation an jedem Freiheitsgrad # least squares damit auch überbestimmte Systeme fkt.
        U = self.truss.supports.astype(float).flatten()
        U[freeDOF] = Uf
        U[supportDOF] = self.truss.Ur
        U = U.reshape(self.NN, self.DOF)
        u = np.concatenate((U[self.truss.bars[:, 0]], U[self.truss.bars[:, 1]]), axis=1) # Verschiebungsvektor für die einzelnen Elemente
        N = E * A / L[:] * (trans[:] * u[:]).sum(axis=1)  # interne Kräfte
        R = (Krf[:] * Uf).sum(axis=1) + (Krr[:] * self.truss.Ur).sum(axis=1)  # Reaktionskräfte
        R = R.reshape(len(self.truss.Ur)//self.DOF, self.DOF)  # TODO Diese Zeile evtl. verbessern
        return np.array(N), np.array(R), U
    
    def computeStiffnessMatrix(self, E, A, L, trans):
        K = np.zeros([self.NDOF, self.NDOF])
        for k in range(self.NE):
            aux = self.DOF * self.truss.bars[k, :]
            # Indices der vom Element betroffenen Knoten für Position der Summierung
            index = np.r_[aux[0]:aux[0] + self.DOF, aux[1]:aux[1] + self.DOF]
            # lokale Steifigkeiten, np.newaxis für 0 Zeile(4x4)
            ES = np.dot(trans[k][np.newaxis].T * E * A, trans[k][np.newaxis]) / L[k]
            # Globale Steifigkeiten durch Summierung der Einzelsteifigkeiten, Position !
            K[np.ix_(index, index)] = K[np.ix_(index, index)] + ES
        return K
    
    def computeWeight(self):
        bars = self.truss.bars
        nOutgoingBars = np.zeros_like(self.truss.nodes[:, 0])
        # compute for each node (index) the number of outgoing bars
        for bar in bars:
            nOutgoingBars[bar] = nOutgoingBars[bar] + 1

        masses = self.truss.mass
        weights = np.zeros_like(self.truss.nodes[:, 0])
        # add a weight to the node equalling the mass of a connected bar
        for i in range(len(bars)):
            weights[bars[i]] = weights[bars[i]] + masses[i]

        # make it correct
        return 9.81*weights/nOutgoingBars

    def Plot(self, nodes, bars, c, lt, lw, lg):
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

    def plotNode(self, node):
        plt.subplot(projection='3d')
        plt.plot(self.truss.nodes[node, 0], self.truss.nodes[node, 1], self.truss.nodes[node, 2], 'bo')

    def plotPoint(self, point):
        plt.subplot(projection='3d')
        plt.plot(point[0], point[1], point[2], 'bo')