import numpy as np

class Truss:
    """
    This class represents a simple truss that can be passed to the fem class to make truss analyses
    """

    def __init__(self, nodes, bars, A, rho, E):
        self.nodes = np.array(nodes).astype(float)
        self.bars = np.array(bars)
        self.F = np.zeros_like(nodes) # externe Kräfte
        self.supports = np.ones_like(nodes).astype(int) # Freiheitsgrade (1 = beweglich, 0 = fest) # evtl. booleans benutzen?
        self.Ur = []
        
        # Material
        self.A = A
        self.rho = rho
        self.E = E

        self._computeLengths()
        self._computeOrientations()
        self._computeMass()
    
    # Geometrie
    def _computeLengths(self):
        d = self.nodes[self.bars[:, 1], :] - self.nodes[self.bars[:, 0], :]  # Endknoten - Anfangsknoten
        self.lengths = np.sqrt((d ** 2).sum(axis=1))  # Länge der Balken (Euklidische Norm)

    def _computeOrientations(self):
        d = self.nodes[self.bars[:, 1], :] - self.nodes[self.bars[:, 0], :]  # Endknoten - Anfangsknoten
        self.orientations = d.T / self.lengths  # Richtungsvektoren der Balken (Transponieren für Dimensionen)

    # Physik
    def _computeMass(self):
        self.mass = self.lengths*self.A*self.rho

    def addSupport(self, node, x, y, z):
        """
        node:   id of node where support should be added
        x:      displacement in x
        y:      displacement in y
        z:      displacement in z
        """
        self.supports[node, 0] = x
        self.supports[node, 1] = y
        self.supports[node, 2] = z

        self.Ur.append(self.supports[node][0])
        self.Ur.append(self.supports[node][1])
        self.Ur.append(self.supports[node][2])

    def addExternalForce(self, node, x, y, z):
        """
        Adds an external force. If there already is an external force at the given node,
        it will add one on top of it.
        node:   id of node where support should be added
        x:      x-component 
        y:      y-component 
        z:      z-component 
        """
        self.F[node, 0] = self.F[node, 0] + x
        self.F[node, 1] = self.F[node, 1] + y
        self.F[node, 2] = self.F[node, 2] + z
    