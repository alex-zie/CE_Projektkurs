import numpy as np

class truss:
    """
    This class represents a simple truss that can be passed to the fem class to make truss analyses
    """

    def __init__(self, nodes, bars, A, rho, E):
        self.nodes = np.array(nodes).astype(float)
        self.bars = np.array(bars)
        self.F = np.zeros_like(nodes) # externe Kräfte
        self.supports = np.ones_like(nodes).astype(int) # Lager
        # self.Ur = ...
        
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