import numpy as np


class Truss:
    """
    This class represents a simple truss that can be passed to the fem class to perfom a truss analyses
    """

    def __init__(self, nodes, bars, A, rho, E):
        """
        :param nodes: nodes of the truss
        :param bars: bars of the truss
        :param A: area of each element of the truss
        :param rho: density of the material of the truss
        :param E: E-Module of each element of the truss

        Constructor of this class initializes the truss with the given parameters and a fixed bearing
        on the bottom nodes of the tower
        """
        self.nodes = np.array(nodes).astype(float)
        self.bars = np.array(bars)
        self.F = np.zeros_like(nodes)  # external forces
        self.force_points = []  # indices of points, where external forces attack
        # degrees of freedom (1 = movable, 0 = fixed), 3 times for (x,y,z)
        self.supports = np.ones_like(nodes).astype(int)
        self.Ur = np.array([]).astype(int)
        # Connection vector of each bar, represented as ending node - starting node
        self.d = self.nodes[self.bars[:, 1], :] - self.nodes[self.bars[:, 0], :]
        # material constants
        self.A = A
        self.rho = rho
        self.E = E
        # second moment of area
        self.I = A ** 2 / 12

        # Compute length mass and orientation of each element
        self._computeLengths(self.d)
        self._computeOrientations(self.d)
        self._computeMass()

    # geometry
    def _computeLengths(self, d):
        """
        :param d: Connection vector of each bar, represented as ending node - starting node.

        Returns the length of each bar, calculated with the euclidian norm, as a matrix shaped like self.bars.
        """
        self.lengths = np.sqrt((d ** 2).sum(axis=1))  # lengths of all bars (euclidian norm)

    def _computeOrientations(self, d):
        """
        :param d: Connection vector of each bar, represented as starting node

        Returns the normalized orientation vector of each bar.

        example:
        d = np.array([1, 1]) -> np.array([1/np.sqrt(2), 1/np.sqrt(2)] '=' [cos(phi),sin(phi)]
        """
        self.orientations = d.T / self.lengths  # orientation vector of bars (transpose for matching dimensions)

    def _computeMass(self):
        """
        Computes mass of each element of thge truss
        """
        self.mass = self.lengths * self.A * self.rho

    def addSupport(self, node, x, y, z):
        """
        node:   id of node where support should be added 
        x:      displacement in x
        y:      displacement in y
        z:      displacement in z
            --> 0: free, 1: blocked
        """
        self.supports[node, 0] = x
        self.supports[node, 1] = y
        self.supports[node, 2] = z

        # appends a 0 per suppressed degree of freedom # TODO can be optimized
        if x == 0:
            self.Ur = np.append(self.Ur, 0)
        if y == 0:
            self.Ur = np.append(self.Ur, 0)
        if z == 0:
            self.Ur = np.append(self.Ur, 0)

    def addExternalForce(self, node, x, y, z):
        """
        Adds an external force. If there already is an external force at the given node,
        it will add one on top of it.
        node:   id of node where support should be added
        x:      x-component 
        y:      y-component 
        z:      z-component 
        """
        self.F[node][0] = self.F[node][0] + x
        self.F[node][1] = self.F[node][1] + y
        self.F[node][2] = self.F[node][2] + z

        self.force_points.append(node)  # save attack point

    def addExternalForces(self, forces):
        """
        Adds multiple external forces. If there already is an external force at a node,
        it will add one on top of it.
        forces: matrix of forces
        """
        self.F = self.F + forces

    def reset(self):
        """
        removes all external forces of the truss
        """
        self.F = np.zeros_like(self.nodes)
        self.force_points = []
