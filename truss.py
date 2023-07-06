import numpy as np


class Truss:
    """
    This class represents a simple truss that can be passed to the fem class to perform a truss analysis
    """

    def __init__(self, nodes, bars, A, rho, E):
        self.nodes = np.array(nodes).astype(float)
        self.bars = np.array(bars)
        self.F = np.zeros_like(nodes)  # external forces
        self.force_points = []  # indices of points, where external forces attack
        self.supports = np.ones_like(nodes).astype(
            int)  # degrees of freedom (1 = movable, 0 = fixed), 3 times for (x,y,z)
        self.Ur = np.array([]).astype(int)

        # material
        self.A = A
        self.rho = rho
        self.E = E
        self.I = A ** 2 / 12
        self.d = self.nodes[self.bars[:, 1], :] - self.nodes[self.bars[:, 0], :]  # direction vector for each bar
        self._computeLengths()
        self._computeOrientations()
        self._computeMass()

    # geometry
    def _computeLengths(self):
        """
        Computes the length of each bar and saves it in the attribute lengths which is shaped like bars
        """
        self.lengths = np.sqrt((self.d ** 2).sum(axis=1))  # length of bars (euclidian norm)

    def _computeOrientations(self):
        """
        Computes the orientation of each bar and saves it in the attribute orientations, it is represented by angles.
        """
        # print(np.sum(self.A ** 0.5 * self.lengths) / 2)
        self.orientations = self.d.T / self.lengths  # orientation vector of bars (transpose for matching dimensions)
    # physics
    def _computeMass(self):
        """
        Computes the mass of each bar and saves it in the attribute mass
        """
        self.mass = np.multiply(self.lengths, self.A) * self.rho

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

        # appends a 0 per suppressed degree of freedom
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
        Removes all external forces
        """
        self.F = np.zeros_like(self.nodes)
        self.force_points = []
