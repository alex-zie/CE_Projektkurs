from truss import Truss
import numpy as np


class crane(Truss):
    """
    Special truss that represents a crane
    """

    def __init__(self, height, length, hs, ls, A, rho,
                 E):  # maybe we can add some unique function e.g. crane_variant_1 to generalize
        """
        :param height:
        :param length:
        :param hs:
        :param ls:
        """
        nodes = []
        bars = []

        # for 10m high tower with segement size hs we need at least
        # height / hs segements. For now i would just ignore the last segement
        # size != hs
        numSegmentsA = np.ceil(length / ls).astype('int')
        numSegmentsT = np.ceil(height / hs).astype('int')
        for i in range(numSegmentsT):
            nodes.append([0, 0, i * hs])  # Left Top
            nodes.append([hs, 0, i * hs])  # Right Top
            nodes.append([0, hs, i * hs])  # Left Bottom
            nodes.append([hs, hs, i * hs])  # Right Bottom
        offset = len(nodes)

        # Erstelle die nodes des Auslegers in positive x Richtung
        for i in range(1, numSegmentsA + 1):
            nodes.append([hs + i * ls, 0, (numSegmentsT - 1) * hs])  # Left Top
            nodes.append([hs + i * ls, hs, (numSegmentsT - 1) * hs])  # Right Top
            nodes.append([hs + i * ls, 0, (numSegmentsT - 2) * hs])  # Left Bottom
            nodes.append([hs + i * ls, hs, (numSegmentsT - 2) * hs])  # Right Bottom

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
            # bars.append([4 * i + 1, 4 * i + 4])  # RT -> RT+1
            bars.append([4 * i + 3, 4 * i + 6])  # RB -> LB+1
            # bars.append([4 * i + 2, 4 * i + 7])  # LB -> RB+1
            bars.append([4 * i + 1, 4 * i + 7])  # RT -> RB+1
            # bars.append([4 * i + 3, 4 * i + 5])  # RB -> RT+1
            bars.append([4 * i + 2, 4 * i + 4])  # LB -> LT+1
            # bars.append([4 * i, 4 * i + 6])  # LT -> LB+1

        # Ausleger
        bars.append([offset - 4, offset])  # LB -> LT
        bars.append([offset - 1, offset + 1])  # RB -> RT
        bars.append([offset - 8, offset + 2])  # LB-1 -> LB
        bars.append([offset - 5, offset + 3])  # RB-1 -> RB
        # x- und y-Richtung (LT für Left Top usw.)
        for i in range(numSegmentsA):
            bars.append([4 * i + offset, 4 * i + offset + 1])  # LT -> RT
            bars.append([4 * i + 2 + offset, 4 * i + 3 + offset])  # LB -> RB
            bars.append([4 * i + offset, 4 * i + 2 + offset])  # LT -> LB
            bars.append([4 * i + 1 + offset, 4 * i + 3 + offset])  # RT -> RB

        # z-Richtung
        for i in range(numSegmentsA - 1):
            bars.append([4 * i + offset, 4 * i + 4 + offset])  # LT
            bars.append([4 * i + 1 + offset, 4 * i + 5 + offset])  # RT
            bars.append([4 * i + 2 + offset, 4 * i + 6 + offset])  # LB
            bars.append([4 * i + 3 + offset, 4 * i + 7 + offset])  # RB

        # convert python list to np.array
        self.nodes = np.array(nodes).astype(float)
        self.bars = np.array(bars)

        # Lager
        self.supports = np.ones_like(self.nodes).astype(int)
        self.Ur = np.array([]).astype('int')
        for i in range(4):
            self.addSupport(i, 0, 0, 0)

        # Externe Kräfte
        self.F = np.zeros_like(self.nodes)
        print(np.array(range(-1, -5, -1)))
        for i in range(-1, -5, -1):
            self.addExternalForce(i, 0, 0, -1e9)  
        # self.addExternalForce(-3, 0, 0, -1e9)
        # self.addExternalForce(-4, 0, 0, -1e9)
        # Material
        self.A = A
        self.rho = rho
        self.E = E

        self._computeLengths()
        self._computeOrientations()
        self._computeMass()
