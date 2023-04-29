import numpy as np

class Crane:

    def __init__(self):
        self.nodes = []
        self.bars = []

    def make_crane(self,height, length, hs, ls):
        """
        maybe we can add some unique function e.g. crane_variant_1 to generalize
        :param height:
        :param length:
        :param hs:
        :param ls:
        :return:
        """
        # for 10m high tower with segement size hs we need at least
        # height / hs segements. For now i would just ignore the last segement
        # size != hs
        numSegmentsA = np.ceil(length / ls).astype('int')
        numSegmentsT = np.ceil(height / hs).astype('int')
        for i in range(numSegmentsT):
            self.nodes.append([0, 0, i * hs])  # Left Top
            self.nodes.append([hs, 0, i * hs])  # Right Top
            self.nodes.append([0, hs, i * hs])  # Left Bottom
            self.nodes.append([hs, hs, i * hs])  # Right Bottom
        offset = len(self.nodes)
        # Erstelle die self.nodes des Auslegers in positive x Richtung
        for i in range(1, numSegmentsA + 1):
            self.nodes.append([hs + i * ls, 0, (numSegmentsT - 1) * hs])  # Left Top
            self.nodes.append([hs + i * ls, hs, (numSegmentsT - 1) * hs])  # Right Top
            self.nodes.append([hs + i * ls, 0, (numSegmentsT - 2) * hs])  # Left Bottom
            self.nodes.append([hs + i * ls, hs, (numSegmentsT - 2) * hs])  # Right Bottom

        # Turm
        # x- und y-Richtung (LT für Left Top usw.)
        for i in range(numSegmentsT):
            self.bars.append([4 * i, 4 * i + 1])  # LT -> RT
            self.bars.append([4 * i + 2, 4 * i + 3])  # LB -> RB
            self.bars.append([4 * i, 4 * i + 2])  # LT -> LB
            self.bars.append([4 * i + 1, 4 * i + 3])  # RT -> RB

        # z-Richtung
        for i in range(numSegmentsT - 1):
            self.bars.append([4 * i, 4 * i + 4])  # LT
            self.bars.append([4 * i + 1, 4 * i + 5])  # RT
            self.bars.append([4 * i + 2, 4 * i + 6])  # LB
            self.bars.append([4 * i + 3, 4 * i + 7])  # RB

        # Kreuzstreben (+1 jeweils für nächste Stufe)
        for i in range(numSegmentsT - 1):
            self.bars.append([4 * i, 4 * i + 5])  # LT -> RT+1
            # self.bars.append([4 * i + 1, 4 * i + 4])  # RT -> RT+1
            self.bars.append([4 * i + 3, 4 * i + 6])  # RB -> LB+1
            # self.bars.append([4 * i + 2, 4 * i + 7])  # LB -> RB+1
            self.bars.append([4 * i + 1, 4 * i + 7])  # RT -> RB+1
            # self.bars.append([4 * i + 3, 4 * i + 5])  # RB -> RT+1
            self.bars.append([4 * i + 2, 4 * i + 4])  # LB -> LT+1
            # self.bars.append([4 * i, 4 * i + 6])  # LT -> LB+1

        # Ausleger
        self.bars.append([offset - 2, offset])  # LB -> LT
        self.bars.append([offset - 1, offset + 1])  # RB -> RT
        self.bars.append([offset - 6, offset + 2])  # LB-1 -> LB
        self.bars.append([offset - 5, offset + 3])  # RB-1 -> RB
        # x- und y-Richtung (LT für Left Top usw.)
        for i in range(numSegmentsA):
            self.bars.append([4 * i + offset, 4 * i + offset + 1])  # LT -> RT
            self.bars.append([4 * i + 2 + offset, 4 * i + 3 + offset])  # LB -> RB
            self.bars.append([4 * i + offset, 4 * i + 2 + offset])  # LT -> LB
            self.bars.append([4 * i + 1 + offset, 4 * i + 3 + offset])  # RT -> RB

        # z-Richtung
        for i in range(numSegmentsA - 1):
            self.bars.append([4 * i + offset, 4 * i + 4 + offset])  # LT
            self.bars.append([4 * i + 1 + offset, 4 * i + 5 + offset])  # RT
            self.bars.append([4 * i + 2 + offset, 4 * i + 6 + offset])  # LB
            self.bars.append([4 * i + 3 + offset, 4 * i + 7 + offset])  # RB


    def getNB(self) -> np.ndarray:
        return np.array(self.nodes).astype(float) , np.array(self.bars)
