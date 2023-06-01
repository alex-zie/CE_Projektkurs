from truss import Truss
import numpy as np


class crane(Truss):
    """
    Special truss that represents a crane
    """

    def __init__(self, variant, height, length, hs, ls, A, rho, E, p=True):
        """
        :param variant:
        :param height:
        :param length:
        :param hs:
        :param ls:
        """
        if hs > height:
            raise Exception("Height of segments cannot be greater than the height of the crane!")
        if ls > length:
            raise Exception("Length of segments cannot be greater than the length of the jib!")

        nodes = []
        bars = []

        self.height = height
        self.length = length
        self.hs = hs
        self.ls = ls

        self.nST = np.ceil(height / hs).astype('int')  # Number of segments of the Tower
        self.nSA = np.ceil(length / ls).astype('int')  # Number of segments of the Ausleger
        self.nSGA = np.ceil((length / 2) / ls).astype('int')  # Number of segments of the Gegenausleger

        # indices of bars on a certain side of the crane
        self.x_side = []
        self.y_side = []

        if (variant == 1):
            if p: print(
                "Creating crane with cuboid tower with " + str(self.nST) + " segments and pyramidal jib with " + str(
                    self.nSA + self.nSGA) + " segments.")
            self.tower_pyramid(nodes, bars)
            offsetT = self.cur_offset(nodes)
            self.gegenausleger_pyramid(nodes, bars, offsetT)
            offsetTG = self.cur_offset(nodes)
            self.ausleger_pyramid(nodes, bars, offsetT, offsetTG)
        else:
            # print("Default Kran")
            self.tower(nodes, bars)
            offsetT = self.cur_offset(nodes)
            #self.gegenausleger(nodes, bars, offsetT)
            #offsetTG = self.cur_offset(nodes)
            #self.ausleger(nodes, bars, offsetT, offsetTG)

        # convert python list to np.array
        self.nodes = np.array(nodes).astype(float)
        self.bars = np.array(bars)

        self.x_side = np.array(self.x_side).astype(int)
        self.y_side = np.array(self.y_side).astype(int)
        # super().setLateralBars(self.x_side, self.y_side)

        # Lager
        self.supports = np.ones_like(self.nodes).astype(int)
        self.Ur = np.array([]).astype('int')
        for i in range(4):
            self.addSupport(i, 0, 0, 0)

        # Externe Kräfte
        self.F = np.zeros_like(self.nodes)

        # Material
        self.A = A
        self.rho = rho
        self.E = E
        self.I = A ** 2 / 12

        self._computeLengths()
        self._computeOrientations()
        self._computeMass()

    def tower(self, nodes, bars):
        # Turm erstellen
        # for 10m high tower with segement size hs we need at least
        # height / hs segements. For now i would just ignore the last segement
        # size != hs

        # Nodes des Turms
        for i in range(self.nST):
            nodes.append([0, 0, i * self.hs])  # Left Top
            nodes.append([self.hs, 0, i * self.hs])  # Right Top
            nodes.append([0, self.hs, i * self.hs])  # Left Bottom
            nodes.append([self.hs, self.hs, i * self.hs])  # Right Bottom

        # Bars des Turms
        # x- und y-Richtung (LT für Left Top usw.)
        for i in range(self.nST):
            bars.append([4 * i, 4 * i + 1])  # LT -> RT
            self.y_side.append(len(bars) - 1)
            bars.append([4 * i + 2, 4 * i + 3])  # LB -> RB
            bars.append([4 * i, 4 * i + 2])  # LT -> LB
            bars.append([4 * i + 1, 4 * i + 3])  # RT -> RB

        # z-Richtung
        for i in range(self.nST - 1):
            bars.append([4 * i, 4 * i + 4])  # LT
            self.y_side.append(len(bars) - 1)
            bars.append([4 * i + 1, 4 * i + 5])  # RT
            self.y_side.append(len(bars) - 1)
            bars.append([4 * i + 2, 4 * i + 6])  # LB
            bars.append([4 * i + 3, 4 * i + 7])  # RB

        # Kreuzstreben (+1 jeweils für nächste Stufe)
        for i in range(self.nST - 1):
            bars.append([4 * i, 4 * i + 5])  # LT -> RT+1
            self.y_side.append(len(bars) - 1)
            # bars.append([4 * i + 1, 4 * i + 4])  # RT -> RT+1
            bars.append([4 * i + 3, 4 * i + 6])  # RB -> LB+1
            # bars.append([4 * i + 2, 4 * i + 7])  # LB -> RB+1
            bars.append([4 * i + 1, 4 * i + 7])  # RT -> RB+1
            # bars.append([4 * i + 3, 4 * i + 5])  # RB -> RT+1
            bars.append([4 * i + 2, 4 * i + 4])  # LB -> LT+1
            # bars.append([4 * i, 4 * i + 6])  # LT -> LB+1

    def gegenausleger(self, nodes, bars, offsetT):
        # Gegenausleger erstellen

        # Nodes des Gegenauslegers in negative x Richtung, x Koordinaten mit +1, weil der Turm auf x= 0 bix x=1 steht
        for i in range(1, self.nSGA + 1):
            nodes.append([-(self.hs + i * self.ls) + 1, 0, (self.nST - 1) * self.hs])  # Left Top
            nodes.append([-(self.hs + i * self.ls) + 1, self.hs, (self.nST - 1) * self.hs])  # Right Top
            nodes.append([-(self.hs + i * self.ls) + 1, 0, (self.nST - 2) * self.hs])  # Left Bottom
            nodes.append([-(self.hs + i * self.ls) + 1, self.hs, (self.nST - 2) * self.hs])  # Right Bottom

        # Bars des Gegenausleger
        bars.append([offsetT - 2, offsetT + 1])  # hinten oben
        bars.append([offsetT - 4, offsetT])  # vorne oben
        bars.append([offsetT - 6, offsetT + 3])  # hinten unten
        bars.append([offsetT - 8, offsetT + 2])  # vorne unten

        # x- und y-Richtung (LT für Left Top usw.)
        for i in range(self.nSGA):
            bars.append([4 * i + offsetT, 4 * i + offsetT + 1])  # LT -> RT ?
            bars.append([4 * i + 2 + offsetT, 4 * i + 3 + offsetT])  # LB -> RB
            bars.append([4 * i + offsetT, 4 * i + 2 + offsetT])  # LT -> LB
            bars.append([4 * i + 1 + offsetT, 4 * i + 3 + offsetT])  # RT -> RB

        # z-Richtung
        for i in range(self.nSGA - 1):
            bars.append([4 * i + offsetT, 4 * i + 4 + offsetT])  # LT
            bars.append([4 * i + 1 + offsetT, 4 * i + 5 + offsetT])  # RT
            bars.append([4 * i + 2 + offsetT, 4 * i + 6 + offsetT])  # LB
            bars.append([4 * i + 3 + offsetT, 4 * i + 7 + offsetT])  # RB

        # Kreuzstreben vorne
        bars.append([offsetT - 8, offsetT])
        for i in range(self.nSGA // 2):
            bars.append([offsetT + 8 * i, offsetT + 6 + 8 * i])
            if self.nSGA % 2 != 0 and i == self.nSGA // 2:  # TODO checken, ob es allgemein funktioniert
                break
            bars.append([offsetT + 6 + 8 * i, offsetT + 8 + 8 * i])
        # Kreusstreben hinten
        bars.append([offsetT - 6, offsetT + 1])
        for i in range(self.nSGA // 2):
            bars.append([offsetT + 1 + 8 * i, offsetT + 7 + 8 * i])
            if self.nSGA % 2 != 0 and i == self.nSA // 2 - 1:
                break
            bars.append([offsetT + 7 + 8 * i, offsetT + 9 + 8 * i])

    def ausleger(self, nodes, bars, offsetT, offsetTG):
        # Nodes des Ausleger in positive x Richtung
        for i in range(1, self.nSA + 1):
            nodes.append([self.hs + i * self.ls, 0, (self.nST - 1) * self.hs])  # Left Top
            nodes.append([self.hs + i * self.ls, self.hs, (self.nST - 1) * self.hs])  # Right Top
            nodes.append([self.hs + i * self.ls, 0, (self.nST - 2) * self.hs])  # Left Bottom
            nodes.append([self.hs + i * self.ls, self.hs, (self.nST - 2) * self.hs])  # Right Bottom

        # Bars des Ausleger
        bars.append([offsetT - 4, offsetTG])  # LB -> LT
        bars.append([offsetT - 1, offsetTG + 1])  # RB -> RT
        bars.append([offsetT - 8, offsetTG + 2])  # LB-1 -> LB
        bars.append([offsetT - 5, offsetTG + 3])  # RB-1 -> RB
        # x- und y-Richtung (LT für Left Top usw.)
        for i in range(self.nSA):
            bars.append([4 * i + offsetTG, 4 * i + offsetTG + 1])  # LT -> RT
            bars.append([4 * i + 2 + offsetTG, 4 * i + 3 + offsetTG])  # LB -> RB
            bars.append([4 * i + offsetTG, 4 * i + 2 + offsetTG])  # LT -> LB
            bars.append([4 * i + 1 + offsetTG, 4 * i + 3 + offsetTG])  # RT -> RB

        # z-Richtung
        for i in range(self.nSA - 1):
            bars.append([4 * i + offsetTG, 4 * i + 4 + offsetTG])  # LT
            bars.append([4 * i + 1 + offsetTG, 4 * i + 5 + offsetTG])  # RT
            bars.append([4 * i + 2 + offsetTG, 4 * i + 6 + offsetTG])  # LB
            bars.append([4 * i + 3 + offsetTG, 4 * i + 7 + offsetTG])  # RB

        # Kreuzstreben vorne
        bars.append([offsetT - 7, offsetTG])
        for i in range(self.nSA // 2):
            bars.append([offsetTG + 8 * i, offsetTG + 6 + 8 * i])
            if self.nSGA % 2 != 0 and i == self.nSA // 2 - 1:
                break
            bars.append([offsetTG + 6 + 8 * i, offsetTG + 8 + 8 * i])

        # Kreusstreben hinten
        bars.append([offsetT - 5, offsetTG + 1])
        for i in range(self.nSA // 2):
            bars.append([offsetTG + 1 + 8 * i, offsetTG + 7 + 8 * i])
            if self.nSGA % 2 != 0 and i == self.nSA / 2 - 1:
                break
            bars.append([offsetTG + 7 + 8 * i, offsetTG + 9 + 8 * i])

    def cur_offset(self, nodes):
        return len(nodes)

    def tower_pyramid(self, nodes, bars):
        # Turm erstellen, der nur height-1 hoch ist und dann oben eine Spitze als Pyramide hat

        # Nodes des Turms
        for i in range(self.nST):
            nodes.append([0, 0, i * self.height / self.nST])  # Left Top
            nodes.append([self.hs, 0, i * self.height / self.nST])  # Right Top
            nodes.append([0, self.hs, i * self.height / self.nST])  # Left Bottom
            nodes.append([self.hs, self.hs, i * self.height / self.nST])  # Right Bottom
        nodes.append([self.hs / 2, self.hs / 2, self.height])

        # Bars des Turms
        # x- und y-Richtung (LT für Left Top usw.)
        for i in range(self.nST):
            bars.append([4 * i, 4 * i + 1])  # LT -> RT
            self.selectYbar(bars)
            bars.append([4 * i + 2, 4 * i + 3])  # LB -> RB
            bars.append([4 * i, 4 * i + 2])  # LT -> LB
            bars.append([4 * i + 1, 4 * i + 3])  # RT -> RB
            self.selectXbar(bars)
        # z-Richtung
        for i in range(self.nST - 1):
            bars.append([4 * i, 4 * i + 4])  # LT
            self.selectYbar(bars)
            bars.append([4 * i + 1, 4 * i + 5])  # RT
            self.selectYbar(bars)
            self.selectXbar(bars)
            bars.append([4 * i + 2, 4 * i + 6])  # LB
            bars.append([4 * i + 3, 4 * i + 7])  # RB
            self.selectXbar(bars)

        # Kreuzstreben (+1 jeweils für nächste Stufe)
        for i in range(self.nST - 1):
            bars.append([4 * i, 4 * i + 5])  # LT -> RT+1
            self.selectYbar(bars)
            # bars.append([4 * i + 1, 4 * i + 4])  # RT -> RT+1
            bars.append([4 * i + 3, 4 * i + 6])  # RB -> LB+1
            self.selectYbar(bars)
            # bars.append([4 * i + 2, 4 * i + 7])  # LB -> RB+1
            bars.append([4 * i + 1, 4 * i + 7])  # RT -> RB+1
            self.selectXbar(bars)
            # bars.append([4 * i + 3, 4 * i + 5])  # RB -> RT+1
            bars.append([4 * i + 2, 4 * i + 4])  # LB -> LT+1
            self.selectXbar(bars)
            # bars.append([4 * i, 4 * i + 6])  # LT -> LB+1

        #aller oberste Spitze
        offsetTO = len(nodes)
        bars.append([offsetTO - 1, offsetTO - 2])
        bars.append([offsetTO - 1, offsetTO - 3])
        bars.append([offsetTO - 1, offsetTO - 4])
        self.selectYbar(bars)
        bars.append([offsetTO - 1, offsetTO - 5])
        self.selectYbar(bars)

    def gegenausleger_pyramid(self, nodes, bars, offsetT):

        # Nodes des Gegenauslegers in negative x Richtung, x Koordinaten mit +1, weil der Turm auf x= 0 bix x=1 steht
        for i in range(1, self.nSGA + 1):  # braucht nur noch die ursprüungliche bottoms
            nodes.append([-(self.hs + i * self.ls) + self.ls, 0, (self.nST-1) * (self.height / self.nST)])  # Left aself.lso y=0
            nodes.append([-(self.hs + i * self.ls) + self.ls, self.hs, (self.nST-1) * (self.height / self.nST)])  # Right aself.lso y=self.hs
            nodes.append([-(0.5 * self.hs + i * self.ls) + self.ls, self.hs / 2,(self.nST) * (self.height / self.nST)])  # nodes der Spitzen --> gleiches problem mit doppelter Höhe??

        # Bars des Gegenausleger

        #sonderfall erste pyramide
        bars.append([offsetT, offsetT - 5])
        self.selectYbar(bars)
        bars.append([offsetT + 1, offsetT - 3])
        bars.append([offsetT + 2, offsetT + 1])
        bars.append([offsetT + 2, offsetT])
        self.selectYbar(bars)
        bars.append([offsetT + 2, offsetT - 5])
        self.selectYbar(bars)
        bars.append([offsetT + 2, offsetT - 3])
        bars.append([offsetT + 2, offsetT - 1])
        self.selectYbar(bars)

        # x- und y-Richtung (LT für Left Top usw.)
        for i in range(self.nSGA - 1):
            bars.append([offsetT + 3 * i, offsetT + 3 + 3 * i])
            self.selectYbar(bars)
            bars.append([offsetT + 1 + 3 * i, offsetT + 1 + 3 + 3 * i])
            bars.append([offsetT + 3 * i, offsetT + 1 + 3 * i])
        offsetGT = len(nodes)
        bars.append([offsetGT - 2, offsetGT - 3])  # aller letzter vorne hinten Strich

        # Pyramiden
        for i in range(self.nSGA - 1):  # hier ab zweite pyramide
            bars.append([offsetT + 5 + 3 * i, offsetT + 5 + 3 * i - 1])
            bars.append([offsetT + 5 + 3 * i, offsetT + 5 + 3 * i - 2])
            self.selectYbar(bars)
            bars.append([offsetT + 5 + 3 * i, offsetT + 5 + 3 * i - 4])
            bars.append([offsetT + 5 + 3 * i, offsetT + 5 + 3 * i - 5])
            self.selectYbar(bars)

        # Linie oben
        for i in range(self.nSGA - 1):
            bars.append([offsetT + 2 + 3 * i, offsetT + 2 + 3 * i + 3])
            self.selectYbar(bars)

    def ausleger_pyramid(self, nodes, bars, offsetT, offsetTG):
        for i in range(1, self.nSA + 1):
            nodes.append([self.hs + i * self.ls, 0, (self.nST-1) * (self.height / self.nST)])  # Left Bottom
            nodes.append([self.hs + i * self.ls, self.hs, (self.nST-1) * (self.height / self.nST)])  # Right Bottom
            nodes.append([(self.hs/2 + i * self.ls), self.hs / 2, (self.nST) * (self.height / self.nST)])  # Top


        # x- und y-Richtung
        for i in range(self.nSA - 1):
            bars.append([offsetTG + i * 3, (offsetTG + 3) + i * 3])
            self.selectYbar(bars)
            bars.append([offsetTG + i * 3, (offsetTG + 1) + i * 3])
            bars.append([(offsetTG + 1) + i * 3, (offsetTG + 4) + i * 3])

        # Bottom nodes with top nodes
        tmp_lastbar1 = 0
        tmp_lastbar2 = 0
        for i in range(self.nSA - 1):
            bars.append([offsetTG + i * 3, (offsetTG + 5) + i * 3])
            self.selectYbar(bars)
            bars.append([offsetTG + 1 + i * 3, (offsetTG + 5) + i * 3])
            bars.append([(offsetTG + 5) + i * 3, offsetTG + 3 + i * 3])
            self.selectYbar(bars)
            tmp_lastbar1 = len(bars) - 1
            bars.append([(offsetTG + 5) + i * 3, offsetTG + 4 + i * 3])
            tmp_lastbar2 = len(bars) - 1
        self.x_side.append(tmp_lastbar1)
        self.x_side.append(tmp_lastbar2)

        # Top Row
        for i in range(self.nSA - 1):
            bars.append([(offsetTG + 2) + i * 3, (offsetTG + 5) + i * 3])
            self.selectYbar(bars)

        # Extra bars
        offsetTO = len(nodes)  # offset after all the nodes

        # Last bar ate the end of the crane
        bars.append([offsetTO - 2, offsetTO - 3])
        self.selectXbar(bars)

        # Top of the Tower with first Node Ausleger
        bars.append([offsetT - 1, offsetTG + 2])
        self.selectYbar(bars)

        # Tower with the base of the Ausleger
        bars.append([offsetT - 4, offsetTG])
        self.selectYbar(bars)
        bars.append([offsetT - 2, offsetTG + 1])

        # Tower with the first Top Node in Ausleger
        bars.append([offsetT - 4, offsetTG + 2])
        self.selectYbar(bars)
        bars.append([offsetT - 2, offsetTG + 2])

        # First Top Node to base Ausleger2
        bars.append([offsetTG + 2, offsetTG + 1])
        bars.append([offsetTG + 2, offsetTG])
        self.selectYbar(bars)

    def selectYbar(self, bars):
        """
        Select the last bar from the bar array and add this to another array to select the Y bars
        """
        self.y_side.append(len(bars) - 1)

    def selectXbar(self, bars):
        """
        Select the last bar from the bar array and add this to another array to select the X bars
        """
        self.x_side.append(len(bars) - 1)

# from truss import Truss
# import numpy as np


# class crane(Truss):
#     """
#     Special truss that represents a crane
#     """

#     def __init__(self, height, length, hs, ls, A, rho, E):  # maybe we can add some unique function e.g. crane_variant_1 to generalize
#         """
#         :param height:
#         :param length:
#         :param hs:
#         :param ls:
#         """

#         nodes = []
#         bars = []


#         #Turm erstellen
#         # for 10m high tower with segement size hs we need at least
#         # height / hs segements. For now i would just ignore the last segement
#         # size != hs
#         numSegmentsT = np.ceil(height / hs).astype('int')
#         numSegmentsA = np.ceil(length / ls).astype('int')
#         numSegmentsT = np.ceil(height / hs).astype('int')
#         numSegmentsGA = 1 #np.ceil((length / 2) / (ls / 2)).astype('int')


#         #Nodes des Turms
#         for i in range(numSegmentsT):
#             nodes.append([0, 0, i * hs])  # Left Top
#             nodes.append([hs, 0, i * hs])  # Right Top
#             nodes.append([0, hs, i * hs])  # Left Bottom
#             nodes.append([hs, hs, i * hs])  # Right Bottom
#         offset = len(nodes)


#         #Bars des Turms
#         # x- und y-Richtung (LT für Left Top usw.)
#         for i in range(numSegmentsT):
#             bars.append([4 * i, 4 * i + 1])  # LT -> RT
#             bars.append([4 * i + 2, 4 * i + 3])  # LB -> RB
#             bars.append([4 * i, 4 * i + 2])  # LT -> LB
#             bars.append([4 * i + 1, 4 * i + 3])  # RT -> RB

#         # z-Richtung
#         for i in range(numSegmentsT - 1):
#             bars.append([4 * i, 4 * i + 4])  # LT
#             bars.append([4 * i + 1, 4 * i + 5])  # RT
#             bars.append([4 * i + 2, 4 * i + 6])  # LB
#             bars.append([4 * i + 3, 4 * i + 7])  # RB

#         # Kreuzstreben (+1 jeweils für nächste Stufe)
#         for i in range(numSegmentsT - 1):
#             bars.append([4 * i, 4 * i + 5])  # LT -> RT+1
#             # bars.append([4 * i + 1, 4 * i + 4])  # RT -> RT+1
#             bars.append([4 * i + 3, 4 * i + 6])  # RB -> LB+1
#             # bars.append([4 * i + 2, 4 * i + 7])  # LB -> RB+1
#             bars.append([4 * i + 1, 4 * i + 7])  # RT -> RB+1
#             # bars.append([4 * i + 3, 4 * i + 5])  # RB -> RT+1
#             bars.append([4 * i + 2, 4 * i + 4])  # LB -> LT+1
#             # bars.append([4 * i, 4 * i + 6])  # LT -> LB+1


#         #Gegenausleger erstellen
#         #Nodes des Gegenauslegers in negative x Richtung, x Koordinaten mit +1, weil der Turm auf x= 0 bix x=1 steht
#         for i in range(1, numSegmentsGA + 1):
#             nodes.append([-(hs + i * ls)+1, 0, (numSegmentsT - 1) * hs])  # Left Top
#             nodes.append([-(hs + i * ls)+1, hs, (numSegmentsT - 1) * hs])  # Right Top
#             nodes.append([-(hs + i * ls)+1, 0, (numSegmentsT - 2) * hs])  # Left Bottom
#             nodes.append([-(hs + i * ls)+1, hs, (numSegmentsT - 2) * hs])  # Right Bottom
#         offsetTG = len(nodes) # Turm und Gegenausleger


#             #Bars des Gegenausleger
#         bars.append([offset-2, offset+1])  # hinten oben
#         bars.append([offset-4, offset])    # vorne oben
#         bars.append([offset-6, offset+3])  # hinten unten
#         bars.append([offset-8, offset+2])  # vorne unten

#         # x- und y-Richtung (LT für Left Top usw.)
#         for i in range(numSegmentsGA):
#             bars.append([4 * i + offset, 4 * i + offset + 1])  # LT -> RT ?
#             bars.append([4 * i + 2 + offset, 4 * i + 3 + offset])  # LB -> RB
#             bars.append([4 * i + offset, 4 * i + 2 + offset])  # LT -> LB
#             bars.append([4 * i + 1 + offset, 4 * i + 3 + offset])  # RT -> RB

#         # z-Richtung
#         for i in range(numSegmentsGA-1):
#             bars.append([4 * i + offset, 4 * i + 4 + offset])  # LT
#             bars.append([4 * i + 1 + offset, 4 * i + 5 + offset])  # RT
#             bars.append([4 * i + 2 + offset, 4 * i + 6 + offset])  # LB
#             bars.append([4 * i + 3 + offset, 4 * i + 7 + offset])  # RB

#         #Kreuzstreben vorne
#         bars.append([offset-8, offset])
#         for i in range(numSegmentsGA//2):
#             bars.append([offset + 8 * i, offset + 6 + 8 * i])
#             if numSegmentsGA % 2 != 0 and i == numSegmentsGA//2: #TODO checken, ob es allgemein funktioniert
#                 break
#             bars.append([ offset + 6 + 8 * i, offset + 8 + 8 * i])

#         #Kreusstreben hinten
#         bars.append([offset-6, offset+1])
#         for i in range(numSegmentsGA//2):
#             bars.append([offset +1 + 8 * i, offset + 7 + 8 * i])
#             if numSegmentsGA % 2 != 0 and i == numSegmentsA//2-1:
#                 break
#             bars.append([ offset + 7 + 8 * i, offset + 9 + 8 * i])


#     #Ausleger erstellen

#         #Nodes des Ausleger in positive x Richtung
#         for i in range(1, numSegmentsA + 1):
#             nodes.append([hs + i * ls, 0, (numSegmentsT - 1) * hs])  # Left Top
#             nodes.append([hs + i * ls, hs, (numSegmentsT - 1) * hs])  # Right Top
#             nodes.append([hs + i * ls, 0, (numSegmentsT - 2) * hs])  # Left Bottom
#             nodes.append([hs + i * ls, hs, (numSegmentsT - 2) * hs])  # Right Bottom


#         #Bars des Ausleger
#         bars.append([offset - 4, offsetTG])  # LB -> LT
#         bars.append([offset - 1, offsetTG + 1])  # RB -> RT
#         bars.append([offset - 8, offsetTG + 2])  # LB-1 -> LB
#         bars.append([offset - 5, offsetTG + 3])  # RB-1 -> RB
#         # x- und y-Richtung (LT für Left Top usw.)
#         for i in range(numSegmentsA):
#             bars.append([4 * i + offsetTG, 4 * i + offsetTG + 1])  # LT -> RT
#             bars.append([4 * i + 2 + offsetTG, 4 * i + 3 + offsetTG])  # LB -> RB
#             bars.append([4 * i + offsetTG, 4 * i + 2 + offsetTG])  # LT -> LB
#             bars.append([4 * i + 1 + offsetTG, 4 * i + 3 + offsetTG])  # RT -> RB

#         # z-Richtung
#         for i in range(numSegmentsA - 1):
#             bars.append([4 * i + offsetTG, 4 * i + 4 + offsetTG])  # LT
#             bars.append([4 * i + 1 + offsetTG, 4 * i + 5 + offsetTG])  # RT
#             bars.append([4 * i + 2 + offsetTG, 4 * i + 6 + offsetTG])  # LB
#             bars.append([4 * i + 3 + offsetTG, 4 * i + 7 + offsetTG])  # RB

#         #Kreuzstreben vorne
#         bars.append([offset-7, offsetTG])
#         for i in range (numSegmentsA//2):
#             bars.append([offsetTG + 8 * i, offsetTG + 6 + 8 * i])
#             if numSegmentsGA % 2 != 0 and i == numSegmentsA//2-1:
#                 break
#             bars.append([ offsetTG + 6 + 8 * i, offsetTG + 8 + 8 * i])

#         #Kreusstreben hinten
#         bars.append([offset-5, offsetTG+1])
#         for i in range (numSegmentsA//2):
#             bars.append([offsetTG +1 + 8 * i, offsetTG + 7 + 8 * i])
#             if numSegmentsGA % 2 != 0 and i == numSegmentsA/2-1:
#                 break
#             bars.append([ offsetTG + 7 + 8 * i, offsetTG + 9 + 8 * i])


#         # convert python list to np.array
#         self.nodes = np.array(nodes).astype(float)
#         self.bars = np.array(bars)

#         # Lager
#         self.supports = np.ones_like(self.nodes).astype(int)
#         self.Ur = np.array([]).astype('int')
#         for i in range(4):
#             self.addSupport(i, 0, 0, 0)

#         # Externe Kräfte
#         self.F = np.zeros_like(self.nodes)
#         # for i in range(-1, -5, -1):
#         #     self.addExternalForce(i, 0, 0, -125e3)
#         # self.addExternalForce(-3, 0, 0, -1e9)
#         # self.addExternalForce(-4, 0, 0, -1e9)
#         # Material
#         self.A = A
#         self.rho = rho
#         self.E = E

#         self._computeLengths()
#         self._computeOrientations()
#         self._computeMass()
