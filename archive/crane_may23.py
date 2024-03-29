# The pyramid crane doesnt work in this version
# This is the version we created on May 18

from truss import Truss
import numpy as np


class crane(Truss):
    """
    Special truss that represents a crane
    """

    def __init__(self, variant=0, height=10, length=10, hs=1, ls=1, A=25e-4, rho=7850, E=210e9):
        """
        :param height:
        :param length:
        :param hs:
        :param ls:
        """

        nodes = []
        bars = []

        self.height = height
        self.length = length
        self.hs = hs
        self.ls = ls

        self.nST = int(height / hs)  # Number of segments of the Tower
        self.nSA = int(length / ls)  # Number of segments of the Ausleger
        self.nSGA = int(length / (2 * ls))  # Number of segments of the Gegenausleger

        # indices of bars on a certain side of the crane
        self.x_side = []
        self.y_side = []

        if variant == 1:
            print("Variante 1: Pyramiden")
            self.tower_pyramid(nodes, bars)
            offsetT = self.cur_offset(nodes)
            self.gegenausleger_pyramid(nodes, bars, offsetT)
            offsetTG = self.cur_offset(nodes)
            self.ausleger_pyramid(nodes, bars, offsetT, offsetTG)
        else:
            print("Default Kran")
            self.tower(nodes, bars)
            offsetT = self.cur_offset(nodes)
            self.gegenausleger(nodes, bars, offsetT)
            offsetTG = self.cur_offset(nodes)
            self.ausleger(nodes, bars, offsetT, offsetTG)

        # convert python list to np.array
        self.nodes = np.array(nodes).astype(float)
        self.bars = np.array(bars)

        self.x_side = np.array(self.x_side).astype(int)
        self.y_side = np.array(self.y_side).astype(int)

        # Lager
        self.supports = np.ones_like(self.nodes).astype(int)
        self.Ur = np.array([]).astype('int')
        for i in range(4):
            self.addSupport(i, 0, 0, 0)

        # Externe Kräfte
        self.F = np.zeros_like(self.nodes)
        # for i in range(-1, -5, -1):
        #     self.addExternalForce(i, 0, 0, -5e4)  
        # self.addExternalForce(-3, 0, 0, -1e9)
        # self.addExternalForce(-4, 0, 0, -1e9)
        # Material
        self.A = A
        self.rho = rho
        self.E = E

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
            self.x_side.append(len(bars) - 1)
            bars.append([4 * i + 1, 4 * i + 3])  # RT -> RB
        self.x_side.remove(self.x_side[-1])
        self.x_side.remove(self.x_side[-1])

        # z-Richtung
        for i in range(self.nST - 1):
            bars.append([4 * i, 4 * i + 4])  # LT
            self.y_side.append(len(bars) - 1)
            self.x_side.append(len(bars) - 1)
            bars.append([4 * i + 1, 4 * i + 5])  # RT
            self.y_side.append(len(bars) - 1)
            bars.append([4 * i + 2, 4 * i + 6])  # LB
            self.x_side.append(len(bars) - 1)
            bars.append([4 * i + 3, 4 * i + 7])  # RB
        self.x_side.remove(self.x_side[-1])
        self.x_side.remove(self.x_side[-1])

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
            self.x_side.append(len(bars) - 1)
            # bars.append([4 * i, 4 * i + 6])  # LT -> LB+1
        self.x_side.remove(self.x_side[-1])

    def gegenausleger(self, nodes, bars, offsetT):
        # Gegenausleger erstellen

        # Nodes des Gegenauslegers in negative x Richtung, x Koordinaten mit +1, weil der Turm auf x= 0 bix x=1 steht
        for i in range(1, self.nSGA + 1):
            nodes.append([-(self.hs + i * self.ls) + self.ls, 0, (self.nST - 1) * self.hs])  # Left Top
            nodes.append([-(self.hs + i * self.ls) + self.ls, self.hs, (self.nST - 1) * self.hs])  # Right Top
            nodes.append([-(self.hs + i * self.ls) + self.ls, 0, (self.nST - 2) * self.hs])  # Left Bottom
            nodes.append([-(self.hs + i * self.ls) + self.ls, self.hs, (self.nST - 2) * self.hs])  # Right Bottom

        # Bars des Gegenausleger
        bars.append([offsetT - 2, offsetT + 1])  # hinten oben
        bars.append([offsetT - 4, offsetT])  # vorne oben
        bars.append([offsetT - 6, offsetT + 3])  # hinten unten
        bars.append([offsetT - 8, offsetT + 2])  # vorne unten

        # z-Richtung (LT für Left Top usw.)
        for i in range(self.nSGA):
            bars.append([4 * i + offsetT, 4 * i + offsetT + 1])  # LT -> RT ?
            bars.append([4 * i + 2 + offsetT, 4 * i + 3 + offsetT])  # LB -> RB
            bars.append([4 * i + offsetT, 4 * i + 2 + offsetT])  # LT -> LB
            bars.append([4 * i + 1 + offsetT, 4 * i + 3 + offsetT])  # RT -> RB
        self.x_side.append(len(bars) - 1)
        self.x_side.append(len(bars) - 2)
        self.x_side.append(len(bars) - 3)
        self.x_side.append(len(bars) - 4)

        # x- und y-Richtung
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
        # Ausleger erstellen

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
            if self.nSGA % 2 != 0 and i == self.nSA / 2 - 1:  # funktioniert nicht immer... für 8
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
            nodes.append([0, 0, i * self.hs])  # Left Top
            nodes.append([self.hs, 0, i * self.hs])  # Right Top
            nodes.append([0, self.hs, i * self.hs])  # Left Bottom
            nodes.append([self.hs, self.hs, i * self.hs])  # Right Bottom
        nodes.append([self.hs / 2, self.hs / 2, self.height])

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
            bars.append([4 * i + 1, 4 * i + 5])  # RT
            bars.append([4 * i + 2, 4 * i + 6])  # LB
            bars.append([4 * i + 3, 4 * i + 7])  # RB

        # Kreuzstreben (+1 jeweils für nächste Stufe)
        for i in range(self.nST - 1):
            bars.append([4 * i, 4 * i + 5])  # LT -> RT+1
            # bars.append([4 * i + 1, 4 * i + 4])  # RT -> RT+1
            bars.append([4 * i + 3, 4 * i + 6])  # RB -> LB+1
            # bars.append([4 * i + 2, 4 * i + 7])  # LB -> RB+1
            bars.append([4 * i + 1, 4 * i + 7])  # RT -> RB+1
            # bars.append([4 * i + 3, 4 * i + 5])  # RB -> RT+1
            bars.append([4 * i + 2, 4 * i + 4])  # LB -> LT+1
            # bars.append([4 * i, 4 * i + 6])  # LT -> LB+1

        # aller oberste Spitze
        offsetTO = len(nodes)
        bars.append([offsetTO - 1, offsetTO - 2])
        bars.append([offsetTO - 1, offsetTO - 3])
        bars.append([offsetTO - 1, offsetTO - 4])
        bars.append([offsetTO - 1, offsetTO - 5])

    def gegenausleger_pyramid(self, nodes, bars, offsetT):

        # Nodes des Gegenauslegers in negative x Richtung, x Koordinaten mit +1, weil der Turm auf x= 0 bix x=1 steht
        for i in range(1, self.nSGA + 1):  # braucht nur noch die ursprüungliche bottoms
            nodes.append([-(self.hs + i * self.ls) + self.ls, 0, (self.nST - 1) * self.hs])  # Left aself.lso y=0
            nodes.append(
                [-(self.hs + i * self.ls) + self.ls, self.hs, (self.nST - 1) * self.hs])  # Right aself.lso y=self.hs
            nodes.append([-(0.5 * self.hs + i * self.ls) + self.ls, self.hs / 2,
                          self.height])  # nodes der Spitzen --> gleiches problem mit doppelter Höhe??

        # Bars des Gegenausleger
        # sonderfall erste pyramide
        bars.append([offsetT, offsetT - 5])
        bars.append([offsetT + 1, offsetT - 3])
        bars.append([offsetT + 2, offsetT + 1])
        bars.append([offsetT + 2, offsetT])
        bars.append([offsetT + 2, offsetT - 5])
        bars.append([offsetT + 2, offsetT - 3])
        bars.append([offsetT + 2, offsetT - 1])

        # x- und y-Richtung (LT für Left Top usw.)
        for i in range(self.nSGA - 1):
            bars.append([offsetT + 3 * i, offsetT + 3 + 3 * i])
            bars.append([offsetT + 1 + 3 * i, offsetT + 1 + 3 + 3 * i])
            bars.append([offsetT + 3 * i, offsetT + 1 + 3 * i])
        offsetGT = len(nodes)
        bars.append([offsetGT - 2, offsetGT - 3])  # aller letzter vorne hinten Strich

        # Pyramiden
        for i in range(self.nSGA - 1):  # hier ab zweite pyramide
            bars.append([offsetT + 5 + 3 * i, offsetT + 5 + 3 * i - 1])
            bars.append([offsetT + 5 + 3 * i, offsetT + 5 + 3 * i - 2])
            bars.append([offsetT + 5 + 3 * i, offsetT + 5 + 3 * i - 4])
            bars.append([offsetT + 5 + 3 * i, offsetT + 5 + 3 * i - 5])

        # Linie oben
        for i in range(self.nSGA - 1):
            bars.append([offsetT + 2 + 3 * i, offsetT + 2 + 3 * i + 3])

    def ausleger_pyramid(self, nodes, bars, offsetT, offsetTG):

        # Nodes des Auslegers
        for i in range(1, self.nSA +1 ):  # braucht nur noch die ursprüungliche bottoms
            nodes.append([(self.hs + i * self.ls), 0, (self.nST - 1) * self.hs])  # Left bottem  also y=0
            nodes.append(
                [(self.hs + i * self.ls), self.hs, (self.nST - 1) * self.hs])  # Right bottom aself.lso y=self.hs
            nodes.append([(0.5 * self.hs + i * self.ls), self.hs / 2, self.height])  # tops

        # Bars des Ausleger
        for i in range(self.nSA - 1):
            bars.append([offsetTG + 3 * i, offsetTG + 3 + 3 * i])
            bars.append([offsetTG + 1 + 3 * i, offsetTG + 1 + 3 + 3 * i])
            bars.append([offsetTG + 3 * i, offsetTG + 1 + 3 * i])
        offsetGTA = len(nodes)
        bars.append([offsetGTA - 2, offsetGTA - 3])  # aller letzter vorne hinten Strich

        # Pyramiden
        for i in range(self.nSA - 1):  # hier ab zweite pyramide
            bars.append([offsetTG + 5 + 3 * i, offsetTG + 5 + 3 * i - 1])
            bars.append([offsetTG + 5 + 3 * i, offsetTG + 5 + 3 * i - 2])
            bars.append([offsetTG + 5 + 3 * i, offsetTG + 5 + 3 * i - 4])
            bars.append([offsetTG + 5 + 3 * i, offsetTG + 5 + 3 * i - 5])

        # Linie oben
        for i in range(self.nSA - 1):
            bars.append([offsetTG + 2 + 3 * i, offsetTG + 2 + 3 * i + 3])

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
