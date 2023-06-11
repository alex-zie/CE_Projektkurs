from truss import Truss
import numpy as np


class crane_1(Truss):
    """
    Special truss that represents a crane
    """

    def __init__(self, height, length, ls, A, rho, E, p=True):
        """
        :param variant:
        :param height:
        :param length:
        :param ls:
        """
        if ls > height:
            raise Exception("Height of segments cannot be greater than the height of the crane!")
        if ls > length:
            raise Exception("Length of segments cannot be greater than the length of the jib!")

        nodes = []
        bars = []

        self.height = height
        self.length = length
        self.ls = ls

        self.nST = np.ceil(height / ls).astype('int')  # Number of segments of the Tower
        self.nSA = np.ceil(length / ls).astype('int')  # Number of segments of the Ausleger
        self.nSGA = np.ceil((length / 2) / ls).astype('int')  # Number of segments of the Gegenausleger

        self.ls = np.min([self.height / self.nST, self.length / self.nSA])

        # indices of bars on a certain side of the crane
        self.x_negative_side = []  # Done
        self.x_positive_side = []  # Done
        self.y_negative_side = []
        self.y_positive_side = []  # Done

        if p:
            print("Creating crane with cuboid tower with " + str(self.nST) + " segments of length "+str(self.ls)+" and pyramidal jib with " + str(self.nSA) + " segments.")
        self.tower_pyramid(nodes, bars)
        offsetT = cur_offset(nodes)
        self.counterweight_nodes = [] # Knoten des Gegenauslegers zur Anbringung der Gegenlast
        self.gegenausleger_pyramid(nodes, bars, offsetT)
        offsetTG = cur_offset(nodes)
        self.ausleger_pyramid(nodes, bars, offsetT, offsetTG)
        self.tip_nodes = [-2, -3] # Knoten an der Spitze des Auslegers zur Lastanbringung

        # convert python list to np.array
        self.nodes = np.array(nodes).astype(float)
        self.bars = np.array(bars)

        self.x_negative_side = np.array(self.x_negative_side).astype(int)
        self.x_positive_side = np.array(self.x_positive_side).astype(int)
        self.y_positive_side = np.array(self.y_positive_side).astype(int)
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
    
    def tower_pyramid(self, nodes, bars):
        # Turm erstellen, der nur height-1 hoch ist und dann oben eine Spitze als Pyramide hat

        # Nodes des Turms
        for i in range(self.nST):
            nodes.append([0, 0, i * self.ls])  # Left Top
            nodes.append([self.ls, 0, i * self.ls])  # Right Top
            nodes.append([0, self.ls, i * self.ls])  # Left Bottom
            nodes.append([self.ls, self.ls, i * self.ls])  # Right Bottom
        nodes.append([self.ls / 2, self.ls / 2, self.ls * self.nST])

        # Bars des Turms
        # x- und y-Richtung (LT für Left Top usw.)
        for i in range(self.nST):
            bars.append([4 * i, 4 * i + 1])  # LT -> RT
            self.selectYPositiveBar(bars)
            bars.append([4 * i + 2, 4 * i + 3])  # LB -> RB
            self.selectYNegativeBar(bars)
            bars.append([4 * i, 4 * i + 2])  # LT -> LB
            self.selectXPositiveBar(bars)
            bars.append([4 * i + 1, 4 * i + 3])  # RT -> RB
            self.selectXNegativeBar(bars)
        # z-Richtung
        for i in range(self.nST - 1):
            bars.append([4 * i, 4 * i + 4])  # LT
            self.selectYPositiveBar(bars)
            self.selectXPositiveBar(bars)
            bars.append([4 * i + 1, 4 * i + 5])  # RT
            self.selectYPositiveBar(bars)
            self.selectXNegativeBar(bars)
            bars.append([4 * i + 2, 4 * i + 6])  # LB
            self.selectXPositiveBar(bars)
            self.selectYNegativeBar(bars)
            bars.append([4 * i + 3, 4 * i + 7])  # RB
            self.selectXNegativeBar(bars)
            self.selectYNegativeBar(bars)

        # Kreuzstreben (+1 jeweils für näclste Stufe)
        for i in range(self.nST - 1):
            bars.append([4 * i, 4 * i + 5])  # LT -> RT+1
            self.selectYPositiveBar(bars)
            self.selectYNegativeBar(bars)
            # bars.append([4 * i + 1, 4 * i + 4])  # RT -> RT+1
            bars.append([4 * i + 3, 4 * i + 6])  # RB -> LB+1
            self.selectYPositiveBar(bars)
            self.selectYNegativeBar(bars)
            # bars.append([4 * i + 2, 4 * i + 7])  # LB -> RB+1
            bars.append([4 * i + 1, 4 * i + 7])  # RT -> RB+1
            self.selectXNegativeBar(bars)
            self.selectXPositiveBar(bars)
            # bars.append([4 * i + 3, 4 * i + 5])  # RB -> RT+1
            bars.append([4 * i + 2, 4 * i + 4])  # LB -> LT+1
            self.selectXNegativeBar(bars)
            self.selectXPositiveBar(bars)
            # bars.append([4 * i, 4 * i + 6])  # LT -> LB+1

        #aller oberste Spitze
        offsetTO = len(nodes)
        bars.append([offsetTO - 1, offsetTO - 2])
        self.selectYNegativeBar(bars)
        bars.append([offsetTO - 1, offsetTO - 3])
        self.selectYNegativeBar(bars)
        bars.append([offsetTO - 1, offsetTO - 4])
        self.selectYPositiveBar(bars)
        bars.append([offsetTO - 1, offsetTO - 5])
        self.selectYPositiveBar(bars)
        
        #Hilfslinie
        bars.append([offsetTO - 3, offsetTO - 4])
        bars.append([offsetTO - 2, offsetTO - 5])

    def gegenausleger_pyramid(self, nodes, bars, offsetT):

        # Nodes des Gegenauslegers in negative x Richtung, x Koordinaten mit +1, weil der Turm auf x= 0 bix x=1 steht
        for i in range(1, self.nSGA + 1):  # braucht nur noch die ursprüngliche bottoms
            nodes.append([-(self.ls + i * self.ls) + self.ls, 0, (self.nST-1) * self.ls])  # Left aself.lso y=0
            nodes.append([-(self.ls + i * self.ls) + self.ls, self.ls, (self.nST-1) * self.ls])  # Right aself.lso y=self.ls
            self.counterweight_nodes.append(len(nodes)-2)
            self.counterweight_nodes.append(len(nodes)-1)
            nodes.append([-(0.5 * self.ls + i * self.ls) + self.ls, self.ls / 2,(self.nST) * self.ls])  # nodes der Spitzen --> gleiches problem mit doppelter Höhe??

        # Bars des Gegenausleger

        #sonderfall erste pyramide
        bars.append([offsetT, offsetT - 5])
        self.selectYPositiveBar(bars)
        bars.append([offsetT + 1, offsetT - 3])
        self.selectYNegativeBar(bars)
        bars.append([offsetT + 2, offsetT + 1])
        self.selectYNegativeBar(bars)
        bars.append([offsetT + 2, offsetT])
        self.selectYPositiveBar(bars)

        bars.append([offsetT + 2, offsetT - 5])
        self.selectYPositiveBar(bars)

        bars.append([offsetT + 2, offsetT - 3])
        self.selectYNegativeBar(bars)
        bars.append([offsetT + 2, offsetT - 1])
        self.selectYPositiveBar(bars)
        self.selectYNegativeBar(bars)

        # x- und y-Richtung (LT für Left Top usw.)
        for i in range(self.nSGA - 1):
            bars.append([offsetT + 3 * i, offsetT + 3 + 3 * i])
            self.selectYPositiveBar(bars)
            bars.append([offsetT + 1 + 3 * i, offsetT + 1 + 3 + 3 * i])
            self.selectYNegativeBar(bars)
            bars.append([offsetT + 3 * i, offsetT + 1 + 3 * i])
        offsetGT = len(nodes)
        bars.append([offsetGT - 2, offsetGT - 3])  # aller letzter vorne hinten Strich
        self.selectXPositiveBar(bars)

        #Hilfslinien unten
        bars.append([offsetT + 1, offsetT -5])        #sonderfall erster Strich

        for i in range(self.nSGA - 1):
            bars.append([offsetT + 3 * i, offsetT + 4 + 3 * i])

        # Pyramiden
        tmp_lastbar1 = 0
        tmp_lastbar2 = 0
        for i in range(self.nSGA - 1):  # hier ab zweite pyramide
            bars.append([offsetT + 5 + 3 * i, offsetT + 5 + 3 * i - 1])
            self.selectYNegativeBar(bars)
            tmp_lastbar1 = len(bars) - 1
            bars.append([offsetT + 5 + 3 * i, offsetT + 5 + 3 * i - 2])
            self.selectYPositiveBar(bars)
            tmp_lastbar2 = len(bars) - 1
            bars.append([offsetT + 5 + 3 * i, offsetT + 5 + 3 * i - 4])
            self.selectYNegativeBar(bars)
            bars.append([offsetT + 5 + 3 * i, offsetT + 5 + 3 * i - 5])
            self.selectYPositiveBar(bars)
        self.x_positive_side.append(tmp_lastbar1)
        self.x_positive_side.append(tmp_lastbar2)

        # Linie oben
        for i in range(self.nSGA - 1):
            bars.append([offsetT + 2 + 3 * i, offsetT + 2 + 3 * i + 3])
            self.selectYPositiveBar(bars)
            self.selectYNegativeBar(bars)

    def ausleger_pyramid(self, nodes, bars, offsetT, offsetTG):
        for i in range(1, self.nSA + 1):
            nodes.append([self.ls + i * self.ls, 0, (self.nST-1) * self.ls])  # Left Bottom
            nodes.append([self.ls + i * self.ls, self.ls, (self.nST-1) * self.ls])  # Right Bottom
            nodes.append([(self.ls/2 + i * self.ls), self.ls / 2, (self.nST) * self.ls])  # Top


        # x- und y-Richtung
        for i in range(self.nSA - 1):
            bars.append([offsetTG + i * 3, (offsetTG + 3) + i * 3])
            self.selectYPositiveBar(bars)
            bars.append([offsetTG + i * 3, (offsetTG + 1) + i * 3])
            bars.append([(offsetTG + 1) + i * 3, (offsetTG + 4) + i * 3])
            self.selectYNegativeBar(bars)

        #Hilfslinie unten
        bars.append([offsetT - 2, offsetTG]) #Sonderfall erste Linie
        for i in range(self.nSA -1):
            bars.append([offsetTG + 1 + i*3, offsetTG + 3 + i*3])

        # Bottom nodes with top nodes
        tmp_lastbar1 = 0
        tmp_lastbar2 = 0
        for i in range(self.nSA - 1):
            bars.append([offsetTG + i * 3, (offsetTG + 5) + i * 3])
            self.selectYPositiveBar(bars)
            bars.append([offsetTG + 1 + i * 3, (offsetTG + 5) + i * 3])
            self.selectYNegativeBar(bars)
            bars.append([(offsetTG + 5) + i * 3, offsetTG + 3 + i * 3])
            self.selectYPositiveBar(bars)
            tmp_lastbar1 = len(bars) - 1
            bars.append([(offsetTG + 5) + i * 3, offsetTG + 4 + i * 3])
            self.selectYNegativeBar(bars)
            tmp_lastbar2 = len(bars) - 1
        self.x_negative_side.append(tmp_lastbar1)
        self.x_negative_side.append(tmp_lastbar2)

        # Top Row
        for i in range(self.nSA - 1):
            bars.append([(offsetTG + 2) + i * 3, (offsetTG + 5) + i * 3])
            self.selectYPositiveBar(bars)
            self.selectYNegativeBar(bars)

        # Extra bars
        offsetTO = len(nodes)  # offset after all the nodes

        # Last bar ate the end of the crane
        bars.append([offsetTO - 2, offsetTO - 3])
        self.selectXNegativeBar(bars)

        # Top of the Tower with first Node Ausleger
        bars.append([offsetT - 1, offsetTG + 2])
        self.selectYPositiveBar(bars)
        self.selectYNegativeBar(bars)

        # Tower with the base of the Ausleger
        bars.append([offsetT - 4, offsetTG])
        self.selectYPositiveBar(bars)
        bars.append([offsetT - 2, offsetTG + 1])
        self.selectYNegativeBar(bars)

        # Tower with the first Top Node in Ausleger
        bars.append([offsetT - 4, offsetTG + 2])
        self.selectYPositiveBar(bars)
        bars.append([offsetT - 2, offsetTG + 2])
        self.selectYNegativeBar(bars)

        # First Top Node to base Ausleger2
        bars.append([offsetTG + 2, offsetTG + 1])
        self.selectYNegativeBar(bars)
        bars.append([offsetTG + 2, offsetTG])
        self.selectYPositiveBar(bars)

    def selectYNegativeBar(self, bars):
        """
        Select the last bar from the bar array and add this to another array to select the Y bars
        """
        self.y_negative_side.append(len(bars) - 1)

    def selectYPositiveBar(self, bars):
        """
        Select the last bar from the bar array and add this to another array to select the Y bars
        """
        self.y_positive_side.append(len(bars) - 1)

    def selectXNegativeBar(self, bars):
        """
        Select the last bar from the bar array and add this to another array to select the X negative bars
        """
        self.x_negative_side.append(len(bars) - 1)

    def selectXPositiveBar(self, bars):
        """
        Select the last bar from the bar array and add this to another array to select the X positive bars
        """
        self.x_positive_side.append(len(bars) - 1)

class crane_2_1(Truss):
    """
    Special truss that represents a crane
    """

    def __init__(self, height, length, ls, A, rho, E, p=True):
        """
        :param variant:
        :param height:
        :param length:
        :param ls:
        """
        if ls > height:
            raise Exception("Height of segments cannot be greater than the height of the crane!")
        if ls > length:
            raise Exception("Length of segments cannot be greater than the length of the jib!")

        nodes = []
        bars = []

        self.height = height
        self.length = length
        self.ls = ls

        self.nST = np.ceil(height / ls).astype('int')  # Number of segments of the Tower
        self.nSA = np.ceil(length / ls).astype('int')  # Number of segments of the Ausleger
        self.nSGA = np.ceil((length / 2) / ls).astype('int')  # Number of segments of the Gegenausleger

        self.ls = np.min([self.height / self.nST, self.length / self.nSA])

        # indices of bars on a certain side of the crane
        self.x_negative_side = []  # Done
        self.x_positive_side = []  # Done
        self.y_negative_side = []
        self.y_positive_side = []  # Done

        if p:
            print("Creating crane with zig-zag tower with " + str(self.nST) + " segments of length "+str(self.ls)+" and pyramidal jib with " + str(self.nSA) + " segments.")
        self.tower_ver2(nodes, bars)
        offsetT = cur_offset(nodes) - 8
        self.ausleger_ver2(nodes, bars, offsetT)
        offsetA = cur_offset(nodes) + 1
        self.tip_nodes = [offsetA-3, offsetA-4] # Knoten an der Spitze des Auslegers zur Lastanbringung
        self.counterweight_nodes = [] # Knoten des Gegenauslegers zur Anbringung der Gegenlast
        self.gegenausleger_ver2(nodes, bars, offsetT, offsetA)

        # convert python list to np.array
        self.nodes = np.array(nodes).astype(float)
        self.bars = np.array(bars)

        self.x_negative_side = np.array(self.x_negative_side).astype(int)
        self.x_positive_side = np.array(self.x_positive_side).astype(int)
        self.y_positive_side = np.array(self.y_positive_side).astype(int)
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

    def selectYNegativeBar(self, bars):
        """
        Select the last bar from the bar array and add this to another array to select the Y bars
        """
        self.y_negative_side.append(len(bars) - 1)

    def selectYPositiveBar(self, bars):
        """
        Select the last bar from the bar array and add this to another array to select the Y bars
        """
        self.y_positive_side.append(len(bars) - 1)

    def selectXNegativeBar(self, bars):
        """
        Select the last bar from the bar array and add this to another array to select the X negative bars
        """
        self.x_negative_side.append(len(bars) - 1)

    def selectXPositiveBar(self, bars):
        """
        Select the last bar from the bar array and add this to another array to select the X positive bars
        """
        self.x_positive_side.append(len(bars) - 1)

    def tower_ver2(self, nodes, bars):
        # Nodes des Turms
        for i in range(self.nST-1):
            nodes.append([0, 0, i * self.ls])  # Left Top
            nodes.append([self.ls, 0, i * self.ls])  # Right Top
            nodes.append([0, self.ls, i * self.ls])  # Left Bottom
            nodes.append([self.ls, self.ls, i * self.ls])  # Right Bottom
        #Spitze
        nodes.append([self.ls*0.25, self.ls*0.25, self.ls * self.nST - self.ls])
        nodes.append([self.ls*0.25, self.ls*0.75, self.ls * self.nST - self.ls])
        nodes.append([self.ls*0.75, self.ls*0.25, self.ls * self.nST - self.ls])
        nodes.append([self.ls*0.75, self.ls*0.75, self.ls * self.nST - self.ls])
        nodes.append([self.ls/2, self.ls/2, self.ls * self.nST])

        # Bars des Turms
        
        # z-Richtung
        for i in range(self.nST - 2):
            bars.append([4 * i, 4 * i + 4])  # LT
            self.selectYPositiveBar(bars)
            self.selectXPositiveBar(bars)
            bars.append([4 * i + 1, 4 * i + 5])  # RT
            self.selectYPositiveBar(bars)
            self.selectXNegativeBar(bars)
            bars.append([4 * i + 2, 4 * i + 6])  # LB
            self.selectXPositiveBar(bars)
            self.selectYNegativeBar(bars)
            bars.append([4 * i + 3, 4 * i + 7])  # RB
            self.selectXNegativeBar(bars)
            self.selectYNegativeBar(bars)

        # Kreuzstreben
        weird_offset = self.nST//2+1
        for i in range(self.nST - int(np.ceil(self.nST/2))):

            if i != self.nST - weird_offset:
                bars.append([8 * i, 8 * i + 5])  # LT -> RT+1
                self.selectYPositiveBar(bars)
                self.selectYNegativeBar(bars)
            if i != 0:
                bars.append([8 * i - 3, 8 * i])  # RT -> RT+1
                self.selectYPositiveBar(bars)
                self.selectYNegativeBar(bars)

            if i != self.nST - weird_offset:
                bars.append([8 * i + 3, 8 * i + 6])  # RB -> LB+1
                self.selectYPositiveBar(bars)
                self.selectYNegativeBar(bars)
            if i != 0:
                bars.append([8 * i - 2, 8 * i + 3])  # LB -> RB+1
                self.selectYPositiveBar(bars)
                self.selectYNegativeBar(bars)

            if i != self.nST - weird_offset:
                bars.append([8 * i + 1, 8 * i + 7])  # RT -> RB+1
                self.selectXNegativeBar(bars)
                self.selectXPositiveBar(bars)
            if i != 0:
                bars.append([8 * i - 1, 8 * i + 1])  # RB -> RT+1
                self.selectXNegativeBar(bars)
                self.selectXPositiveBar(bars)

            if i != self.nST - weird_offset:
                bars.append([8 * i + 2, 8 * i + 4])  # LB -> LT+1
                self.selectXNegativeBar(bars)
                self.selectXPositiveBar(bars)
            if i != 0:
                bars.append([8 * i - 4, 8 * i + 2])  # LT -> LB+1
                self.selectXNegativeBar(bars)
                self.selectXPositiveBar(bars)

        offset = len(nodes)
        #aller oberste Spitze
        bars.append([offset-2, offset-6])
        self.selectYNegativeBar(bars)
        bars.append([offset-3, offset-8])
        self.selectYNegativeBar(bars)
        bars.append([offset-4, offset-7])
        self.selectYPositiveBar(bars)
        bars.append([offset-5, offset-9])
        self.selectYPositiveBar(bars)
        bars.append([offset-2, offset-1])
        self.selectYNegativeBar(bars)
        bars.append([offset-3, offset-1])
        self.selectYNegativeBar(bars)
        bars.append([offset-4, offset-1])
        self.selectYPositiveBar(bars)
        bars.append([offset-5, offset-1])
        self.selectYPositiveBar(bars)

        bars.append([offset-6, offset-7]) # x-Richtung
        bars.append([offset-8, offset-9]) # x-Richtung
        bars.append([offset-6, offset-9]) # Diagonale unten

        # spitze verstärkung mitte 
        # bars.append([offset-2, offset-3]) # y-Richtung
        # bars.append([offset-4, offset-5]) # y-Richtung
        # bars.append([offset-2, offset-4]) # y-Richtung
        # bars.append([offset-3, offset-5])

        #kreuzstreben der Spitze
        bars.append([offset-6, offset-4])
        self.selectYNegativeBar(bars)
        bars.append([offset-7, offset-5])
        self.selectYNegativeBar(bars)
        bars.append([offset-8, offset-2])
        self.selectYNegativeBar(bars)
        bars.append([offset-9, offset-3])
        self.selectYNegativeBar(bars)

    def ausleger_ver2(self, nodes, bars, offsetT):
        print("last index tower:", len(nodes)-1)
        for i in range(1, self.nSA + 1):
            nodes.append([self.ls + i * self.ls, 0, (self.nST-2) * self.ls])  # Left Bottom
            nodes.append([self.ls + i * self.ls, self.ls, (self.nST-2) * self.ls])  # Right Bottom
            nodes.append([(self.ls/2 + i * self.ls), self.ls / 2, (self.nST-1) * self.ls])  # Top

        # x- und y-Richtung
        bars.append([offsetT, offsetT+2])
        bars.append([offsetT, offsetT + 8])
        bars.append([offsetT+2, offsetT + 9])
        for i in range(self.nSA - 1):
            bars.append([offsetT+8 + 3*i, offsetT+9 + 3*i])
            self.selectYPositiveBar(bars)
            bars.append([offsetT+8 + 3*i, offsetT + 11 + 3*i])
            self.selectYNegativeBar(bars)
            bars.append([offsetT+9 + 3*i, offsetT + 12 + 3*i])
            self.selectYNegativeBar(bars)
        bars.append([offsetT+8 + (self.nSA - 1)*3, offsetT+9 + (self.nSA - 1)*i])

        # making the pyramids
        bars.append([offsetT, offsetT + 10])
        self.selectYPositiveBar(bars)
        bars.append([offsetT + 2, offsetT + 10])
        self.selectYNegativeBar(bars)
        bars.append([offsetT + 8, offsetT + 10])
        self.selectYPositiveBar(bars)
        bars.append([offsetT + 9, offsetT + 10])
        self.selectYNegativeBar(bars)
        tmp_lastbar1 = 0
        tmp_lastbar2 = 0
        for i in range(self.nSA - 1):
            bars.append([offsetT + 8 + i * 3, (offsetT + 8 + 5) + i * 3])
            self.selectYPositiveBar(bars)
            bars.append([offsetT + 8 + 1 + i * 3, (offsetT + 8 + 5) + i * 3])
            self.selectYNegativeBar(bars)
            bars.append([(offsetT + 8 + 5) + i * 3, offsetT + 8 + 3 + i * 3])
            self.selectYPositiveBar(bars)
            tmp_lastbar1 = len(bars) - 1
            bars.append([(offsetT + 8 + 5) + i * 3, offsetT + 8 + 4 + i * 3])
            self.selectYNegativeBar(bars)
            tmp_lastbar2 = len(bars) - 1
        self.x_negative_side.append(tmp_lastbar1)
        self.x_negative_side.append(tmp_lastbar2)

        # Top Row
        for i in range(self.nSA - 1):
            bars.append([(offsetT + 10) + i * 3, (offsetT + 13) + i * 3])
            self.selectYPositiveBar(bars)
            self.selectYNegativeBar(bars)

        # diagonals
        bars.append([offsetT, offsetT+9])
        for i in range(self.nSA - 1):
            bars.append([offsetT+8 + 3*i, offsetT+12 + 3*i])

        # rope
        #bars.append([offsetT+7, offsetT + (self.nSA//2)*3 + 10])

    def gegenausleger_ver2(self, nodes, bars, offsetT, offsetA):
        nodes.append([(self.ls/2 - self.ls), self.ls / 2, (self.nST-1) * self.ls]) # first top
        for i in range(2, self.nSA//2 + 3):
            nodes.append([self.ls - i * self.ls, 0, (self.nST-2) * self.ls])  # Left Bottom
            self.counterweight_nodes.append(len(nodes) - 1)
            nodes.append([self.ls - i * self.ls, self.ls, (self.nST-2) * self.ls])  # Right Bottom
            self.counterweight_nodes.append(len(nodes) - 1)
            if i!=self.nSA//2 + 2:
                nodes.append([(self.ls/2 - i * self.ls), self.ls / 2, (self.nST-1) * self.ls])  # Top

        # x- und y-Richtung
        bars.append([offsetT-1, offsetA])
        bars.append([offsetT+1, offsetA+1])
        bars.append([offsetT-1, offsetT+1]) # y
        for i in range(self.nSA//2):
            bars.append([offsetA + 3*i, offsetA + 1 + 3*i]) # y
            self.selectYPositiveBar(bars)
            bars.append([offsetA + 3*i, offsetA + 3 + 3*i]) # x vorne
            self.selectYNegativeBar(bars)
            bars.append([offsetA + 3*i + 1, offsetA + 4 + 3*i]) # x hinten
            self.selectYNegativeBar(bars)
        bars.append([offsetA + 3*(self.nSA//2), offsetA + 1 + 3*(self.nSA//2)])

        # making the pyramids
        tmp_lastbar1 = 0
        tmp_lastbar2 = 0
        bars.append([offsetT-1, offsetA-1])
        bars.append([offsetT+1, offsetA-1])
        bars.append([offsetA, offsetA-1])
        bars.append([offsetA+1, offsetA-1])
        for i in range(self.nSA//2):
            bars.append([offsetA + 3*i, offsetA + 2 + 3*i])
            self.selectYPositiveBar(bars)
            bars.append([offsetA + 3*i + 1, offsetA + 2 + 3*i])
            self.selectYNegativeBar(bars)
            bars.append([offsetA + 3*i + 4, offsetA + 2 + 3*i])
            self.selectYNegativeBar(bars)
            bars.append([offsetA + 3*i + 3, offsetA + 2 + 3*i])
            self.selectYNegativeBar(bars)
        self.x_negative_side.append(tmp_lastbar1)
        self.x_negative_side.append(tmp_lastbar2)

        # Top Row
        for i in range(self.nSA//2):
            bars.append([offsetA - 1 + i * 3, offsetA + 2 + i * 3])
            self.selectYPositiveBar(bars)
            self.selectYNegativeBar(bars)

        # diagonals
        bars.append([offsetT+1, offsetA])
        for i in range(self.nSA//2):
            bars.append([offsetA+1 + 3*i, offsetA+3 + 3*i])

        # rope
        #bars.append([offsetT+7, -3])

class crane_2_2(Truss):
    """
    Special truss that represents a crane
    """

    def __init__(self, height, length, ls, A, rho, E, p=True):
        """
        :param variant:
        :param height:
        :param length:
        :param ls:
        """
        if ls > height:
            raise Exception("Height of segments cannot be greater than the height of the crane!")
        if ls > length:
            raise Exception("Length of segments cannot be greater than the length of the jib!")

        nodes = []
        bars = []

        self.height = height
        self.length = length
        self.ls = ls

        self.nST = np.ceil(height / ls).astype('int')  # Number of segments of the Tower
        self.nSA = np.ceil(length / ls).astype('int')  # Number of segments of the Ausleger
        self.nSGA = np.ceil((length / 2) / ls).astype('int')  # Number of segments of the Gegenausleger

        self.ls = np.min([self.height / self.nST, self.length / self.nSA])

        # indices of bars on a certain side of the crane
        self.x_negative_side = []  # Done
        self.x_positive_side = []  # Done
        self.y_negative_side = []
        self.y_positive_side = []  # Done

        if p:
            print("Creating crane with zig-zag tower with " + str(self.nST) + " segments of length "+str(self.ls)+" and pyramidal jib with " + str(self.nSA) + " segments.")
        self.tower_ver2(nodes, bars)
        offsetT = cur_offset(nodes)
        self.ausleger_ver2(nodes, bars, offsetT)
        offsetA = cur_offset(nodes) + 1
        self.tip_nodes = [offsetA-3, offsetA-4, offsetA-6, offsetA-7] # Knoten an der Spitze des Auslegers zur Lastanbringung
        self.counterweight_nodes = [] # Knoten des Gegenauslegers zur Anbringung der Gegenlast
        self.gegenausleger_ver2(nodes, bars, offsetT, offsetA)

        # convert python list to np.array
        self.nodes = np.array(nodes).astype(float)
        self.bars = np.array(bars)

        self.x_negative_side = np.array(self.x_negative_side).astype(int)
        self.x_positive_side = np.array(self.x_positive_side).astype(int)
        self.y_positive_side = np.array(self.y_positive_side).astype(int)
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

    def selectYNegativeBar(self, bars):
        """
        Select the last bar from the bar array and add this to another array to select the Y bars
        """
        self.y_negative_side.append(len(bars) - 1)

    def selectYPositiveBar(self, bars):
        """
        Select the last bar from the bar array and add this to another array to select the Y bars
        """
        self.y_positive_side.append(len(bars) - 1)

    def selectXNegativeBar(self, bars):
        """
        Select the last bar from the bar array and add this to another array to select the X negative bars
        """
        self.x_negative_side.append(len(bars) - 1)

    def selectXPositiveBar(self, bars):
        """
        Select the last bar from the bar array and add this to another array to select the X positive bars
        """
        self.x_positive_side.append(len(bars) - 1)

    def tower_ver2(self, nodes, bars):
        # Nodes des Turms
        for i in range(self.nST):
            nodes.append([0, 0, i * self.ls])  # Left Top
            nodes.append([self.ls, 0, i * self.ls])  # Right Top
            nodes.append([0, self.ls, i * self.ls])  # Left Bottom
            nodes.append([self.ls, self.ls, i * self.ls])  # Right Bottom
        nodes.append([self.ls/2, self.ls/2, self.ls * self.nST]) #Spitze

        # Bars des Turms
        
        # z-Richtung
        for i in range(self.nST - 1):
            bars.append([4 * i, 4 * i + 4])  # LT
            self.selectYPositiveBar(bars)
            self.selectXPositiveBar(bars)
            bars.append([4 * i + 1, 4 * i + 5])  # RT
            self.selectYPositiveBar(bars)
            self.selectXNegativeBar(bars)
            bars.append([4 * i + 2, 4 * i + 6])  # LB
            self.selectXPositiveBar(bars)
            self.selectYNegativeBar(bars)
            bars.append([4 * i + 3, 4 * i + 7])  # RB
            self.selectXNegativeBar(bars)
            self.selectYNegativeBar(bars)

        # Kreuzstreben
        weird_offset = self.nST//2+1
        for i in range(self.nST + 1 - int(np.ceil(self.nST/2))):

            if i != self.nST - weird_offset:
                bars.append([8 * i, 8 * i + 5])  # LT -> RT+1
                self.selectYPositiveBar(bars)
                self.selectYNegativeBar(bars)
            if i != 0:
                bars.append([8 * i - 3, 8 * i])  # RT -> RT+1
                self.selectYPositiveBar(bars)
                self.selectYNegativeBar(bars)

            if i != self.nST - weird_offset:
                bars.append([8 * i + 3, 8 * i + 6])  # RB -> LB+1
                self.selectYPositiveBar(bars)
                self.selectYNegativeBar(bars)
            if i != 0:
                bars.append([8 * i - 2, 8 * i + 3])  # LB -> RB+1
                self.selectYPositiveBar(bars)
                self.selectYNegativeBar(bars)

            if i != self.nST - weird_offset:
                bars.append([8 * i + 1, 8 * i + 7])  # RT -> RB+1
                self.selectXNegativeBar(bars)
                self.selectXPositiveBar(bars)
            if i != 0:
                bars.append([8 * i - 1, 8 * i + 1])  # RB -> RT+1
                self.selectXNegativeBar(bars)
                self.selectXPositiveBar(bars)

            if i != self.nST - weird_offset:
                bars.append([8 * i + 2, 8 * i + 4])  # LB -> LT+1
                self.selectXNegativeBar(bars)
                self.selectXPositiveBar(bars)
            if i != 0:
                bars.append([8 * i - 4, 8 * i + 2])  # LT -> LB+1
                self.selectXNegativeBar(bars)
                self.selectXPositiveBar(bars)

        offset = len(nodes)
        #aller oberste Spitze
        bars.append([offset-2, offset-1])
        self.selectYNegativeBar(bars)
        bars.append([offset-3, offset-1])
        self.selectYNegativeBar(bars)
        bars.append([offset-4, offset-1])
        self.selectYPositiveBar(bars)
        bars.append([offset-5, offset-1])
        self.selectYPositiveBar(bars)

        bars.append([offset-2, offset-3]) # x-Richtung
        bars.append([offset-4, offset-5]) # x-Richtung
        bars.append([offset-2, offset-5]) # Diagonale unten

    def ausleger_ver2(self, nodes, bars, offsetT):

        for i in range(1, self.nSA + 1):
            nodes.append([self.ls + i * self.ls, 0, (self.nST-1) * self.ls])  # Left Bottom
            nodes.append([self.ls + i * self.ls, self.ls, (self.nST-1) * self.ls])  # Right Bottom
            nodes.append([(self.ls/2 + i * self.ls), self.ls / 2, self.nST * self.ls])  # Top

        # x- und y-Richtung
        bars.append([offsetT-2, offsetT-4])
        bars.append([offsetT-4, offsetT + 3])
        bars.append([offsetT-2, offsetT + 4])
        for i in range(self.nSA - 1):
            bars.append([offsetT + 3*i, offsetT+1 + 3*i]) # y
            self.selectYPositiveBar(bars)
            bars.append([offsetT + 3*i, offsetT + 3 + 3*i])
            self.selectYNegativeBar(bars)
            bars.append([offsetT+1 + 3*i, offsetT + 4 + 3*i])
            self.selectYNegativeBar(bars)
        bars.append([offsetT + (self.nSA - 1)*3, offsetT+1 + (self.nSA - 1)*i])

        # making the pyramids
        bars.append([offsetT, offsetT + 2])
        self.selectYPositiveBar(bars)
        bars.append([offsetT + 1, offsetT + 2])
        self.selectYNegativeBar(bars)
        bars.append([offsetT - 2, offsetT + 2])
        self.selectYPositiveBar(bars)
        bars.append([offsetT - 4, offsetT + 2])
        self.selectYNegativeBar(bars)
        tmp_lastbar1 = 0
        tmp_lastbar2 = 0
        for i in range(self.nSA - 1):
            bars.append([offsetT + i * 3, (offsetT + 5) + i * 3])
            self.selectYPositiveBar(bars)
            bars.append([offsetT + 1 + i * 3, (offsetT + 5) + i * 3])
            self.selectYNegativeBar(bars)
            bars.append([(offsetT + 5) + i * 3, offsetT + 3 + i * 3])
            self.selectYPositiveBar(bars)
            tmp_lastbar1 = len(bars) - 1
            bars.append([(offsetT + 5) + i * 3, offsetT + 4 + i * 3])
            self.selectYNegativeBar(bars)
            tmp_lastbar2 = len(bars) - 1
        self.x_negative_side.append(tmp_lastbar1)
        self.x_negative_side.append(tmp_lastbar2)

        # Top Row
        #bars.append([offsetT-1, offsetT + 2])
        for i in range(self.nSA - 1):
            bars.append([(offsetT + 2) + i * 3, (offsetT + 5) + i * 3])
            self.selectYPositiveBar(bars)
            self.selectYNegativeBar(bars)

        # diagonals
        bars.append([offsetT-4, offsetT+1])
        for i in range(self.nSA - 1):
            bars.append([offsetT + 3*i, offsetT+4 + 3*i])


    def gegenausleger_ver2(self, nodes, bars, offsetT, offsetA):
        nodes.append([(self.ls/2 - self.ls), self.ls / 2, self.nST * self.ls]) # first top
        for i in range(2, self.nSA//2 + 3):
            nodes.append([self.ls - i * self.ls, 0, (self.nST-1) * self.ls])  # Left Bottom
            self.counterweight_nodes.append(len(nodes) - 1)
            nodes.append([self.ls - i * self.ls, self.ls, (self.nST-1) * self.ls])  # Right Bottom
            self.counterweight_nodes.append(len(nodes) - 1)
            if i!=self.nSA//2 + 2:
                nodes.append([(self.ls/2 - i * self.ls), self.ls / 2, self.nST * self.ls])  # Top

        # x- und y-Richtung
        bars.append([offsetT-3, offsetA+1]) # erstes x
        bars.append([offsetT-5, offsetA]) # erstes x
        bars.append([offsetT-3, offsetT-5]) # y
        for i in range(self.nSA//2):
            bars.append([offsetA + 3*i, offsetA + 1 + 3*i]) # y
            self.selectYPositiveBar(bars)
            bars.append([offsetA + 3*i, offsetA + 3 + 3*i]) # x vorne
            self.selectYNegativeBar(bars)
            bars.append([offsetA + 3*i + 1, offsetA + 4 + 3*i]) # x hinten
            self.selectYNegativeBar(bars)
        bars.append([offsetA + 3*(self.nSA//2), offsetA + 1 + 3*(self.nSA//2)])

        # making the pyramids
        tmp_lastbar1 = 0
        tmp_lastbar2 = 0
        bars.append([offsetT-5, offsetA-1])
        bars.append([offsetT-3, offsetA-1])
        bars.append([offsetA, offsetA-1])
        bars.append([offsetA+1, offsetA-1])
        for i in range(self.nSA//2):
            bars.append([offsetA + 3*i, offsetA + 2 + 3*i])
            self.selectYPositiveBar(bars)
            bars.append([offsetA + 3*i + 1, offsetA + 2 + 3*i])
            self.selectYNegativeBar(bars)
            bars.append([offsetA + 3*i + 4, offsetA + 2 + 3*i])
            self.selectYNegativeBar(bars)
            bars.append([offsetA + 3*i + 3, offsetA + 2 + 3*i])
            self.selectYNegativeBar(bars)
        self.x_negative_side.append(tmp_lastbar1)
        self.x_negative_side.append(tmp_lastbar2)

        # Top Row
        #bars.append([offsetA - 1, offsetT-1]) # verbindung zur Spitze
        for i in range(self.nSA//2):
            bars.append([offsetA - 1 + i * 3, offsetA + 2 + i * 3])
            self.selectYPositiveBar(bars)
            self.selectYNegativeBar(bars)

        # diagonals
        bars.append([offsetT-3, offsetA])
        for i in range(self.nSA//2):
            bars.append([offsetA+1 + 3*i, offsetA+3 + 3*i])

        # rope
        #bars.append([offsetT+7, -3])

def cur_offset(nodes):
    return len(nodes)
    