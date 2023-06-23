from truss import Truss
import numpy as np


class crane_1(Truss):
    """
    Simple crane with cubical tower and jib made of pyramids
    """

    def __init__(self, height, length, ls, A, rho, E, p=True, max_bar_length=-1):
        """
        height : float
            Maximum crane height. Not exact, because height depends on the length of the segments.
        length : float
            Maximum length of jib. Not exact, because jib length depends on the length of the segments.
        ls : float
            the length of the segments
        A : float (or vector of float)
            crossection area(s) of bars
        rho: float
            desity of material
        E : float
            Young's modulus of material
        p : bool
            print information after creation
        max_bar_length : float
            Maximal allowed length of beams. Throws exception if surpassed.
        optimized_crossections : bool
            will increase crossections of particularly strained beams up to ten times if true
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

        self.nST = np.ceil(height / ls).astype('int')  # number of segments of the tower
        self.nSA = np.ceil(length / ls).astype('int')  # number of segments of the Jib
        self.nSGA = np.ceil((length / 2) / ls).astype('int')  # number of segments of the counter Jib
        self.ls = np.min([self.height / self.nST, self.length / self.nSA])

        # indices of bars on a certain side of the crane
        self.x_positive_side = []  
        self.x_negative_side = []  
        self.y_positive_side = []
        self.y_negative_side = [] 

        if p:
            print("Creating crane with cuboid tower with " + str(self.nST) + " segments of length "+str(self.ls)+" and pyramidal jib with " + str(self.nSA) + " segments.")
        self.make_tower(nodes, bars)
        offsetT = cur_offset(nodes)
        self.counterweight_nodes = [] # nodes of counter jib for counter weight 
        self.make_counterjib(nodes, bars, offsetT)
        offsetTG = cur_offset(nodes)
        self.make_jib(nodes, bars, offsetT, offsetTG)
        self.tip_nodes = [-2, -3, -5, -6] # nodes at front of jib where weight is applied
        super().__init__(nodes, bars, A, rho, E)

        # indices of bars on a certain side of the crane
        self.x_positive_side = np.array(self.x_positive_side).astype(int)
        self.x_negative_side = np.array(self.x_negative_side).astype(int)
        self.y_positive_side = np.array(self.y_positive_side).astype(int)
        self.y_negative_side = np.array(self.y_negative_side).astype(int)

        # supports
        self.supports = np.ones_like(self.nodes).astype(int)
        self.Ur = np.array([]).astype('int')
        for i in range(4):
            self.addSupport(i, 0, 0, 0)

        self._computeLengths()
        self._computeOrientations()
        self._computeMass()

        if max_bar_length != -1 and max_bar_length > np.max(self.lengths):
            raise Exception("Maximum bar length exceeded!\n Maximum length: "+str(np.max(self.lengths)))

    def __str__(self):
        return self.to_string()

    def to_string(self):
        return "Kran Modell 1"
    
    def make_tower(self, nodes, bars):
        # nodes
        for i in range(self.nST):
            nodes.append([0, 0, i * self.ls])  
            nodes.append([self.ls, 0, i * self.ls])  
            nodes.append([0, self.ls, i * self.ls])  
            nodes.append([self.ls, self.ls, i * self.ls])  
        nodes.append([self.ls / 2, self.ls / 2, self.ls * self.nST])

        # bars 
        # x- and y-direction (LT for Left Top usw.)
        for i in range(self.nST):
            bars.append([4 * i, 4 * i + 1])  
            selectYNegativeBar(self, bars)
            bars.append([4 * i + 2, 4 * i + 3])  
            selectYPositiveBar(self, bars)
            bars.append([4 * i, 4 * i + 2]) 
            selectXNegativeBar(self, bars)
            bars.append([4 * i + 1, 4 * i + 3])  
            selectXPositiveBar(self, bars)
            
        # z-direction
        for i in range(self.nST - 1):
            bars.append([4 * i, 4 * i + 4])  # LT
            selectYNegativeBar(self, bars)
            selectXNegativeBar(self, bars)
            bars.append([4 * i + 1, 4 * i + 5])  # RT
            selectYNegativeBar(self, bars)
            selectXPositiveBar(self, bars)
            bars.append([4 * i + 2, 4 * i + 6])  # LB
            selectXNegativeBar(self, bars)
            selectYPositiveBar(self, bars)
            bars.append([4 * i + 3, 4 * i + 7])  # RB
            selectXPositiveBar(self, bars)
            selectYPositiveBar(self, bars)

        # crossbars 
        for i in range(self.nST - 1):
            bars.append([4 * i, 4 * i + 5])  
            selectYNegativeBar(self, bars)
            selectYPositiveBar(self, bars)
            
            bars.append([4 * i + 3, 4 * i + 6])  
            selectYNegativeBar(self, bars)
            selectYPositiveBar(self, bars)

            bars.append([4 * i + 1, 4 * i + 7]) 
            selectXPositiveBar(self, bars)
            selectXNegativeBar(self, bars)

            bars.append([4 * i + 2, 4 * i + 4])  
            selectXPositiveBar(self, bars)
            selectXNegativeBar(self, bars)

        # top
        offsetTO = len(nodes)
        bars.append([offsetTO - 1, offsetTO - 2])
        selectYPositiveBar(self, bars)
        bars.append([offsetTO - 1, offsetTO - 3])
        selectYPositiveBar(self, bars)
        bars.append([offsetTO - 1, offsetTO - 4])
        selectYNegativeBar(self, bars)
        bars.append([offsetTO - 1, offsetTO - 5])
        selectYNegativeBar(self, bars)
        
        # # diagonal lines at the top (optional)
        # bars.append([offsetTO - 3, offsetTO - 4])
        # bars.append([offsetTO - 2, offsetTO - 5])

    def make_counterjib(self, nodes, bars, offsetT):
        # nodes 
        for i in range(1, self.nSGA + 1):  
            nodes.append([-(self.ls + i * self.ls) + self.ls, 0, (self.nST - 1) * self.ls])  
            nodes.append([-(self.ls + i * self.ls) + self.ls, self.ls, (self.nST - 1) * self.ls])  
            self.counterweight_nodes.append(len(nodes) - 2)
            self.counterweight_nodes.append(len(nodes) - 1)
            nodes.append([-(0.5 * self.ls + i * self.ls) + self.ls, self.ls / 2,(self.nST) * self.ls])  

        # bars 
        # first pyramid
        bars.append([offsetT, offsetT - 5])
        selectYNegativeBar(self, bars)
        bars.append([offsetT + 1, offsetT - 3])
        selectYPositiveBar(self, bars)
        bars.append([offsetT + 2, offsetT + 1])
        selectYPositiveBar(self, bars)
        bars.append([offsetT + 2, offsetT])
        selectYNegativeBar(self, bars)
        bars.append([offsetT + 2, offsetT - 5])
        selectYNegativeBar(self, bars)
        bars.append([offsetT + 2, offsetT - 3])
        selectYPositiveBar(self, bars)

        # connection to peak of tower
        bars.append([offsetT + 2, offsetT - 1])
        selectYNegativeBar(self, bars)
        selectYPositiveBar(self, bars)

        # x- and y-direction 
        for i in range(self.nSGA - 1):
            bars.append([offsetT + 3 * i, offsetT + 3 + 3 * i])
            selectYNegativeBar(self, bars)
            bars.append([offsetT + 1 + 3 * i, offsetT + 1 + 3 + 3 * i])
            selectYPositiveBar(self, bars)
            bars.append([offsetT + 3 * i, offsetT + 1 + 3 * i])
        offsetGT = len(nodes)
        bars.append([offsetGT - 2, offsetGT - 3])  # last bar
        selectXNegativeBar(self, bars)

        # # diagonal lines at the bottom (optional) 
        # bars.append([offsetT + 1, offsetT -5])  #first line
        # for i in range(self.nSGA - 1):
        #     bars.append([offsetT + 3 * i, offsetT + 4 + 3 * i])

        # making the pyramids
        tmp_lastbar1 = 0
        tmp_lastbar2 = 0
        for i in range(self.nSGA - 1):  # starting from second pyramid
            bars.append([offsetT + 5 + 3 * i, offsetT + 5 + 3 * i - 1])
            selectYPositiveBar(self, bars)
            tmp_lastbar1 = len(bars) - 1
            bars.append([offsetT + 5 + 3 * i, offsetT + 5 + 3 * i - 2])
            selectYNegativeBar(self, bars)
            tmp_lastbar2 = len(bars) - 1
            bars.append([offsetT + 5 + 3 * i, offsetT + 5 + 3 * i - 4])
            selectYPositiveBar(self, bars)
            bars.append([offsetT + 5 + 3 * i, offsetT + 5 + 3 * i - 5])
            selectYNegativeBar(self, bars)
        self.x_negative_side.append(tmp_lastbar1)
        self.x_negative_side.append(tmp_lastbar2)

        # line on top
        for i in range(self.nSGA - 1):
            bars.append([offsetT + 2 + 3 * i, offsetT + 2 + 3 * i + 3])
            selectYNegativeBar(self, bars)
            selectYPositiveBar(self, bars)


    def make_jib(self, nodes, bars, offsetT, offsetTG):
        # nodes
        for i in range(1, self.nSA + 1):
            nodes.append([self.ls + i * self.ls, 0, (self.nST-1) * self.ls])  # Left Bottom
            nodes.append([self.ls + i * self.ls, self.ls, (self.nST-1) * self.ls])  # Right Bottom
            nodes.append([(self.ls/2 + i * self.ls), self.ls / 2, (self.nST) * self.ls])  # Top
        offsetTO = len(nodes)  

        # bars
        # x- und y-dircetion
        for i in range(self.nSA - 1):
            bars.append([offsetTG + i * 3, (offsetTG + 3) + i * 3])
            selectYNegativeBar(self, bars)
            bars.append([offsetTG + i * 3, (offsetTG + 1) + i * 3])
            bars.append([(offsetTG + 1) + i * 3, (offsetTG + 4) + i * 3])
            selectYPositiveBar(self, bars)

        # # diagonal lines at the bottom (optional)
        # bars.append([offsetT - 2, offsetTG]) #first line
        # for i in range(self.nSA -1):
        #     bars.append([offsetTG + 1 + i*3, offsetTG + 3 + i*3])

        # bottom nodes with top nodes
        tmp_lastbar1 = 0
        tmp_lastbar2 = 0
        for i in range(self.nSA - 1):
            bars.append([offsetTG + i * 3, (offsetTG + 5) + i * 3])
            selectYNegativeBar(self, bars)
            bars.append([offsetTG + 1 + i * 3, (offsetTG + 5) + i * 3])
            selectYPositiveBar(self, bars)
            bars.append([(offsetTG + 5) + i * 3, offsetTG + 3 + i * 3])
            selectYNegativeBar(self, bars)
            tmp_lastbar1 = len(bars) - 1
            bars.append([(offsetTG + 5) + i * 3, offsetTG + 4 + i * 3])
            selectYPositiveBar(self, bars)
            tmp_lastbar2 = len(bars) - 1
        self.x_positive_side.append(tmp_lastbar1)
        self.x_positive_side.append(tmp_lastbar2)

        # top row
        for i in range(self.nSA - 1):
            bars.append([(offsetTG + 2) + i * 3, (offsetTG + 5) + i * 3])
            selectYNegativeBar(self, bars)
            selectYPositiveBar(self, bars)


        bars.append([offsetTO - 2, offsetTO - 3]) # last bar at the end of the crane
        selectXPositiveBar(self, bars)

        bars.append([offsetT - 1, offsetTG + 2]) # top of the Tower with first jib node
        selectYNegativeBar(self, bars)
        selectYPositiveBar(self, bars)

        # tower with the base of the jib
        bars.append([offsetT - 4, offsetTG])
        selectYNegativeBar(self, bars)
        bars.append([offsetT - 2, offsetTG + 1])
        selectYPositiveBar(self, bars)

        # tower with the first Top Node in jib
        bars.append([offsetT - 4, offsetTG + 2])
        selectYNegativeBar(self, bars)
        bars.append([offsetT - 2, offsetTG + 2])
        selectYPositiveBar(self, bars)

        # first Top Node to base jib
        bars.append([offsetTG + 2, offsetTG + 1])
        selectYPositiveBar(self, bars)
        bars.append([offsetTG + 2, offsetTG])
        selectYNegativeBar(self, bars)

class crane_2_1(Truss):
    """
    Crane with zig-zag tower, peak, jib made of pyramids and ropes
    """

    def __init__(self, height, length, ls, A, rho, E, p=True):
        """
        height : float
            Maximum crane height. Not exact, because height depends on the length of the segments.
        length : float
            Maximum length of jib. Not exact, because jib length depends on the length of the segments.
        ls : float
            the length of the segments
        A : float
            crossection area of bars
        rho: float
            desity of material
        E : float
            Young's modulus of material
        p : bool
            print information after creation
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

        self.nST = np.ceil(height / ls).astype('int')  # number of segments of the tower
        self.nSA = np.ceil(length / ls).astype('int')  # number of segments of the jib
        self.nSGA = np.ceil((length / 2) / ls).astype('int')  # number of segments of the counter Jib

        self.ls = np.min([self.height / self.nST, self.length / self.nSA])

        # indices of bars on a certain side of the crane
        self.x_positive_side = []  
        self.x_negative_side = []  
        self.y_positive_side = []
        self.y_negative_side = []

        if p:
            print("Creating crane with zig-zag tower with " + str(self.nST) + " segments of length "+str(self.ls)+" and pyramidal jib with " + str(self.nSA) + " segments.")
        self.make_tower(nodes, bars)
        offsetT = cur_offset(nodes) - 8
        self.make_jib(nodes, bars, offsetT)
        offsetA = cur_offset(nodes) + 1
        self.tip_nodes = [offsetA-3, offsetA-4, offsetA-6, offsetA-7] # nodes at front of jib for weights
        self.counterweight_nodes = [] # nodes of counter jib for counter weight
        self.make_counterjib(nodes, bars, offsetT, offsetA)

        super().__init__(nodes, bars, A, rho, E)

        self.x_positive_side = np.array(self.x_positive_side).astype(int)
        self.x_negative_side = np.array(self.x_negative_side).astype(int)
        self.y_positive_side = np.array(self.y_positive_side).astype(int)
        self.y_negative_side = np.array(self.y_negative_side).astype(int)

        # supports
        self.supports = np.ones_like(self.nodes).astype(int)
        self.Ur = np.array([]).astype('int')
        for i in range(4):
            self.addSupport(i, 0, 0, 0)

        # external Forces
        self.F = np.zeros_like(self.nodes)

        self._computeLengths()
        self._computeOrientations()
        self._computeMass()

    def __str__(self):
        return self.to_string()

    def to_string(self):
        return "Kran Modell 2"

    def make_tower(self, nodes, bars):
        # nodes 
        for i in range(self.nST-1):
            nodes.append([0, 0, i * self.ls])  # Left Top
            nodes.append([self.ls, 0, i * self.ls])  # Right Top
            nodes.append([0, self.ls, i * self.ls])  # Left Bottom
            nodes.append([self.ls, self.ls, i * self.ls])  # Right Bottom
        # top
        nodes.append([self.ls*0.25, self.ls*0.25, self.ls * self.nST - self.ls])
        nodes.append([self.ls*0.25, self.ls*0.75, self.ls * self.nST - self.ls])
        nodes.append([self.ls*0.75, self.ls*0.25, self.ls * self.nST - self.ls])
        nodes.append([self.ls*0.75, self.ls*0.75, self.ls * self.nST - self.ls])
        nodes.append([self.ls/2, self.ls/2, self.ls * self.nST])

        # bars 
        # z-direction
        for i in range(self.nST - 2):
            bars.append([4 * i, 4 * i + 4])  # LT
            selectYNegativeBar(self, bars)
            selectXNegativeBar(self, bars)
            bars.append([4 * i + 1, 4 * i + 5])  # RT
            selectYNegativeBar(self, bars)
            selectXPositiveBar(self, bars)
            bars.append([4 * i + 2, 4 * i + 6])  # LB
            selectXNegativeBar(self, bars)
            selectYPositiveBar(self, bars)
            bars.append([4 * i + 3, 4 * i + 7])  # RB
            selectXPositiveBar(self, bars)
            selectYPositiveBar(self, bars)

        # crossbars
        weird_offset = self.nST//2+1
        for i in range(self.nST - int(np.ceil(self.nST/2))):

            if i != self.nST - weird_offset:
                bars.append([8 * i, 8 * i + 5])  
                selectYNegativeBar(self, bars)
                selectYPositiveBar(self, bars)
            if i != 0:
                bars.append([8 * i - 3, 8 * i])  
                selectYNegativeBar(self, bars)
                selectYPositiveBar(self, bars)

            if i != self.nST - weird_offset:
                bars.append([8 * i + 3, 8 * i + 6])  
                selectYNegativeBar(self, bars)
                selectYPositiveBar(self, bars)
            if i != 0:
                bars.append([8 * i - 2, 8 * i + 3])  
                selectYNegativeBar(self, bars)
                selectYPositiveBar(self, bars)

            if i != self.nST - weird_offset:
                bars.append([8 * i + 1, 8 * i + 7])  
                selectXPositiveBar(self, bars)
                selectXNegativeBar(self, bars)
            if i != 0:
                bars.append([8 * i - 1, 8 * i + 1])  
                selectXPositiveBar(self, bars)
                selectXNegativeBar(self, bars)

            if i != self.nST - weird_offset:
                bars.append([8 * i + 2, 8 * i + 4])  
                selectXPositiveBar(self, bars)
                selectXNegativeBar(self, bars)
            if i != 0:
                bars.append([8 * i - 4, 8 * i + 2])  
                selectXPositiveBar(self, bars)
                selectXNegativeBar(self, bars)

        # peak
        offset = len(nodes)
        bars.append([offset-2, offset-6])
        selectYPositiveBar(self, bars)
        bars.append([offset-3, offset-8])
        selectYNegativeBar(self, bars)
        bars.append([offset-4, offset-7])
        selectYPositiveBar(self, bars)
        bars.append([offset-5, offset-9])
        selectYNegativeBar(self, bars)
        bars.append([offset-2, offset-1])
        selectXPositiveBar(self, bars)
        selectYPositiveBar(self, bars)
        bars.append([offset-3, offset-1])
        selectXPositiveBar(self, bars)
        selectYNegativeBar(self, bars)
        bars.append([offset-4, offset-1])
        selectXNegativeBar(self, bars)
        selectYPositiveBar(self, bars)
        bars.append([offset-5, offset-1])
        selectXNegativeBar(self, bars)
        selectYNegativeBar(self, bars)

        bars.append([offset-6, offset-7]) # x-direction
        bars.append([offset-8, offset-9]) # x-direction
        bars.append([offset-6, offset-9]) # diagonal bottom

        # top support in the middle 
        # bars.append([offset-2, offset-3]) # y-direction
        # bars.append([offset-4, offset-5]) # y-direction
        # bars.append([offset-2, offset-4]) # y-direction
        # bars.append([offset-3, offset-5])

        # crossbars of top
        bars.append([offset-6, offset-4])
        selectYPositiveBar(self, bars)
        bars.append([offset-7, offset-5])
        bars.append([offset-8, offset-2])
        bars.append([offset-9, offset-3])
        selectYNegativeBar(self, bars)

    def make_jib(self, nodes, bars, offsetT):
        # nodes
        for i in range(1, self.nSA + 1):
            nodes.append([self.ls + i * self.ls, 0, (self.nST-2) * self.ls])  # Left Bottom
            nodes.append([self.ls + i * self.ls, self.ls, (self.nST-2) * self.ls])  # Right Bottom
            nodes.append([(self.ls/2 + i * self.ls), self.ls / 2, (self.nST-1) * self.ls])  # Top

        # bars
        # x- and y-dircetion
        bars.append([offsetT, offsetT+2])
        bars.append([offsetT, offsetT + 8])
        selectYNegativeBar(self, bars)
        bars.append([offsetT+2, offsetT + 9])
        selectYPositiveBar(self, bars)
        for i in range(self.nSA - 1):
            bars.append([offsetT+8 + 3*i, offsetT+9 + 3*i]) # y
            bars.append([offsetT+8 + 3*i, offsetT + 11 + 3*i]) # x front
            selectYNegativeBar(self, bars)
            bars.append([offsetT+9 + 3*i, offsetT + 12 + 3*i]) # x back
            selectYPositiveBar(self, bars)
        bars.append([offsetT+8 + (self.nSA - 1)*3, offsetT+9 + (self.nSA - 1)*3])

        # making the pyramids
        bars.append([offsetT, offsetT + 10])
        selectYNegativeBar(self, bars)
        bars.append([offsetT + 2, offsetT + 10])
        selectYPositiveBar(self, bars)
        bars.append([offsetT + 8, offsetT + 10])
        selectYNegativeBar(self, bars)
        bars.append([offsetT + 9, offsetT + 10])
        selectYPositiveBar(self, bars)
        tmp_lastbar1 = 0
        tmp_lastbar2 = 0
        for i in range(self.nSA - 1):
            bars.append([offsetT + 8 + i * 3, (offsetT + 8 + 5) + i * 3])
            selectYNegativeBar(self, bars)
            bars.append([offsetT + 8 + 1 + i * 3, (offsetT + 8 + 5) + i * 3])
            selectYPositiveBar(self, bars)
            bars.append([(offsetT + 8 + 5) + i * 3, offsetT + 8 + 3 + i * 3])
            selectYNegativeBar(self, bars)
            tmp_lastbar1 = len(bars) - 1
            bars.append([(offsetT + 8 + 5) + i * 3, offsetT + 8 + 4 + i * 3])
            selectYPositiveBar(self, bars)
            tmp_lastbar2 = len(bars) - 1
        self.x_positive_side.append(tmp_lastbar1)
        self.x_positive_side.append(tmp_lastbar2)

        # top row
        for i in range(self.nSA - 1):
            bars.append([(offsetT + 10) + i * 3, (offsetT + 13) + i * 3])
            selectYNegativeBar(self, bars)
            selectYPositiveBar(self, bars)

        # diagonals
        bars.append([offsetT, offsetT+9])
        for i in range(self.nSA - 1):
            bars.append([offsetT+8 + 3*i, offsetT+12 + 3*i])

        # rope
        bars.append([offsetT+7, offsetT + (self.nSA//2)*3 + 10])

    def make_counterjib(self, nodes, bars, offsetT, offsetA):
        # nodes
        nodes.append([(self.ls/2 - self.ls), self.ls / 2, (self.nST-1) * self.ls]) # first top
        for i in range(2, self.nSA//2 + 3):
            nodes.append([self.ls - i * self.ls, 0, (self.nST-2) * self.ls])  # Left Bottom
            self.counterweight_nodes.append(len(nodes) - 1)
            nodes.append([self.ls - i * self.ls, self.ls, (self.nST-2) * self.ls])  # Right Bottom
            self.counterweight_nodes.append(len(nodes) - 1)
            if i!=self.nSA//2 + 2:
                nodes.append([(self.ls/2 - i * self.ls), self.ls / 2, (self.nST-1) * self.ls])  # Top

        #bars
        # x- and y-dircetion
        bars.append([offsetT-1, offsetA]) # x
        selectYNegativeBar(self, bars)
        bars.append([offsetT+1, offsetA+1]) # x
        selectYPositiveBar(self, bars)
        bars.append([offsetT-1, offsetT+1]) # y
        for i in range(self.nSA//2):
            bars.append([offsetA + 3*i, offsetA + 1 + 3*i]) # y
            bars.append([offsetA + 3*i, offsetA + 3 + 3*i]) # x front
            selectYNegativeBar(self, bars)
            bars.append([offsetA + 3*i + 1, offsetA + 4 + 3*i]) # x back
            selectYPositiveBar(self, bars)
        bars.append([offsetA + 3*(self.nSA//2), offsetA + 1 + 3*(self.nSA//2)])

        # making the pyramids
        bars.append([offsetT-1, offsetA-1])
        selectYNegativeBar(self, bars)
        bars.append([offsetT+1, offsetA-1])
        selectYPositiveBar(self, bars)
        bars.append([offsetA, offsetA-1])
        selectYNegativeBar(self, bars)
        bars.append([offsetA+1, offsetA-1])
        selectYPositiveBar(self, bars)
        for i in range(self.nSA//2):
            bars.append([offsetA + 3*i, offsetA + 2 + 3*i])
            selectYNegativeBar(self, bars)
            bars.append([offsetA + 3*i + 1, offsetA + 2 + 3*i])
            selectYPositiveBar(self, bars)
            bars.append([offsetA + 3*i + 4, offsetA + 2 + 3*i])
            selectYPositiveBar(self, bars)
            bars.append([offsetA + 3*i + 3, offsetA + 2 + 3*i])
            selectYNegativeBar(self, bars)
        self.x_negative_side.append(len(bars)-2)
        self.x_negative_side.append(len(bars)-1)

        # top row
        for i in range(self.nSA//2):
            bars.append([offsetA - 1 + i * 3, offsetA + 2 + i * 3])
            selectYNegativeBar(self, bars)
            selectYPositiveBar(self, bars)

        # diagonals
        bars.append([offsetT+1, offsetA])
        for i in range(self.nSA//2):
            bars.append([offsetA+1 + 3*i, offsetA+3 + 3*i])

        # rope
        bars.append([offsetT+7, -3])

class crane_2_2(Truss):
    """
    crane with almost cubicals in the tower,  jib made from pyramids and rope

    """

    def __init__(self, height, length, ls, A, rho, E, p=True, max_bar_length=-1):
        """
         height : float
            Maximum crane height. Not exact, because height depends on the length of the segments.
        length : float
            Maximum length of jib. Not exact, because jib length depends on the length of the segments.
        ls : float
            the length of the segments
        A : float
            crossection area
        rho: float
            desity of material
        E : float
            Young's modulus of material
        p : bool
            print information after creation
        max_bar_length : float
            Maximal allowed length of beams. Throws exception if surpassed.
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

        self.nST = np.ceil(height / ls).astype('int')  # number of segments of the tower
        self.nSA = np.ceil(length / ls).astype('int')  # number of segments of the jib
        self.nSGA = np.ceil((length / 2) / ls).astype('int')  # number of segments of the counter jib

        self.ls = np.min([self.height / self.nST, self.length / self.nSA])

        # indices of bars on a certain side of the crane
        self.x_positive_side = []  
        self.x_negative_side = []  
        self.y_positive_side = []
        self.y_negative_side = []  

        if p:
            print("Creating crane with zig-zag tower with " + str(self.nST) + " segments of length "+str(self.ls)+" and pyramidal jib with " + str(self.nSA) + " segments.")
        self.make_tower(nodes, bars)
        offsetT = cur_offset(nodes)
        self.make_jib(nodes, bars, offsetT)
        offsetA = cur_offset(nodes) + 1
        self.tip_nodes = [offsetA-3, offsetA-4, offsetA-6, offsetA-7] # nodes at front of jib for weights
        self.counterweight_nodes = [] # nodes of counter jib for counter weight
        self.make_counterjib(nodes, bars, offsetT, offsetA)

        super().__init__(nodes, bars, A, rho, E)

        self.x_positive_side = np.array(self.x_positive_side).astype(int)
        self.x_negative_side = np.array(self.x_negative_side).astype(int)
        self.y_positive_side = np.array(self.y_positive_side).astype(int)
        self.y_negative_side = np.array(self.y_negative_side).astype(int)

        # supports
        self.supports = np.ones_like(self.nodes).astype(int)
        self.Ur = np.array([]).astype('int')
        for i in range(4):
            self.addSupport(i, 0, 0, 0)

        # external forces
        self.F = np.zeros_like(self.nodes)

        # Connection vector of each bar, as ending node - starting node
        self.d = self.nodes[self.bars[:, 1], :] - self.nodes[self.bars[:, 0], :]

        self._computeLengths(self.d)
        self._computeOrientations(self.d)
        self._computeMass()

        if max_bar_length != -1 and max_bar_length > np.max(self.lengths):
            raise Exception("Maximum bar length exceeded!\n Maximum length: "+str(np.max(self.lengths)))

    def __str__(self):
        return self.to_string()

    def to_string(self):
        return "Kran Modell 3"


    def make_tower(self, nodes, bars):
        # nodes 
        for i in range(self.nST):
            nodes.append([0, 0, i * self.ls])  # Left Top
            nodes.append([self.ls, 0, i * self.ls])  # Right Top
            nodes.append([0, self.ls, i * self.ls])  # Left Bottom
            nodes.append([self.ls, self.ls, i * self.ls])  # Right Bottom
        nodes.append([self.ls/2, self.ls/2, self.ls * self.nST]) #Spitze

        # bars des Turms
        # z-dircetion
        for i in range(self.nST - 1):
            bars.append([4 * i, 4 * i + 4])  # LT
            selectYNegativeBar(self, bars)
            selectXNegativeBar(self, bars)
            bars.append([4 * i + 1, 4 * i + 5])  # RT
            selectYNegativeBar(self, bars)
            selectXPositiveBar(self, bars)
            bars.append([4 * i + 2, 4 * i + 6])  # LB
            selectXNegativeBar(self, bars)
            selectYPositiveBar(self, bars)
            bars.append([4 * i + 3, 4 * i + 7])  # RB
            selectXPositiveBar(self, bars)
            selectYPositiveBar(self, bars)

        # crossbars
        if self.nST%2 == 0:
            weird_offset = self.nST//2 - 1
            nseg = self.nST
        else:
            nseg = self.nST+1
            weird_offset = self.nST//2
        
        for i in range(nseg - int(np.ceil(nseg/2))):

            if i != nseg - weird_offset:
                bars.append([8 * i, 8 * i + 5])  # LT -> RT+1
                selectYNegativeBar(self, bars)
                selectYPositiveBar(self, bars)
            if i != 0:
                bars.append([8 * i - 3, 8 * i])  # RT -> RT+1
                selectYNegativeBar(self, bars)
                selectYPositiveBar(self, bars)

            if i != nseg - weird_offset:
                bars.append([8 * i + 3, 8 * i + 6])  # RB -> LB+1
                selectYNegativeBar(self, bars)
                selectYPositiveBar(self, bars)
            if i != 0:
                bars.append([8 * i - 2, 8 * i + 3])  # LB -> RB+1
                selectYNegativeBar(self, bars)
                selectYPositiveBar(self, bars)

            if i != nseg - weird_offset:
                bars.append([8 * i + 1, 8 * i + 7])  # RT -> RB+1
                selectXPositiveBar(self, bars)
                selectXNegativeBar(self, bars)
            if i != 0:
                bars.append([8 * i - 1, 8 * i + 1])  # RB -> RT+1
                selectXPositiveBar(self, bars)
                selectXNegativeBar(self, bars)

            if i != nseg - weird_offset:
                bars.append([8 * i + 2, 8 * i + 4])  # LB -> LT+1
                selectXPositiveBar(self, bars)
                selectXNegativeBar(self, bars)
            if i != 0:
                bars.append([8 * i - 4, 8 * i + 2])  # LT -> LB+1
                selectXPositiveBar(self, bars)
                selectXNegativeBar(self, bars)

        offset = len(nodes)
        
        # top
        bars.append([offset-2, offset-1])
        selectYPositiveBar(self, bars)
        bars.append([offset-3, offset-1])
        selectYPositiveBar(self, bars)
        bars.append([offset-4, offset-1])
        selectYNegativeBar(self, bars)
        bars.append([offset-5, offset-1])
        selectYNegativeBar(self, bars)

        bars.append([offset-2, offset-3]) # x-direction
        selectYPositiveBar(self, bars)
        bars.append([offset-4, offset-5]) # x-direction
        selectYNegativeBar(self, bars)
        bars.append([offset-2, offset-5]) # digonal bottom 

    def make_jib(self, nodes, bars, offsetT):

        for i in range(1, self.nSA + 1):
            nodes.append([self.ls + i * self.ls, 0, (self.nST-1) * self.ls])  # Left Bottom
            nodes.append([self.ls + i * self.ls, self.ls, (self.nST-1) * self.ls])  # Right Bottom
            nodes.append([(self.ls/2 + i * self.ls), self.ls / 2, self.nST * self.ls])  # Top

        # x- and y-dircetion
        bars.append([offsetT-2, offsetT-4])
        bars.append([offsetT-4, offsetT + 3])
        selectYNegativeBar(self, bars)
        bars.append([offsetT-2, offsetT + 4])
        selectYPositiveBar(self, bars)
        for i in range(self.nSA - 1):
            bars.append([offsetT + 3*i, offsetT+1 + 3*i]) # y
            bars.append([offsetT + 3*i, offsetT + 3 + 3*i]) # front
            selectYNegativeBar(self, bars)
            bars.append([offsetT+1 + 3*i, offsetT + 4 + 3*i]) # back
            selectYPositiveBar(self, bars)
        bars.append([offsetT + (self.nSA - 1)*3, offsetT+1 + (self.nSA - 1)*3])

        # making the pyramids
        bars.append([offsetT, offsetT + 2])
        selectYNegativeBar(self, bars)
        bars.append([offsetT + 1, offsetT + 2])
        selectYPositiveBar(self, bars)
        bars.append([offsetT - 2, offsetT + 2])
        selectYPositiveBar(self, bars)
        bars.append([offsetT - 4, offsetT + 2])
        selectYNegativeBar(self, bars)
        for i in range(self.nSA - 1):
            bars.append([offsetT + i * 3, (offsetT + 5) + i * 3])
            selectYNegativeBar(self, bars)
            bars.append([offsetT + 1 + i * 3, (offsetT + 5) + i * 3])
            selectYPositiveBar(self, bars)
            bars.append([(offsetT + 5) + i * 3, offsetT + 3 + i * 3])
            selectYNegativeBar(self, bars)
            bars.append([(offsetT + 5) + i * 3, offsetT + 4 + i * 3])
            selectYPositiveBar(self, bars)
        self.x_positive_side.append(len(bars)-2)
        self.x_positive_side.append(len(bars)-1)
    

        # top row
        bars.append([offsetT-1, offsetT + 2])
        selectYPositiveBar(self, bars)
        selectYNegativeBar(self, bars)
        for i in range(self.nSA - 1):
            bars.append([(offsetT + 2) + i * 3, (offsetT + 5) + i * 3])
            selectYNegativeBar(self, bars)
            selectYPositiveBar(self, bars)

        # diagonals
        bars.append([offsetT-4, offsetT+1])
        for i in range(self.nSA - 1):
            bars.append([offsetT + 3*i, offsetT+4 + 3*i])


    def make_counterjib(self, nodes, bars, offsetT, offsetA):
        nodes.append([(self.ls/2 - self.ls), self.ls / 2, self.nST * self.ls]) # First top
        for i in range(2, self.nSA//2 + 3):
            nodes.append([self.ls - i * self.ls, 0, (self.nST-1) * self.ls])  # Left Bottom
            self.counterweight_nodes.append(len(nodes) - 1)
            nodes.append([self.ls - i * self.ls, self.ls, (self.nST-1) * self.ls])  # Right Bottom
            self.counterweight_nodes.append(len(nodes) - 1)
            if i!=self.nSA//2 + 2:
                nodes.append([(self.ls/2 - i * self.ls), self.ls / 2, self.nST * self.ls])  # Top

        # x- and y-dircetion
        bars.append([offsetT-3, offsetA+1]) # First x back
        selectYPositiveBar(self, bars)
        bars.append([offsetT-5, offsetA]) # First x front
        selectYNegativeBar(self, bars)
        bars.append([offsetT-3, offsetT-5]) # y
        for i in range(self.nSA//2):
            bars.append([offsetA + 3*i, offsetA + 1 + 3*i]) # y
            bars.append([offsetA + 3*i, offsetA + 3 + 3*i]) # x front
            selectYNegativeBar(self, bars)
            bars.append([offsetA + 3*i + 1, offsetA + 4 + 3*i]) # x back
            selectYPositiveBar(self, bars)
        bars.append([offsetA + 3*(self.nSA//2), offsetA + 1 + 3*(self.nSA//2)])

        # making the pyramids
        bars.append([offsetT-5, offsetA-1])
        selectYNegativeBar(self, bars)
        bars.append([offsetT-3, offsetA-1])
        selectYPositiveBar(self, bars)
        bars.append([offsetA, offsetA-1])
        selectYNegativeBar(self, bars)
        bars.append([offsetA+1, offsetA-1])
        selectYPositiveBar(self, bars)
        for i in range(self.nSA//2):
            bars.append([offsetA + 3*i, offsetA + 2 + 3*i])
            selectYNegativeBar(self, bars)
            bars.append([offsetA + 3*i + 1, offsetA + 2 + 3*i])
            selectYPositiveBar(self, bars)
            bars.append([offsetA + 3*i + 4, offsetA + 2 + 3*i])
            selectYPositiveBar(self, bars)
            bars.append([offsetA + 3*i + 3, offsetA + 2 + 3*i])
            selectYNegativeBar(self, bars)
        self.x_negative_side.append(len(bars)-2)
        self.x_negative_side.append(len(bars)-1)

        # top row
        bars.append([offsetA - 1, offsetT-1]) # connect to top
        selectYPositiveBar(self, bars)
        selectYNegativeBar(self, bars)
        for i in range(self.nSA//2):
            bars.append([offsetA - 1 + i * 3, offsetA + 2 + i * 3])
            selectYNegativeBar(self, bars)
            selectYPositiveBar(self, bars)

        # diagonals
        bars.append([offsetT-3, offsetA])
        for i in range(self.nSA//2):
            bars.append([offsetA+1 + 3*i, offsetA+3 + 3*i])

        # rope
        #bars.append([offsetT+7, -3])

def cur_offset(nodes):
    return len(nodes)

def selectXPositiveBar(crane, bars):
    """
    Select the index from the last bar added in bars array and add to X positive bars
    """
    crane.x_positive_side.append(len(bars) - 1)

def selectXNegativeBar(crane, bars):
    """
    Select the index from the last bar added in bars array and add to X negative bars
    """
    crane.x_negative_side.append(len(bars) - 1)
    
def selectYPositiveBar(crane, bars):
    """
    Select the index from the last bar added in bars array and add to Y positive bars
    """
    crane.y_positive_side.append(len(bars) - 1)

def selectYNegativeBar(crane, bars):
    """
    Select the index from the last bar added in bars array and add to Y negative bars
    """
    crane.y_negative_side.append(len(bars) - 1)



    