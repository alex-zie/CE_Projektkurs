# manually enforced crossections of certain tower and jib beams

class crane_1(Truss):
    """
    Simple crane with cubical tower and jib made of pyramids
    """

    def __init__(self, height, length, ls, A, rho, E, p=True, max_bar_length=-1, optimized_cossections=False):
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

        # indices of bars with thicker crossection
        self.thick_bars = []

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
        if optimized_cossections and (isinstance(A, float) or isinstance(A, int)):
            self.A = A*np.ones(len(bars))
            self.A[self.thick_bars] = 10*A

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
            selectThickBar(self, bars)
            bars.append([4 * i + 1, 4 * i + 5])  # RT
            selectYNegativeBar(self, bars)
            selectXPositiveBar(self, bars)
            selectThickBar(self, bars)
            bars.append([4 * i + 2, 4 * i + 6])  # LB
            selectXNegativeBar(self, bars)
            selectYPositiveBar(self, bars)
            selectThickBar(self, bars)
            bars.append([4 * i + 3, 4 * i + 7])  # RB
            selectXPositiveBar(self, bars)
            selectYPositiveBar(self, bars)
            selectThickBar(self, bars)

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
        selectThickBar(self, bars)
        bars.append([offsetTO - 1, offsetTO - 3])
        selectYPositiveBar(self, bars)
        selectThickBar(self, bars)
        bars.append([offsetTO - 1, offsetTO - 4])
        selectYNegativeBar(self, bars)
        selectThickBar(self, bars)
        bars.append([offsetTO - 1, offsetTO - 5])
        selectYNegativeBar(self, bars)
        selectThickBar(self, bars)
        
        # support line (optional)
        bars.append([offsetTO - 3, offsetTO - 4])
        bars.append([offsetTO - 2, offsetTO - 5])

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
        selectThickBar(self, bars)

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

        # support lines at the bottom (optional) 
        bars.append([offsetT + 1, offsetT -5])  #first line
        for i in range(self.nSGA - 1):
            bars.append([offsetT + 3 * i, offsetT + 4 + 3 * i])

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
            selectThickBar(self, bars)


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

        # support line at bottom
        bars.append([offsetT - 2, offsetTG]) #first line
        for i in range(self.nSA -1):
            bars.append([offsetTG + 1 + i*3, offsetTG + 3 + i*3])

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
            selectThickBar(self, bars)

        bars.append([offsetTO - 2, offsetTO - 3]) # last bar at the end of the crane
        selectXPositiveBar(self, bars)

        bars.append([offsetT - 1, offsetTG + 2]) # top of the Tower with first jib node
        selectYNegativeBar(self, bars)
        selectYPositiveBar(self, bars)
        selectThickBar(self, bars)

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