import numpy as np

class truss:
    """
    This class represents a simple truss that can be passed to the fem class to make truss analyses
    """

    def __init__(self, nodes, bars):
        self.nodes = np.array(nodes).astype(float)
        self.bars = np.array(bars)
        self.F = np.zeros_like(nodes)
        self.DOFCON = np.ones_like(nodes).astype(int)
