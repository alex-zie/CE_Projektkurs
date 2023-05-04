from fem import FEM
from truss import Truss
import numpy as np
import matplotlib.pyplot as plt

def tetrahedron():
    global myTruss

    nodes.append([0,0,-3])
    nodes.append([20,0,-3])
    nodes.append([10, 9, -3])
    nodes.append([10, 6, 3])

    bars.append([0,1])
    bars.append([1,2])
    bars.append([0,2])
    bars.append([0,3])
    bars.append([1,3])
    bars.append([2,3])

    myTruss = Truss(nodes, bars, 0.111, 7850, 1e4)
    myTruss.addSupport(0, 0, 0, 1)
    myTruss.addSupport(1, 0, 1, 0)
    myTruss.addSupport(2, 1, 0, 0)

    myTruss.addExternalForce(3, 0, 0, -1)

def triangle():
    global myTruss

    nodes.append([0, 0, 0])
    nodes.append([0, 1, 0])
    nodes.append([1, 0, 0])

    bars.append([0,1])
    bars.append([1,2])
    bars.append([0,2])

    myTruss = Truss(nodes, bars, 0.111, 7850, 1e4)
    myTruss.addSupport(0, 0, 0, 0)
    #myTruss.addSupport(1, 0, 0, 0)
    #myTruss.addSupport(2, 0, 0, 0)

    myTruss.addExternalForce(2, 0, 0, -10)

def bridge():
    global myTruss

    nodes.append([0, 0, 0])
    nodes.append([1, 0, 0])
    nodes.append([2, 0, 0])
    nodes.append([3, 0, 0])
    nodes.append([4, 0, 0])

    bars.append([0,1])
    bars.append([1,2])
    bars.append([2,3])
    bars.append([3,4])

    myTruss = Truss(nodes, bars, 1, 7850, 1e4)
    myTruss.addSupport(0, 0, 0, 0)
    myTruss.addSupport(4, 0, 0, 0)

    myTruss.addExternalForce(2, 0, 0, -10) # funktioniert nicht

nodes = []
bars = []
myTruss = None

tetrahedron()

fem = FEM(myTruss)

N, R, U = fem.TrussAnalysis()
print(myTruss.F)
print('\nAxial Forces (positive = tension, negative = compression)')
# Anschaulichkeit
print(N[np.newaxis].T)
print('\nReaction Forces (positive = upward, negative = downward')
print(R)
print('\nDeformation at nodes')
print(U)

fem.Plot(np.array(nodes).astype(float), np.array(bars), 'gray', '--', 1, 'Undeformed')
scale = 5
# Berechne die neue Position der Knoten
Dnodes = U * scale + nodes
fem.Plot(np.array(Dnodes), np.array(bars), 'red', '-', 2, 'Deformed')
plt.show()


